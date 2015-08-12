/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "polyMoleculeCloud.H"
#include "polyAllConfigurations.H"
#include "fvMesh.H"
#include "polyMolsToDelete.H"
#include "polyMappingModels.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<polyMolecule>, 0);
};

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMoleculeCloud::buildConstProps()
{
    Info<< nl << "Reading moleculeProperties dictionary." << endl;

    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size()); 

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    dictionary moleculeProperties
    (
        moleculePropertiesDict.subDict("moleculeProperties")
    );

    forAll(idList, i)
    {
        const word& id(idList[i]);

        const dictionary& molDict(moleculeProperties.subDict(id));

        const word cloudType = molDict.lookup("cloudType");

        if(cloudType == "polyMoleculeCloud")
        {
            List<word> siteIdNames = molDict.lookup("siteIds");
    
            List<label> siteIds(siteIdNames.size());
    
            forAll(siteIdNames, sI)
            {
                const word& siteId = siteIdNames[sI];
    
                siteIds[sI] = findIndex(siteIdList, siteId);
    
                if (siteIds[sI] == -1)
                {
                    FatalErrorIn("polyMoleculeCloud.C") << nl
                        << siteId << " site not found."
                        << nl << abort(FatalError);
                }
            }
    
            polyMolecule::constantProperties& constProp = constPropList_[i];
    
            constProp = polyMolecule::constantProperties(molDict, redUnits_, siteIds);
        }
    }
}


void Foam::polyMoleculeCloud::setSiteSizesAndPositions()
{
    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        const polyMolecule::constantProperties& cP = constProps(mol().id());

        mol().setSiteSizes(cP.nSites());

        mol().setSitePositions(cP);
    }
}

void Foam::polyMoleculeCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    iterator mol(this->begin());

    for
    (
        mol = this->begin();
        mol != this->end();
        ++mol
    )
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }
    
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }
}

void Foam::polyMoleculeCloud::removeHighEnergyOverlaps()
{
    Info<< nl << "Removing high energy overlaps, limit = "
        << pot_.potentialEnergyLimit()
        << nl << "Removal order:";

    forAll(pot_.removalOrder(), rO)
    {
        if(pot_.removalOrder()[rO] != -1)
        {
            Info<< ' ' << pot_.idList()[pot_.removalOrder()[rO]];
        }
    }

    Info<< nl ;

    label initialSize = this->size();

    if (Pstream::parRun())
    {
        reduce(initialSize, sumOp<label>());
    }

    buildCellOccupancy();
    
    label nMolsDeleted = 0;

    prepareInteractions();

    setIPL();
    
    iL_.setRIPL();
    
    label nMolsInt = 0;
    label nMolsExt = 0;
    label nMolsRef = 0;

    {
        DynamicList<polyMolecule*> molsToDelete;

        polyMolecule* molI = NULL;
        polyMolecule* molJ = NULL;

        forAll(ipl_, c)
        {
            nMolsExt = ipl_[c].size();
            nMolsInt = cellOccupancy_[c].size();
            nMolsRef = iL_.ripl()[c].size();
                 
            for (int i = 0; i < nMolsInt; i++)
            {
                molI = cellOccupancy_[c][i];

				label idI = molI->id();
				bool molIDeleted = false;

				for (int j = 0; j < nMolsInt; j++)
				{
					if(j > i)
					{
						molJ = cellOccupancy_[c][j];

						label molJDeleted = findIndex(molsToDelete, molJ);

						if(!molIDeleted && (molJDeleted == -1))
						{
							if
							(
								evaluatePotentialLimit
								(
									molI,
									molJ,
									pot_.potentialEnergyLimit()
								)
							)
							{
								label idJ = molJ->id();

								label removeIdI = findIndex(pot_.removalOrder(), idI);
								label removeIdJ = findIndex(pot_.removalOrder(), idJ);

								if(removeIdI < removeIdJ)
								{
									molsToDelete.append(molI);
									molIDeleted = true;
								}
								else
								{
									molsToDelete.append(molJ);
								}
							}
						}
					}
				}

				for (int j = 0; j < nMolsExt; j++)
				{
					molJ = ipl_[c][j];

					label molJDeleted = findIndex(molsToDelete, molJ);

					if(!molIDeleted && (molJDeleted == -1))
					{
						if
						(
							evaluatePotentialLimit
							(
								molI,
								molJ,
								pot_.potentialEnergyLimit()
							)
						)
						{
							label idJ = molJ->id();

							label removeIdI = findIndex(pot_.removalOrder(), idI);
							label removeIdJ = findIndex(pot_.removalOrder(), idJ);

							if(removeIdI < removeIdJ)
							{
								molsToDelete.append(molI);
								molIDeleted = true;
							}
							else
							{
								molsToDelete.append(molJ);
							}
						}
					}
				}

				for (int j = 0; j < nMolsRef; j++)
				{
					molJ = iL_.ripl()[c][j];

					if(!molIDeleted)
					{
						if
						(
							evaluatePotentialLimit
							(
								molI,
								molJ,
								pot_.potentialEnergyLimit()
							)
						)
						{
							label idJ = molJ->id();

							label removeIdI = findIndex(pot_.removalOrder(), idI);
							label removeIdJ = findIndex(pot_.removalOrder(), idJ);

							if(removeIdI < removeIdJ)
							{
								molsToDelete.append(molI);
								molIDeleted = true;
							}
							else if
							(
								(removeIdI == removeIdJ)
							)
							{
								if (molI->trackingNumber() > molJ->trackingNumber())
								{
									molsToDelete.append(molI);
									molIDeleted = true;
								}
							}
						}
					}
				}
            }
        }

        forAll (molsToDelete, mTD)
        {
            nMolsDeleted++;
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    label newSize = this->size();

    if (Pstream::parRun())
    {
        reduce(newSize, sumOp<label>());
    }
    
    if (Pstream::parRun())
    {
        reduce(nMolsDeleted, sumOp<label>());
    }    

    if(nMolsDeleted > 0)
    {
        // to make sure the user sees this
        for (int j = 0; j < 50; j++)
        {
            Info << nl << "WARNING: molecules removed due to overlaps = "
                << nMolsDeleted <<  endl;
        }
    }
    
}


Foam::label Foam::polyMoleculeCloud::nSites() const
{
    label n = 0;

    const_iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        n += constProps(mol().id()).nSites();
    }

    return n;
}


void Foam::polyMoleculeCloud::checkMoleculesInMesh()
{
    Info << nl << "checking cell-molecule addressing" << endl;

    DynamicList<polyMolecule*> molsToDelete;

    label initialSize = this->size();

    iterator mol(this->begin());

    label noOfModifiedMols = 0;

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        
        label cell = -1;
        label tetFace = -1;
        label tetPt = -1;

        mesh_.findCellFacePt
        (
            mol().position(),
            cell,
            tetFace,
            tetPt
        );
        
        if(cell != -1)
        {
            if(mol().cell() != cell)
            {
                mol().cell() = cell;
                mol().tetFace() = tetFace;
                mol().tetPt() = tetPt;
                noOfModifiedMols++;
            }
        }
        else
        {
            Pout<< "WARNING - molecule outside mesh at position = "
                << mol().position()
                << endl;
                
            polyMolecule* molI = &mol();
            molsToDelete.append(molI);
        }
    }

    if(noOfModifiedMols > 0)
    {
        Pout<< tab << " molecules that changed cell = " 
            << noOfModifiedMols
            << endl;
    }

    forAll (molsToDelete, mTD)
    {
        deleteParticle(*(molsToDelete[mTD]));
    }

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info<< tab <<" molecules removed from outside mesh = " 
        << molsRemoved 
        << endl;
}

// destructor

Foam::polyMoleculeCloud::~polyMoleculeCloud()
{
     if(staticDil_!=NULL){
	for(int i = 0; i < meshSize_; i++){
             delete [] staticDil_[i];
	}
	delete [] staticDil_;
     }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



//- Use for running MD (mdFoam)
Foam::polyMoleculeCloud::polyMoleculeCloud
(
    Time& t,
    const polyMesh& mesh,
    const potential& pot,
    const reducedUnits& rU,
    cachedRandomMD& rndGen
)
:
    Cloud<polyMolecule>(mesh, "polyMoleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    redUnits_(rU),
    rndGen_(rndGen),
    cellOccupancy_(mesh_.nCells()),
    constPropList_(),
    fields_(t, mesh_, *this),
    boundaries_(t, mesh, *this),
    controllers_(t, mesh, *this),
    trackingInfo_(mesh, *this),
    moleculeTracking_(),
    cyclics_(t, mesh_, -1), 
    iL_(mesh, rU, cyclics_, pot_.pairPotentials().rCutMax(), "poly"),
    ipl_(mesh.nCells()),
    staticDil_(0),
    dilSizeList_(mesh_.nCells()),
    meshSize_(mesh_.nCells())
{
    polyMolecule::readFields(*this);

    rndGen.initialise(this->size() != 0 ? this->size() : 10000); //Initialise the random number cache (initialise to 10000 if size is zero)

    buildConstProps();
    
    // copy all dil information
    staticDil_ = new int*[mesh_.nCells()];
    forAll(iL_.dil(),d){
	int size = iL_.dil()[d].size();
	dilSizeList_[d] = size;
	staticDil_[d] = new int[size];
	forAll(iL_.dil()[d],j){
	    staticDil_[d][j] = iL_.dil()[d][j];
	}
    }

    setSiteSizesAndPositions();

    checkMoleculesInMesh();

    // set tracking numbers
    setTrackingNumbers();

    removeHighEnergyOverlaps();
    
    buildCellOccupancy();
    clearLagrangianFields();
    calculateForce();
    updateAcceleration();
    
    fields_.createFields();
    boundaries_.setInitialConfig();
    controllers_.initialConfig();
    
    // TESTS
    writeReferredCloud();
}



//- general constructor
Foam::polyMoleculeCloud::polyMoleculeCloud
(
    Time& t,
    const polyMesh& mesh,
    const potential& pot,
    const reducedUnits& rU,
    cachedRandomMD& rndGen, 
    const word& option,
    const bool& clearFields
)
    :
    Cloud<polyMolecule>(mesh, "polyMoleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    redUnits_(rU),
    rndGen_(rndGen),    
    cellOccupancy_(mesh_.nCells()),
    constPropList_(),
    fields_(t, mesh_),
    boundaries_(t, mesh),
    controllers_(t, mesh),
    trackingInfo_(mesh, *this),
    moleculeTracking_(),
    cyclics_(t, mesh_, -1),
    iL_(mesh, rU, cyclics_, pot_.pairPotentials().rCutMax(), "poly"),
    ipl_(mesh.nCells())
{
    polyMolecule::readFields(*this);

    label initialMolecules = this->size();

    rndGen.initialise(initialMolecules != 0 ? initialMolecules : 10000); //Initialise the random number cache (initialise to 10000 if size is zero)

    if (Pstream::parRun())
    {
        reduce(initialMolecules, sumOp<label>());
    }
   
    if(clearFields)
    {
        Info << "clearing existing field of molecules " << endl;

        clear();

        initialMolecules = 0;
    }
    
    buildConstProps();
    setSiteSizesAndPositions();

    if((option == "mdInitialise") && clearFields)
    {
        polyAllConfigurations conf(mesh, *this);
        conf.setInitialConfig();
        buildCellOccupancy();
    }
    else if((option == "mdInitialise") && !clearFields)
    {
        checkMoleculesInMesh();
        buildCellOccupancy();
        polyAllConfigurations conf(mesh, *this);
        conf.setInitialConfig();
    }
    else if(option == "delete")
    {
        checkMoleculesInMesh();
        setTrackingNumbers();
        buildCellOccupancy();
        prepareInteractions();
        polyMolsToDelete molsDel(mesh_, *this);
    }
    else if(option == "mapping")
    {
        polyMappingModels molsToMap(mesh_, *this);
        buildCellOccupancy();
    }
    else if(option == "quickMapping")
    {
        checkMoleculesInMesh();
        buildCellOccupancy();
    }
    else if(option == "NULL")
    {
        buildCellOccupancy();
    }
    else 
    {
        Info << "ERROR" << endl;
    }

    label finalMolecules = this->size();
    
    if (Pstream::parRun())
    {
        reduce(finalMolecules, sumOp<label>());
    }

    Info << nl << "Initial molecules = " << initialMolecules 
         << ", modified molecules = " << finalMolecules - initialMolecules
         << ", total molecules: " << finalMolecules 
         << endl;
}

// * * * * * * * * * * * * * * * * Static Constructors  * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::polyMoleculeCloud> Foam::polyMoleculeCloud::New
(
    Time& t,
    const polyMesh& mesh,
    const potential& pot,
    const reducedUnits& rU,
    cachedRandomMD& rndGen
)
{
    return autoPtr<polyMoleculeCloud>
    (
        new polyMoleculeCloud(t, mesh, pot, rU, rndGen)
    );
}

Foam::autoPtr<Foam::polyMoleculeCloud> Foam::polyMoleculeCloud::New
(
    Time& t,
    const polyMesh& mesh,
    const potential& pot,
    const reducedUnits& rU,
    cachedRandomMD& rndGen,
    const word& option,
    const bool& clearFields
)
{
    return autoPtr<polyMoleculeCloud>
    (
        new polyMoleculeCloud(t, mesh, pot, rU, rndGen, option, clearFields)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void  Foam::polyMoleculeCloud::createMolecule
(
    const vector& position,
    const label cell,
    const label tetFace,
    const label tetPt,
    const tensor& Q,
    const vector& v,
    const vector& a,
    const vector& pi,
    const vector& tau,
    const vector& specialPosition,
    const label special,
    const label id,
    const scalar& fraction,
    const label trackingNumber
)
{
    addParticle
    (
        new polyMolecule
        (
            mesh_,
            position,
            cell,
            tetFace,
            tetPt,
            Q,
            v,
            a,
            pi,
            tau,
            specialPosition,
            constProps(id),
            special,
            id,
            fraction,
            trackingNumber
        )
    );
}

void Foam::polyMoleculeCloud::preliminaries()
{
    // control
    controllers_.controlVelocitiesI();
}

// update initial half velocity
void Foam::polyMoleculeCloud::initialHalfVelocity()
{
    velocityUpdate(mesh_.time().deltaT().value());
}

// move molecules (tracking)
void Foam::polyMoleculeCloud::move()
{
    polyMolecule::trackingData td1(*this, 1);
    Cloud<polyMolecule>::move(td1, mesh_.time().deltaTValue());

    updateAfterMove(mesh_.time().deltaT().value());
}

// control
void Foam::polyMoleculeCloud::controlBeforeForces()
{
    controllers_.controlPriorToForces();
}

// update acceleration from net forces
void Foam::polyMoleculeCloud::updateAcceleration()
{
    accelerationUpdate();
}

// control
void Foam::polyMoleculeCloud::controlAfterForces()
{
    boundaries_.controlAfterForces();
    controllers_.controlState();
}

// update final half velocity
void Foam::polyMoleculeCloud::finalHalfVelocity()
{
    velocityUpdate(mesh_.time().deltaT().value());
}

void Foam::polyMoleculeCloud::postPreliminaries()
{
    // control
    controllers_.controlVelocitiesII();

    fields_.calculateFields();
    fields_.writeFields();

    boundaries_.calculateProps();
    boundaries_.outputResults();

    controllers_.calculateStateProps();
    controllers_.outputStateResults();

    if(mesh_.time().outputTime())
    {
        writeReferredCloud();
    }

    trackingInfo_.clean(); 
}

void Foam::polyMoleculeCloud::evolve()
{
    evolveBeforeForces();
    calculateForce();
    evolveAfterForces();
}

void Foam::polyMoleculeCloud::evolveBeforeForces()
{
    preliminaries();
    initialHalfVelocity();
    move();
    boundaries_.controlAfterMove();
    buildCellOccupancy();
    controlBeforeForces();
    clearLagrangianFields();
}

void Foam::polyMoleculeCloud::evolveAfterForces()
{
    updateAcceleration();
    controlAfterForces();
    finalHalfVelocity();
    postPreliminaries();
}

void Foam::polyMoleculeCloud::updateAfterMove(const scalar& trackTime)
{
    forAllIter(polyMoleculeCloud, *this, mol)
    {
        if(!mol().frozen())
        {
            const polyMolecule::constantProperties& cP = constProps(mol().id());
            mol().updateAfterMove(cP, trackTime);
        }
    }
}

void Foam::polyMoleculeCloud::velocityUpdate(const scalar& trackTime)
{
    forAllIter(polyMoleculeCloud, *this, mol)
    {
        if(!mol().frozen())
        {
            const polyMolecule::constantProperties& cP = constProps(mol().id());
            mol().updateHalfVelocity(cP, trackTime);
        }
    }
}

void Foam::polyMoleculeCloud::accelerationUpdate()
{
    forAllIter(polyMoleculeCloud, *this, mol)
    {
        if(!mol().frozen())
        {
            const polyMolecule::constantProperties& cP = constProps(mol().id());
            mol().updateAcceleration(cP);
        }
    }
}

//- used if you want to read a new field at every time-step from an input file
//- e.g. to be used in a utility that computes measurements
// Used by reconstructXmol utility - reconstructPar does not produce the XMOL
// files after parallel processing
void Foam::polyMoleculeCloud::readNewField()
{
    label initialSize = this->size();

    clear();

    IOPosition<Cloud<polyMolecule> > ioP(*this);

    if (ioP.headerOk())
    {
        ioP.readData(*this, false);
        ioP.close();
    }
    else
    {
        // WARNING
        WarningIn("readNewField()")
            << "Cannot read particle positions file " << nl
            << "    " << ioP.objectPath() << nl
            << "    assuming the initial cloud contains 0 particles." << endl;        
    }    
    
    particle::readFields(*this);

    polyMolecule::readFields(*this);

    if (this->size() != initialSize)
    {
        Info << "Changed polyMoleculeCloud size, from: " 
                << initialSize << ", to: " << this->size() << endl;
    }

    setSiteSizesAndPositions();    
}

void Foam::polyMoleculeCloud::clearLagrangianFields()
{
    iterator mol(this->begin());

    // Set accumulated quantities to zero
    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().a() = vector::zero;

        mol().tau() = vector::zero;

        mol().siteForces() = vector::zero;

        mol().potentialEnergy() = 0.0;

        mol().rf() = tensor::zero;

        mol().R() = GREAT;
    }
}

void Foam::polyMoleculeCloud::calculateForce()
{
    ompCalculatePairForces();
}

void Foam::polyMoleculeCloud::setIPL()
{
   forAll(ipl_, c)
   {
        ipl_[c].clear();
   }

    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        forAll(iL_.inverseDIL()[mol().cell()], c)
        {
            ipl_[iL_.inverseDIL()[mol().cell()][c]].append(&mol());
        }
    }
}

void Foam::polyMoleculeCloud::rebuildCellOccupancy()
{
    buildCellOccupancy();
}

void Foam::polyMoleculeCloud::prepareInteractions()
{
    iL_.setReferredParticles(cellOccupancy());
}

void Foam::polyMoleculeCloud::evaluatePairCritical
(
    polyMolecule* molI,
    polyMolecule* molJ    
)
{
    const pairPotentialList& pairPot = pot_.pairPotentials();
    const pairPotential& electrostatic = pairPot.electrostatic();

    label idI = molI->id();
    label idJ = molJ->id();

    const polyMolecule::constantProperties& constPropI(constProps(idI));
    const polyMolecule::constantProperties& constPropJ(constProps(idJ));

    // pair potential interactions
    controllers_.controlDuringForceComputation(molI, molJ); 

    fields_.measurementsDuringForceComputation
    (
        molI,
        molJ
    );

    if(!molI->frozen() || !molJ->frozen())
    {
        // fraction
        scalar f = molI->fraction();
    
        if(molJ->fraction() < f)
        {
            f = molJ->fraction();
        }

        forAll(constPropI.pairPotSites(), pI)
        {
        	label sI = constPropI.pairPotSites()[pI];

            forAll(constPropJ.pairPotSites(), pJ)
            {
                label sJ = constPropJ.pairPotSites()[pJ];

                Foam::vector rsIsJ = molI->sitePositions()[sI] - molJ->sitePositions()[sJ];

                scalar rsIsJMagSq = magSqr(rsIsJ);

                label idsI = constPropI.sites()[sI].siteId();
                label idsJ = constPropJ.sites()[sJ].siteId();

                if(pairPot.rCutSqr(idsI, idsJ, rsIsJMagSq))
                {
                    scalar rsIsJMag = mag(rsIsJ);

                    Foam::vector fsIsJ = (rsIsJMag == 0.0? Foam::vector::zero : (f * (rsIsJ/rsIsJMag) * pairPot.force(idsI, idsJ, rsIsJMag)));


                    scalar potentialEnergy
                    (
                        f*pairPot.energy(idsI, idsJ, rsIsJMag)
                    );
                    Foam::vector rIJ = molI->position() - molJ->position();
                    Foam::tensor virialContribution = (rsIsJMagSq == 0.0? (rsIsJ*fsIsJ)*(rsIsJ & rIJ): (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq);

		    molI->siteForces()[sI] += fsIsJ;
                    molI->potentialEnergy() += 0.5*potentialEnergy;
                    molI->rf() += virialContribution;
#pragma omp critical (forces)
		    {
		    //molI->siteForces()[sI] += fsIsJ;

                    molJ->siteForces()[sJ] += -fsIsJ;
		    
                    //molI->potentialEnergy() += 0.5*potentialEnergy;
        
                    molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                    //Foam::vector rIJ = molI->position() - molJ->position();
        
                    //Foam::tensor virialContribution = (rsIsJMagSq == 0.0? (rsIsJ*fsIsJ)*(rsIsJ & rIJ): (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq);
        
                   // molI->rf() += virialContribution;
                    molJ->rf() += virialContribution;
		    }
                    fields_.measurementsDuringForceComputationSite
                    (
                        molI,
                        molJ,
                        sI,
                        sJ
                    );
                }
            }
        }

        {
            Foam::vector rIJ = molI->position() - molJ->position();
    
            scalar rIJMag = mag(rIJ);
    
            if(molI->R() > rIJMag)
            {
                molI->R() = rIJMag;
            }
    
            if(molJ->R() > rIJMag)
            {
                molJ->R() = rIJMag;
            }
        }

        forAll(constPropI.electrostaticSites(), pI)
        {
            label sI = constPropI.electrostaticSites()[pI];
    
            forAll(constPropJ.electrostaticSites(), pJ)
            {
                label sJ = constPropJ.electrostaticSites()[pJ];
    
                vector rsIsJ =
                molI->sitePositions()[sI] - molJ->sitePositions()[sJ];
    
                scalar rsIsJMagSq = magSqr(rsIsJ);
        
                if(rsIsJMagSq < electrostatic.rCutSqr())
                {
                    scalar rsIsJMag = mag(rsIsJ);
        
                    scalar chargeI = constPropI.sites()[sI].siteCharge();
                    scalar chargeJ = constPropJ.sites()[sJ].siteCharge();
        
                    vector fsIsJ =
                        f*(rsIsJ/rsIsJMag)
                        *chargeI*chargeJ*electrostatic.force(rsIsJMag);
        
                    //molI->siteForces()[sI] += fsIsJ;
                    //molJ->siteForces()[sJ] += -fsIsJ;
        
                    scalar potentialEnergy =
                        f*chargeI*chargeJ
                        *electrostatic.energy(rsIsJMag);
                    
		    vector rIJ = molI->position() - molJ->position();
                    tensor virialContribution =
                        (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
        
                    molI->siteForces()[sI] += fsIsJ;
                    molI->potentialEnergy() += 0.5*potentialEnergy;
                    molI->rf() += virialContribution;
                   
#pragma omp critical (force2)
{ 
		    molJ->siteForces()[sJ] += -fsIsJ;
        
                    molJ->potentialEnergy() += 0.5*potentialEnergy;
        
                    //vector rIJ = molI->position() - molJ->position();
        
                    //tensor virialContribution =
                    //    (rsIsJ*fsIsJ)*(rsIsJ & rIJ)/rsIsJMagSq;
        
                    molJ->rf() += virialContribution;
}    
                    fields_.measurementsDuringForceComputationSite
                    (
                        molI,
                        molJ,
                        sI,
                        sJ
                    );
                }
            }
        }
    }  
}

void Foam::polyMoleculeCloud::ompCalculatePairForces()
{
   prepareInteractions();
   
   std::vector<std::vector<polyMolecule*> > tempCellOcc(mesh_.nCells()); 
   forAll(cellOccupancy_,c){
       int size = cellOccupancy_[c].size();
       tempCellOcc[c].resize(size);
       forAll(cellOccupancy_[c],j){
	   tempCellOcc[c][j] = cellOccupancy_[c][j]; 
       }
   }
   
   const labelListList& dil = iL_.dil();
   int** tempdil = &staticDil_[0];

   polyMolecule* cellOcc = tempCellOcc[0][0];
   //pointer to dilsizelist and celloccupancysizelist
   int* dilsize = &dilSizeList_[0];
   int sizedil = dilSizeList_.size();

#pragma omp parallel 
{
#pragma omp single nowait
{
   for(int d = 0; d < sizedil; d++)
#pragma omp task firstprivate(tempCellOcc)
   {
	std::vector<polyMolecule*> templist = tempCellOcc[d];
	int sizecellocc = templist.size();
	for(int i = 0; i < sizecellocc; i++)
	{
	    polyMolecule* molI = templist[i]; 
	    for(int j = i+1; j < sizecellocc; j++)
	    {
		polyMolecule* molJ = templist[j];
		evaluatePair(molI,molJ);
	    }	    
	}
   }
}//end single region
}//end parallel region 

    for(int d = 0; d < sizedil; d++)
    {
#pragma omp parallel
      {
#pragma omp single nowait
	{
	std::vector<polyMolecule*> templistI = tempCellOcc[d];
	int size = templistI.size();
	for(int cellI = 0; cellI < size; cellI++)
#pragma omp task firstprivate(templistI,tempCellOcc,dilsize)
	{
	    polyMolecule* molI = templistI[cellI];
	    int sizeDilD = dilsize[d];
	    for(int dj = 0; dj < sizeDilD; dj++)
		{
		    std::vector<polyMolecule*> cellJList = tempCellOcc[tempdil[d][dj]];
		    int sizeCellJ = cellJList.size();
		    for(int cellJ = 0; cellJ < sizeCellJ; cellJ++)
		    {
			polyMolecule* molJ = cellJList[cellJ];
			evaluatePairCritical(molI,molJ);
		    }//cellJ ends
		}// dj ends
	}
	}//end single
      }//end parallel region
    }
   
   // Real-Referred interactions
   for(int r = 0; r < iL_.refCellsParticles().size(); r++)
   {
       const List<label>& realCells = iL_.refCells()[r].neighbouringCells();
       for(int i = 0; i < iL_.refCellsParticles()[r].size(); i++)
       {
	   polyMolecule* molJ = iL_.refCellsParticles()[r][i];
	   for(int rc = 0; rc < realCells.size(); rc++)
	   {
	       List<polyMolecule*> cellIlist = cellOccupancy_[realCells[rc]];
	       for(int i = 0; i < cellIlist.size(); i++)
	       {
		    polyMolecule* molI = cellIlist[i];
		    evaluatePair(molI,molJ);  
	       }
	   }
       }
   }

}// ompCalculatePairForces
void Foam::polyMoleculeCloud::calculatePairForces()
{

    prepareInteractions();

    polyMolecule* molI = NULL;
    polyMolecule* molJ = NULL;

    {
        // Real-Real interactions
        const labelListList& dil = iL_.dil();

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

				forAll(dil[d], interactingCells)
				{
					List<polyMolecule*> cellJ =
						cellOccupancy_[dil[d][interactingCells]];

					forAll(cellJ, cellJMols)
					{
						molJ = cellJ[cellJMols];

						evaluatePair(molI, molJ);
					}
				}

				forAll(cellOccupancy_[d], cellIOtherMols)
				{
					molJ = cellOccupancy_[d][cellIOtherMols];

					if (molJ > molI)
					{
						evaluatePair(molI, molJ);
					}
				}
            }
        }
    }
    {
        // Real-Referred interactions
        forAll(iL_.refCellsParticles(), r)
        {
            const List<label>& realCells = iL_.refCells()[r].neighbouringCells();

            forAll(iL_.refCellsParticles()[r], i)
            {
            	molJ = iL_.refCellsParticles()[r][i];

				forAll(realCells, rC)
				{
					List<polyMolecule*> molsInCell = cellOccupancy_[realCells[rC]];

					forAll(molsInCell, j)
					{
						molI = molsInCell[j];
						evaluatePair(molI, molJ);
					}
				}
            }
        }
    }

}

void Foam::polyMoleculeCloud::writeXYZ(const fileName& fName) const
{
    OFstream os(fName);

    os << nSites() << nl << "polyMoleculeCloud site positions in angstroms" << nl;

    const_iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        const polyMolecule::constantProperties& cP = constProps(mol().id());

        forAll(mol().sitePositions(), i)
        {
            const point& sP = mol().sitePositions()[i];

            os << pot_.siteIdList()[cP.sites()[i].siteId()]
                << ' ' << sP.x()*redUnits_.refLength()*1.0e10
                << ' ' << sP.y()*redUnits_.refLength()*1.0e10
                << ' ' << sP.z()*redUnits_.refLength()*1.0e10
                << nl;
        }
    }
}

// new function added to write the referred cloud 
// - to visualise the particles in ParaFOAM/VMD
void Foam::polyMoleculeCloud::writeReferredCloud()
{
    if(iL_.write())
    {
        const Time& runTime = mesh_.time();
        
        Info << "Writing out referred cloud" << endl;

        fileName timePath(runTime.path()/runTime.timeName()/"lagrangian");
        
        if (!isDir(timePath))
        {
            mkDir(timePath);
        }
        
        fileName fName1(timePath/"referredCloud.xmol"); // VMD          
        fileName fName2(timePath/"referredCloud_RU.xmol");  //ParaFOAM
        
        label nParticles = iL_.referredCloud().size();
        
        Info << "number of particles = " << nParticles << endl;

        label nSites = 0;
        
        forAllIter
        (
            IDLList<polyMolecule>,
            iL_.referredCloud(),
            mol
        )
        {
            nSites += mol().sitePositions().size();
        }

        Info << "number of sites = " << nSites << endl;
        
        OFstream os1(fName1);
        OFstream os2(fName2);
        
        os1 << nSites << nl << "referred polyMoleculeCloud site positions in angstroms" << nl;
        os2 << nSites << nl << "referred polyMoleculeCloud site positions in reduced units" << nl;    
        
        forAllIter
        (
            IDLList<polyMolecule>,
            iL_.referredCloud(),
            mol
        )
        {
            const polyMolecule::constantProperties& cP = constProps(mol().id());

            forAll(mol().sitePositions(), j)
            {            
            	const point& sP = mol().sitePositions()[j];

                os1 << pot_.siteIdList()[cP.sites()[j].siteId()]
                        << ' ' << sP.x()*redUnits_.refLength()*1e10
                        << ' ' << sP.y()*redUnits_.refLength()*1e10
                        << ' ' << sP.z()*redUnits_.refLength()*1e10
                        << nl;
                        
                os2 << pot_.siteIdList()[cP.sites()[j].siteId()]
                        << ' ' << sP.x()
                        << ' ' << sP.y()
                        << ' ' << sP.z()
                        << nl;
            }
        }
    }
}

void Foam::polyMoleculeCloud::testTrackingNumbers()
{
    moleculeTracking_.resetTrackingNumbers();

    if(moleculeTracking_.resetTracking())
    {
        setTrackingNumbers();  
    }
}

void Foam::polyMoleculeCloud::setTrackingNumbers()
{
    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().trackingNumber() = getTrackingNumber();
    }
}

void Foam::polyMoleculeCloud::insertMolInCellOccupancy(polyMolecule* mol)
{
    cellOccupancy_[mol->cell()].append(mol);
}

void Foam::polyMoleculeCloud::removeMolFromCellOccupancy
(
    polyMolecule* molI
)
{
    DynamicList<polyMolecule*> updatedMolsInCell(0);

    const label& cellI = molI->cell();

    {
        const List<polyMolecule*>& molsInCell = cellOccupancy_[cellI];
    
        forAll(molsInCell, m)
        {
            polyMolecule* molJ = molsInCell[m];
    
            if(molI != molJ)
            {
                updatedMolsInCell.append(molJ);
            }
        }
    }

    cellOccupancy_[cellI].clear();
    cellOccupancy_[cellI].transfer(updatedMolsInCell);
}


void Foam::polyMoleculeCloud::removeMolFromCellOccupancy
(
    const label& cellMolId,
    const label& cell
)
{
    DynamicList<polyMolecule*> molsInCell(0);

    forAll(cellOccupancy_[cell], c)
    {
        if(c != cellMolId)
        {
            molsInCell.append(cellOccupancy_[cell][c]);
        }
    }

    cellOccupancy_[cell].clear();
    cellOccupancy_[cell].transfer(molsInCell);
}

Foam::label Foam::polyMoleculeCloud::getTrackingNumber()
{
    return moleculeTracking_.getTrackingNumber();
}

// ************************************************************************* //
