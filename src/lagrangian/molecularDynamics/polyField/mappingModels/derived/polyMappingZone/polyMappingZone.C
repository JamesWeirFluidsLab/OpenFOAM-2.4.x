/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "polyMappingZone.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "graph.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(polyMappingZone, 0);

addToRunTimeSelectionTable(polyMappingModel, polyMappingZone, dictionary);



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
polyMappingZone::polyMappingZone
(
    polyMoleculeCloud& molCloud,
    const dictionary& dict
)
:
    polyMappingModel(molCloud, dict),
    propsDict_(dict.subDict(typeName + "Properties")),
//     regionName_(propsDict_.lookup("zoneName")),
//     regionId_(-1),
    translation_(propsDict_.lookup("translationalVector")),
//     oneSpecie_(false),
//     molId_(-1)
//     deleteMols_(false),
    molIds_()
{
//     const cellZoneMesh& cellZones = mesh_.cellZones();
//     regionId_ = cellZones.findZoneID(regionName_);
// 
//     if(regionId_ == -1)
//     {
//         FatalErrorIn("polyMappingZone::polyMappingZone()")
//             << "Cannot find region: " << regionName_ << nl << "in: "
//             << mesh_.time().system()/"molsToDeleteDict"
//             << exit(FatalError);
//     }

    molIds_.clear();

    selectIds ids
    (
        molCloud_.pot(),
        propsDict_
    );

    molIds_ = ids.molIds();

//     deleteMols_ = Switch(propsDict_.lookup("deleteMolsShiftedOutOfZone"));

    findMolsToMap();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

polyMappingZone::~polyMappingZone()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void polyMappingZone::findMolsToMap()
{

    DynamicList<polyMolecule*> molsToDelete;
    
    {
        IDLList<polyMolecule>::iterator molI(molCloud_.begin());

        for
        (
            molI = molCloud_.begin();
            molI != molCloud_.end();
            ++molI
        )
        {
            if(findIndex(molIds_, molI().id()) != -1)
            {
                vector newPos = molI().position() + translation_;

                label cell = mesh_.findCell(newPos);
                
                if(cell == -1)
                {
                    polyMolecule* mol = &molI();
                    molsToDelete.append(mol);                        
                }
                
                molI().position() = newPos;
                molI().cell() = cell;
            }
        }
    }

    //molsToDelete.shrink();

    Info<< " deleting " << molsToDelete.size() << " molecules "
        << endl;

    forAll (molsToDelete, mTD)
    {
        molCloud_.deleteParticle(*(molsToDelete[mTD]));
    }


}


} // End namespace Foam

// ************************************************************************* //
