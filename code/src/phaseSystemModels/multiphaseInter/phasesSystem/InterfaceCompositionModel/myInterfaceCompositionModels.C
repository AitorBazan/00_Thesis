/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceCompositionModel.H"
#include "InterfaceCompositionModel.H"
//#include "thermoPhysicsTypes.H"
#include "myThermoPhysicsTypes.H"

#include "bPolynomial.H"

#include "pureMixture.H"

#include "rhoThermo.H"

#include "heRhoThermo.H"

#include "solidThermo.H"
#include "heSolidThermo.H"
#include "solidThermoPhysicsTypes.H"

#include "Lee.H"
#include "LeeCNT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInterfacePureType(Type, Thermo, Comp, Mix, Phys, OtherThermo, OtherComp, OtherMix, OtherPhys)\
                                                                               \
    typedef Thermo<Comp, Mix<Phys>>                                            \
        Type##Thermo##Comp##Mix##Phys;                                         \
                                                                               \
    typedef OtherThermo<OtherComp, OtherMix<OtherPhys>>                        \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys;              \
                                                                               \
    addInterfaceCompositionToRunTimeSelectionTable                             \
    (                                                                          \
        Type,                                                                  \
        Type##Thermo##Comp##Mix##Phys,                                         \
        Type##Other##OtherThermo##OtherComp##OtherMix##OtherPhys               \
    )

// Addition to the run-time selection table
#define addInterfaceCompositionToRunTimeSelectionTable(Type, Thermo, OtherThermo)\
                                                                               \
    typedef Type<Thermo, OtherThermo>                                          \
        Type##Thermo##OtherThermo;                                             \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Type##Thermo##OtherThermo,                                             \
        (                                                                      \
            word(Type##Thermo##OtherThermo::typeName_()) + "<"                 \
          + word(Thermo::typeName) + ","                                       \
          + word(OtherThermo::typeName) + ">"                                  \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        interfaceCompositionModel,                                             \
        Type##Thermo##OtherThermo,                                             \
        dictionary                                                             \
    )
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    using namespace meltingEvaporationModels;

    //NOTE: First thermo (from) and second otherThermo (to)


    // Lee model definitions

        // From pure phase (poly) to phase (solidThermo)
        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constbPolFluidHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );
        
        makeInterfacePureType
        (
            LeeCNT,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            constbPolFluidHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );

        makeInterfacePureType
        (
            Lee,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            polybPolFluidHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );
        
        makeInterfacePureType
        (
            LeeCNT,
            heRhoThermo,
            rhoThermo,
            pureMixture,
            polybPolFluidHThermoPhysics,
            heSolidThermo,
            solidThermo,
            pureMixture,
            hConstSolidThermoPhysics
        );        
        // interfaceHeatResistance model definitions

        // From pure phase (poly) to phase (solidThermo)
        // makeInterfacePureType
        // (
        //     interfaceHeatResistance,
        //     heRhoThermo,
        //     rhoThermo,
        //     pureMixture,
        //     constbPolFluidHThermoPhysics,
        //     heSolidThermo,
        //     solidThermo,
        //     pureMixture,
        //     hConstSolidThermoPhysics
        // );

}


// ************************************************************************* //
