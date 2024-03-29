/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "makeReactionThermo.H"

#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "specie.H"
// #include "perfectGas.H"
// #include "incompressiblePerfectGas.H"
#include "hConstThermo.H"
// #include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "thermo.H"
// #include "rhoConst.H"
// #include "rPolynomial.H"
// #include "perfectFluid.H"
// #include "adiabaticPerfectFluid.H"
// #include "Boussinesq.H"

#include "constTransport.H"
#include "polynomialTransport.H"
// #include "homogeneousMixture.H"
// #include "inhomogeneousMixture.H"
// #include "veryInhomogeneousMixture.H"
// #include "multiComponentMixture.H"
// #include "reactingMixture.H"
// #include "singleStepReactingMixture.H"
#include "singleComponentMixture.H"

#include "myThermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    bPolynomial,
    specie
);

makeReactionThermo
(
    rhoReactionThermo,
    heRhoThermo,
    singleComponentMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hConstThermo,
    bPolynomial,
    specie
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
