/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

Class
    Foam::fv::mySolidificationMeltingSource

\*---------------------------------------------------------------------------*/

#ifndef mySolidificationMeltingSource_H
#define mySolidificationMeltingSource_H

#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class mySolidificationMeltingSource Declaration
\*---------------------------------------------------------------------------*/

class mySolidificationMeltingSource
:
    public cellSetOption
{

    // Private Data

        //- Temperature at which melting occurs [K]
        scalar Tmelt_;

        //- Latent heat of fusion [J/kg]
        scalar L_;

        //- Phase fraction under-relaxation coefficient
        scalar relax_;

        //- Name of operand temperature field
        word TName_;

        //- Name of specific heat capacity field
        word CpName_;

        //- Name of operand velocity field
        word UName_;

        //- Name of operand flux field
        word rhoCpPhiName_;

        //- Phase fraction indicator field for VOF
        volScalarField alpha_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        void update();

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn);


public:

    //- Runtime type information
    TypeName("mySolidificationMeltingSource");


    // Constructors

        //- Construct from explicit source name and mesh
        mySolidificationMeltingSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        mySolidificationMeltingSource
        (
            const mySolidificationMeltingSource&
        ) = delete;

        //- No copy assignment
        void operator=(const mySolidificationMeltingSource&) = delete;


    //- Destructor
    ~mySolidificationMeltingSource() = default;


    // Member Functions

        //- Add explicit contribution to enthalpy equation
        virtual void addSup(fvMatrix<scalar>& eqn, const label fieldi);

        //- Add explicit contribution to compressible enthalpy equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Read source dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mySolidificationMeltingSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
