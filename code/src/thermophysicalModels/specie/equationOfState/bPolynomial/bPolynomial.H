/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::bPolynomial

Group
    grpSpecieEquationOfState

Description
    Incompressible, polynomial form of equation of state, using a polynomial
    function for density.

Usage
    \table
        Property     | Description
        rhoCoeffs<8> | Density polynomial coefficients
    \endtable

    Example of the specification of the equation of state:
    \verbatim
    equationOfState
    {
        rhoCoeffs<8>    ( 1000 -0.05 0.003 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expression is evaluated as so:

        \f[
            \rho = 1000 - 0.05 T + 0.003 T^2
        \f]

Note
    Input in [kg/m3], but internally uses [kg/m3/kmol].

SourceFiles
    bPolynomialI.H
    bPolynomial.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef bPolynomial_H
#define bPolynomial_H

#include "autoPtr.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Specie, int PolySize>
class bPolynomial;

template<class Specie, int PolySize>
bPolynomial<Specie, PolySize> operator+
(
    const bPolynomial<Specie, PolySize>&,
    const bPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
bPolynomial<Specie, PolySize> operator*
(
    const scalar,
    const bPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
bPolynomial<Specie, PolySize> operator==
(
    const bPolynomial<Specie, PolySize>&,
    const bPolynomial<Specie, PolySize>&
);

template<class Specie, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const bPolynomial<Specie, PolySize>&
);


/*---------------------------------------------------------------------------*\
                        Class bPolynomial Declaration
\*---------------------------------------------------------------------------*/

template<class Specie, int PolySize=8>
class bPolynomial
:
    public Specie
{
    // Private Data

        //- Reference density
        scalar rhoRef_;

        //- Reference temperature
        scalar Tref_;

        //- Density polynomial coefficients
        Polynomial<PolySize> rhoCoeffs_;

    // Private Member Functions

        //- Coeffs name. Eg, "rhoCoeffs<10>"
        inline static word coeffsName(const char* name)
        {
            return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        }


public:

    // Generated Methods: copy construct, copy assignment

    // Constructors

        //- Construct from components
        inline bPolynomial
        (
            const Specie& sp,
            const scalar rhoRef,
            const scalar Tref,            
            const Polynomial<PolySize>& rhoPoly
        );

        //- Construct from dictionary
        explicit bPolynomial(const dictionary& dict);

        //- Construct as named copy
        inline bPolynomial(const word& name, const bPolynomial&);

        //- Construct and return a clone
        inline autoPtr<bPolynomial> clone() const;

        // Selector from dictionary
        inline static autoPtr<bPolynomial> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "bPolynomial<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const bPolynomial&);
        inline void operator*=(const scalar);


    // Friend Operators

        friend bPolynomial operator+ <Specie, PolySize>
        (
            const bPolynomial&,
            const bPolynomial&
        );

        friend bPolynomial operator* <Specie, PolySize>
        (
            const scalar s,
            const bPolynomial&
        );

        friend bPolynomial operator== <Specie, PolySize>
        (
            const bPolynomial&,
            const bPolynomial&
        );


    // IOstream Operators

        friend Ostream& operator<< <Specie, PolySize>
        (
            Ostream&,
            const bPolynomial&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makebPolynomial(PolySize)                                            \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
(                                                                              \
    bPolynomial<Specie, PolySize>,                                           \
    "bPolynomial<"#PolySize">",                                              \
    0                                                                          \
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "bPolynomialI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "bPolynomial.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
