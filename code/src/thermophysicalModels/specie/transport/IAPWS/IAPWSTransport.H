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
    Foam::IAPWSTransport

Group
    grpSpecieTransport

Description
    Transport package using polynomial functions for \c mu and \c kappa.

Usage

    \table
        Property        | Description
        muCoeffs<8>     | Dynamic viscosity polynomial coefficients
        kappaCoeffs<8>  | Thermal conductivity polynomial coefficients
    \endtable

    Example of the specification of the transport properties:
    \verbatim
    transport
    {
        muCoeffs<8>     ( 1000 -0.05 0.003 0 0 0 0 0 );
        kappaCoeffs<8>  ( 2000 -0.15 0.023 0 0 0 0 0 );
    }
    \endverbatim

    The polynomial expressions are evaluated as so:

        \f[
            \mu    = 1000 - 0.05 T + 0.003 T^2
        \f]

        \f[
            \kappa = 2000 - 0.15 T + 0.023 T^2
        \f]

Note
    - Dynamic viscosity polynomial coefficients evaluate to an expression in
      [Pa.s], but internally uses [Pa.s/kmol].
    - Thermal conductivity polynomial coefficients evaluate to an expression in
      [W/m/K], but internally uses [W/m/K/kmol].

SourceFiles
    IAPWSTransportI.H
    IAPWSTransport.C

See also
    Foam::Polynomial

\*---------------------------------------------------------------------------*/

#ifndef IAPWSTransport_H
#define IAPWSTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo, int PolySize> class IAPWSTransport;

template<class Thermo, int PolySize>
inline IAPWSTransport<Thermo, PolySize> operator+
(
    const IAPWSTransport<Thermo, PolySize>&,
    const IAPWSTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline IAPWSTransport<Thermo, PolySize> operator*
(
    const scalar,
    const IAPWSTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const IAPWSTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                     Class IAPWSTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=8>
class IAPWSTransport
:
    public Thermo
{
    // Private Data

        //- Dynamic viscosity polynomial coefficients
        Polynomial<PolySize> muCoeffs_;

        //- Thermal conductivity polynomial coefficients
        Polynomial<PolySize> kappaCoeffs_;


    // Private Member Functions

        //- Coeffs name. Eg, "muCoeffs<10>"
        inline static word coeffsName(const char* name)
        {
            return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        }

        //- Construct from components
        inline IAPWSTransport
        (
            const Thermo& t,
            const Polynomial<PolySize>& muPoly,
            const Polynomial<PolySize>& kappaPoly
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct as named copy
        inline IAPWSTransport(const word&, const IAPWSTransport&);

        //- Construct from dictionary
        explicit IAPWSTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<IAPWSTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<IAPWSTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "polynomial<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const IAPWSTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend IAPWSTransport operator+ <Thermo, PolySize>
        (
            const IAPWSTransport&,
            const IAPWSTransport&
        );

        friend IAPWSTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const IAPWSTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const IAPWSTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IAPWSTransportI.H"

#ifdef NoRepository
    #include "IAPWSTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
