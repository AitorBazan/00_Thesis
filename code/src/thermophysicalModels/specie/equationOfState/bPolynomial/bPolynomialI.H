/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "bPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::bPolynomial<Specie, PolySize>::bPolynomial
(
    const Specie& sp,
    const scalar rhoRef, 
    const scalar Tref,
    const Polynomial<PolySize>& rhoCoeffs
)
:
    Specie(sp),
    rhoRef_(rhoRef),
    Tref_(Tref),    
    rhoCoeffs_(rhoCoeffs)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::bPolynomial<Specie, PolySize>::bPolynomial
(
    const word& name,
    const bPolynomial<Specie, PolySize>& ip
)
:
    Specie(name, ip),
    rhoRef_(ip.rhoRef_),
    Tref_(ip.Tref_),
    rhoCoeffs_(ip.rhoCoeffs_)
{}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::bPolynomial<Specie, PolySize>>
Foam::bPolynomial<Specie, PolySize>::clone() const
{
    return autoPtr<bPolynomial<Specie, PolySize>>::New(*this);
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::bPolynomial<Specie, PolySize>>
Foam::bPolynomial<Specie, PolySize>::New(const dictionary& dict)
{
    return autoPtr<bPolynomial<Specie, PolySize>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::rho
(
    scalar p,
    scalar T
) const
{
    return rhoCoeffs_.value(T-Tref_);
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::H
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::E
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::Cv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::bPolynomial<Specie, PolySize>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline void Foam::bPolynomial<Specie, PolySize>::operator+=
(
    const bPolynomial<Specie, PolySize>& ip
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ip);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = ip.Y()/this->Y();
        rhoRef_ = Y1*rhoRef_ + Y2*ip.rhoRef_;
        Tref_ = Y1*Tref_ + Y2*ip.Tref_;
        rhoCoeffs_ = Y1*rhoCoeffs_ + Y2*ip.rhoCoeffs_;
    }
}


template<class Specie, int PolySize>
inline void Foam::bPolynomial<Specie, PolySize>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::bPolynomial<Specie, PolySize> Foam::operator+
(
    const bPolynomial<Specie, PolySize>& ip1,
    const bPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
      + static_cast<const Specie&>(ip2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return bPolynomial<Specie, PolySize>
        (
            sp,
            ip1.rhoRef_,
            ip1.Tref_,            
            ip1.rhoCoeffs_
        );
    }
    else
    {
        const scalar Y1 = ip1.Y()/sp.Y();
        const scalar Y2 = ip2.Y()/sp.Y();

        return bPolynomial<Specie, PolySize>
        (
            sp,
            Y1*ip1.rhoRef_ + Y2*ip2.rhoRef_,
            Y1*ip1.Tref_ + Y2*ip2.Tref_,            
            Y1*ip1.rhoCoeffs_ + Y2*ip2.rhoCoeffs_
        );
    }
}


template<class Specie, int PolySize>
Foam::bPolynomial<Specie, PolySize> Foam::operator*
(
    const scalar s,
    const bPolynomial<Specie, PolySize>& ip
)
{
    return bPolynomial<Specie, PolySize>
    (
        s*static_cast<const Specie&>(ip),
        ip.rhoRef_,
        ip.Tref_,        
        ip.rhoCoeffs_
    );
}


template<class Specie, int PolySize>
Foam::bPolynomial<Specie, PolySize> Foam::operator==
(
    const bPolynomial<Specie, PolySize>& ip1,
    const bPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
     == static_cast<const Specie&>(ip2)
    );

    const scalar Y1 = ip1.Y()/sp.Y();
    const scalar Y2 = ip2.Y()/sp.Y();

    return bPolynomial<Specie, PolySize>
    (
        sp,
        Y2*ip2.rhoRef_ - Y1*ip1.rhoRef_,
        Y2*ip2.Tref_   - Y1*ip1.Tref_,        
        Y2*ip2.rhoCoeffs_ - Y1*ip1.rhoCoeffs_
    );
}


// ************************************************************************* //
