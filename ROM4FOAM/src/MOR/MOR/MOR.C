/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
#include "fvCFD.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace MOR
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


scalar projection
(
    const volVectorField& a,
    const volVectorField& b
)
{

    return fvc::domainIntegrate(a&b).value();
  
}


scalar projection
(
    const volScalarField& a,
    const volScalarField& b
)
{
    return fvc::domainIntegrate(a*b).value();
  
}


scalar L2norm
     (
         const volVectorField& a
     )
{
  
    return Foam::sqrt(fvc::domainIntegrate(a&a).value());
    
}


scalar L2norm
     (
         const volScalarField& a
     )
{
   
    return Foam::sqrt(fvc::domainIntegrate(a*a).value());
   
}

scalar H1norm
     (
         const volScalarField& a
     )
{
    return Foam::sqrt(fvc::domainIntegrate(a*a).value()+ fvc::domainIntegrate((fvc::grad(a))&(fvc::grad(a))).value()) ; 
}

scalar H1norm
     (
         const volVectorField& a
     )
{
    return Foam::sqrt((fvc::domainIntegrate(a&a).value())+ fvc::domainIntegrate((fvc::grad(a))&&(fvc::grad(a))).value()) ; 
}


template<class Type>
scalar L_infinity_norm
     (
         const GeometricField<Type, fvPatchField, volMesh>&  a
     )
{
    return max(mag(a.internalField())).value(); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace MOR

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
