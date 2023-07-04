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

Class
    Foam::MOR

Description

SourceFiles
    MORI.H
    MOR.C
    MORIO.C

\*---------------------------------------------------------------------------*/

#ifndef MOR_H
#define MOR_H

#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Namespace MOR functions Declaration
\*---------------------------------------------------------------------------*/
namespace MOR
{
   
    //- L2 projection of snapshots
    scalar projection
    (
        const volVectorField& a,
        const volVectorField& b
    );
   

    scalar projection
    (
        const volScalarField& a,
        const volScalarField& b
    );



    // L2norms  
    scalar L2norm
    (
        const volVectorField& a
    );

    scalar L2norm
    (
        const volScalarField& a
    );
    
    // H1_norm
    
    scalar H1norm
    (
        const volScalarField& a
    );

    scalar H1norm
    (
        const volVectorField& a
    );
     
    // L_infinity_norm
    template<class Type> 
    scalar L_infinity_norm
    (
        const GeometricField<Type, fvPatchField, volMesh>&  a
    );
    
     
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
#   include "MOR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
