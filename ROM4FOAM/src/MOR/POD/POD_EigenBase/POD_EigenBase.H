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
    Foam::POD_EigenBase

Description

SourceFiles
    POD_EigenBaseI.H
    POD_EigenBase.C
    POD_EigenBaseIO.C

\*---------------------------------------------------------------------------*/

#ifndef POD_EigenBase_H
#define POD_EigenBase_H

#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "EigenSolver.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class POD_EigenBase Declaration
\*---------------------------------------------------------------------------*/

class POD_EigenBase
{

    // Private data

        //- Eigen-values
        scalarField eigenValues_;

        //- Cumulative relative eigen-values.  Renormalised to sum to 1
        scalarField cumEigenValues_;

        //- Eigen-vectors
        FieldField<Field, scalar> eigenVectors_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        POD_EigenBase(const POD_EigenBase&);

        //- Disallow default bitwise assignment
        void operator=(const POD_EigenBase&);


        //- Calculate eigen base
        void calcEigenBase(const scalarSquareMatrix& orthMatrix);


public:

    // Constructors

        //- Construct given a list of fields
        template<class Type>
        POD_EigenBase(const PtrList<GeometricField<Type, fvPatchField, volMesh> >& snapshots);

       


    // Destructor - default


    // Member Functions

        //- Return eigen-values sorted in decreasing order
        const scalarField& eigenValues() const
        {
            return eigenValues_;
        }

        //- Return cumulative eigen-values
        const scalarField& cumulativeEigenValues() const
        {
            return cumEigenValues_;
        }

        //- Return eigen-vectors
        const FieldField<Field, scalar>& eigenVectors() const
        {
            return eigenVectors_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "POD_EigenBase.C"
#endif

#endif

// ************************************************************************* //


