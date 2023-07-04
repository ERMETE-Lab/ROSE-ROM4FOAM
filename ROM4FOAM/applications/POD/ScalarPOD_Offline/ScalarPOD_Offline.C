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

Application
    ScalarPOD_Offline

Author
	Stefano Riva and Carolina Introini

\*---------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"
#include "POD_EigenBase.H"

#include "regionProperties.H"
#include "turbulentFluidThermoModel.H" // necessary for the BC used in fluid-solid cases

// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //

#include "ReadPODsolverDict.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    auto start_s=std::clock();
    argList::addOption
	(
	    "region",
	    "name",
	    "Specify the mesh region" 
	);
#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"

    //get POD parameters

    PODparameters POD_parameters = getPODparameters(args);

    List<fileName> folders_list = POD_parameters.folders_list ;

    word fieldName = POD_parameters.fieldName;

    label maxBasis = POD_parameters.maxBasis;

    scalar accuracy = POD_parameters.accuracy;

    // Initialise all Lists necessary to store the snapshots and PODBasis

    PtrList<volScalarField> scalarSnapshotsList;
    PtrList<volScalarField> scalarPODOrthoNormalBasis;

    // get ScalarSnapshots

    for (label folderI=0; folderI < folders_list.size(); ++ folderI)
    {
        chDir(args.rootPath()/folders_list[folderI]);

        Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), folders_list[folderI]);

#include "CreateMesh.H"

        Info << "Reading snapshots in "<< folders_list[folderI]<<"\n"<< endl;

        // Get times list
        instantList timeDirs = timeSelector::select0(runTime, args);

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            Info<< "Time = " << runTime.timeName() << endl;

            if (mesh.readUpdate() != polyMesh::UNCHANGED)
            {

                FatalErrorInFunction
                        << "polyMesh has changed. POD can be performed only on unchanged polyMesh"
                        << abort(FatalError);
            }

            scalarSnapshotsList.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

    }

    // create new folder in which all POD Basis will be stored

#include"createFolderResults.H"

#include "CreateMesh.H"

    // performPOD

    POD_EigenBase eigenBase( scalarSnapshotsList );

    label Basisize = 0;

    const scalarField& EigenValues = eigenBase.eigenValues();
    const scalarField& cumEigenValues = eigenBase.cumulativeEigenValues();

    forAll (cumEigenValues, i)
    {
        Basisize++;

        if ((cumEigenValues[i] > accuracy) | (Basisize==maxBasis))
        {
            break;
        }
    }

    Info<< "Cumulative eigen-values: "
        << setprecision(14) << cumEigenValues << nl
        << "Base size: " << Basisize << endl;

    /// Establish orthonormal base

    scalarPODOrthoNormalBasis.setSize(Basisize);

    for (label baseI = 0; baseI < Basisize; baseI++)
    {
        const scalarField& eigenVector = eigenBase.eigenVectors()[baseI];

        GeometricField<scalar, fvPatchField, volMesh> * onBasePtr
        (
            new GeometricField<scalar, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldName + "POD" + name(baseI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                mesh,
                dimensioned<scalar>
                (
                    "zero",
                    scalarSnapshotsList[baseI].dimensions(),
                    pTraits<scalar>::zero
                )
            )
        );

        GeometricField<scalar, fvPatchField, volMesh>& onBase = *onBasePtr;

        forAll (eigenVector, eigenI)
        {
            onBase += eigenVector[eigenI]*scalarSnapshotsList[eigenI];
        }

        /// Re-normalise ortho-normal vector

        scalar SqrtEigenValue = Foam::sqrt(EigenValues[baseI]);
        if (SqrtEigenValue > SMALL)
        {
            onBase /= SqrtEigenValue;
        }

        scalarPODOrthoNormalBasis.set(baseI, onBasePtr);
    }

    // Calculate interpolation coefficients

    scalarRectangularMatrix coeffs(scalarSnapshotsList.size(), scalarPODOrthoNormalBasis.size());

    forAll (scalarSnapshotsList, snapshotI)
    {
        forAll (scalarPODOrthoNormalBasis, baseI)
        {
            coeffs[snapshotI][baseI] =
                MOR::projection
                (
                    scalarSnapshotsList[snapshotI],
                    scalarPODOrthoNormalBasis[baseI]
                );
        }
    }

    // Check all snapshots
	scalarRectangularMatrix L2AbsErrorMatrix (scalarSnapshotsList.size(), scalarPODOrthoNormalBasis.size());
	scalarRectangularMatrix L2RelErrorMatrix (scalarSnapshotsList.size(), scalarPODOrthoNormalBasis.size());
    
    forAll (scalarSnapshotsList, snapI)
    {

        volScalarField reconstruct
        (
            IOobject
            (
                fieldName + "PODreconstruct",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensioned<scalar>
            (
                "zero",
                scalarSnapshotsList[snapI].dimensions(),
                pTraits<scalar>::zero
            )
        );

        for (label baseI = 0; baseI < scalarPODOrthoNormalBasis.size(); baseI++)
        {
            reconstruct +=
                coeffs[snapI][baseI]*scalarPODOrthoNormalBasis[baseI];
        
            scalar sumFieldError =
                MOR::L2norm
                (
                    reconstruct - scalarSnapshotsList[snapI]
                );

            scalar measure =
                MOR::L2norm(scalarSnapshotsList[snapI]);

            if (measure< SMALL)
            {
                measure = SMALL;
            }

            L2AbsErrorMatrix[snapI][baseI] = sumFieldError;
            L2RelErrorMatrix[snapI][baseI] = sumFieldError/measure;

        }
    }

    List<scalar> L2AbsErrorList ( scalarPODOrthoNormalBasis.size() );
    List<scalar> L2RelErrorList ( scalarPODOrthoNormalBasis.size() );


    forAll (L2AbsErrorList, baseI)
    {
        scalarField AbsTmpVector ( scalarSnapshotsList.size() );
        scalarField RelTmpVector ( scalarSnapshotsList.size() );

        forAll(scalarSnapshotsList, snapI)
        {
            AbsTmpVector[snapI] = L2AbsErrorMatrix[snapI][baseI];
            RelTmpVector[snapI] = L2RelErrorMatrix[snapI][baseI];
        }

        L2AbsErrorList[baseI] = max(AbsTmpVector);
        L2RelErrorList[baseI] = max(RelTmpVector);
    }
    // Write out all fields

    Info << "Writing POD scalar base for volScalarField "
         << fieldName << " in folder '"
         << runTime.timeName()<<"'"  << endl;

    for (label i = 0; i < scalarPODOrthoNormalBasis.size(); i++)
    {
        scalarPODOrthoNormalBasis(i)->write();
    }

    // Write out all fields

#include "printAllFiles.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    auto stop_s=std::clock();

    Info<< nl << "ExecutionTime = " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
