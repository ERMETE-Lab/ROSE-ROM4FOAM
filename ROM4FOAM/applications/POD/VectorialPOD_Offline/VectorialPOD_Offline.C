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
    VectorialPOD_Offline

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

    PtrList<volVectorField> vectorSnapshotsList;
    PtrList<volVectorField> vectorPODOrthoNormalBasis;

    // get ScalarSnapshots

    for (label folderI=0; folderI < folders_list.size(); ++ folderI)
    {
        chDir(args.rootPath()/folders_list[folderI]);

        Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), folders_list[folderI]);

#include "CreateMesh.H"

        Info << "\nReading snapshots in "<< folders_list[folderI]<< endl;

        // Get times list
        instantList timeDirs = timeSelector::select0(runTime, args);

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            Info<< "    Time = " << runTime.timeName() << endl;

            if (mesh.readUpdate() != polyMesh::UNCHANGED)
            {

                FatalErrorInFunction
                        << "polyMesh has changed. POD can be performed only on unchanged polyMesh"
                        << abort(FatalError);
            }

            vectorSnapshotsList.append
            (
                new volVectorField
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

    POD_EigenBase eigenBase( vectorSnapshotsList );

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

    vectorPODOrthoNormalBasis.setSize(Basisize);

    for (label baseI = 0; baseI < Basisize; baseI++)
    {
        const scalarField& eigenVector = eigenBase.eigenVectors()[baseI];

        GeometricField<vector, fvPatchField, volMesh> * onBasePtr
        (
            new GeometricField<vector, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldName + "POD" + name(baseI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                mesh,
                dimensioned<vector>
                (
                    "zero",
                    vectorSnapshotsList[baseI].dimensions(),
                    pTraits<vector>::zero
                )
            )
        );

        GeometricField<vector, fvPatchField, volMesh>& onBase = *onBasePtr;

        forAll (eigenVector, eigenI)
        {
            onBase += eigenVector[eigenI]*vectorSnapshotsList[eigenI];
        }

        /// Re-normalise ortho-normal vector

        scalar SqrtEigenValue = Foam::sqrt(EigenValues[baseI]);
        if (SqrtEigenValue > SMALL)
        {
            onBase /= SqrtEigenValue;
        }

        vectorPODOrthoNormalBasis.set(baseI, onBasePtr);
    }

    // Calculate interpolation coefficients

    scalarRectangularMatrix coeffs(vectorSnapshotsList.size(), vectorPODOrthoNormalBasis.size());

    forAll (vectorSnapshotsList, snapshotI)
    {
        forAll (vectorPODOrthoNormalBasis, baseI)
        {
            coeffs[snapshotI][baseI] =
                MOR::projection
                (
                    vectorSnapshotsList[snapshotI],
                    vectorPODOrthoNormalBasis[baseI]
                );
        }
    }

    // Check all snapshots
    scalarRectangularMatrix L2AbsErrorMatrix (vectorSnapshotsList.size(), vectorPODOrthoNormalBasis.size());
	scalarRectangularMatrix L2RelErrorMatrix (vectorSnapshotsList.size(), vectorPODOrthoNormalBasis.size());

    forAll (vectorSnapshotsList, snapI)
    {

        volVectorField reconstruct
        (
            IOobject
            (
                fieldName + "PODreconstruct",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensioned<vector>
            (
                "zero",
                vectorSnapshotsList[snapI].dimensions(),
                pTraits<vector>::zero
            )
        );

        for (label baseI = 0; baseI < vectorPODOrthoNormalBasis.size(); baseI++)
        {
            reconstruct +=
                coeffs[snapI][baseI]*vectorPODOrthoNormalBasis[baseI];
        


        scalar sumFieldError =
            MOR::L2norm
            (
                reconstruct - vectorSnapshotsList[snapI]
            );

        scalar measure =
            MOR::L2norm(vectorSnapshotsList[snapI]);

        if (measure< SMALL)
        {
            measure = SMALL;
        }

        L2AbsErrorMatrix[snapI][baseI] = sumFieldError;
        L2RelErrorMatrix[snapI][baseI] = sumFieldError/measure;
        }

    }

    // Write out all fields
List<scalar> L2AbsErrorList ( vectorPODOrthoNormalBasis.size() );
List<scalar> L2RelErrorList ( vectorPODOrthoNormalBasis.size() );


forAll (L2AbsErrorList, baseI)
{
    scalarField AbsTmpVector ( vectorSnapshotsList.size() );
    scalarField RelTmpVector ( vectorSnapshotsList.size() );

    forAll(vectorSnapshotsList, snapI)
    {
        AbsTmpVector[snapI] = L2AbsErrorMatrix[snapI][baseI];
        RelTmpVector[snapI] = L2RelErrorMatrix[snapI][baseI];
    }

    L2AbsErrorList[baseI] = max(AbsTmpVector);
    L2RelErrorList[baseI] = max(RelTmpVector);
}

    Info << "Writing POD vector base for volVectorField "
         << fieldName << " in folder '"
         << runTime.timeName()<<"'"  << endl;

    for (label i = 0; i < vectorPODOrthoNormalBasis.size(); i++)
    {
        vectorPODOrthoNormalBasis(i)->write();
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
