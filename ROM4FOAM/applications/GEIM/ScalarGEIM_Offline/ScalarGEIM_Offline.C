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
    ScalarGEIM_Offline

Author
	Stefano Riva and Carolina Introini

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"
#include "EigenSolver.H"
#include "SortableList.H"

#include "regionProperties.H"
#include "turbulentFluidThermoModel.H" // necessary for the BC used in fluid-solid cases

// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //
#include "ReadGEIMsolverDict.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- update B matrix (interpolation matrix)
void updateBmatrix
(
    scalarSquareMatrix& BMatrix,
    const PtrList<GeometricField<scalar, fvPatchField, volMesh>>& magic_functions,
    const PtrList<GeometricField<scalar, fvPatchField, volMesh>>& MagicSensors
)
{
    label mf_size = magic_functions.size();
    label ms_size = MagicSensors.size();
    if ( ms_size!=mf_size)
    {
        FatalErrorInFunction
                << "Error:magic functions and magic points must have same size"
                << abort(FatalError);
    }

    BMatrix.setSize(mf_size);

    for (label ii=0; ii< mf_size-1; ++ii)
    {
        BMatrix(ii,mf_size-1)=0;
    }

    for (label jj=0; jj< mf_size; ++jj)
    {
        BMatrix(mf_size-1,jj)=MOR::projection(magic_functions[jj], MagicSensors[mf_size-1]);

    }

}


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


    //get GEIM_Offline parameters from GEIMsolverDict

    GEIMparameters GEIM_parameters = getGEIMOfflineParameters(args);

    word fieldName = GEIM_parameters.fieldName ;

    List<fileName> foldersList = GEIM_parameters.folders_list ;

    label MaxSensorsNumber = GEIM_parameters.MaxSensorsNumber ;

    pointField sensorsPositions = GEIM_parameters.sensorsPositions;

    scalar SensorsVariance = GEIM_parameters.SensorsVariance;

    scalar desidered_error = GEIM_parameters.error;

    // Initialise List necessary to store the snapshots 

    PtrList<volScalarField> scalarSnapshotsList;

    // get ScalarSnapshots

    for (label folderI=0; folderI < foldersList.size(); ++ folderI)
    {
        chDir(args.rootPath()/foldersList[folderI]);

        Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), foldersList[folderI]);
        const wordList regionNames(selectRegionNames(args, runTime));
        if (regionNames.size() > 1)
        {
            Info << "Error: more than one region selected!" << endl;
            abort();
        }
        const word& region = regionNames[0];
#include "CreateMesh.H"

        Info << "\nReading snapshots in "<< foldersList[folderI]<< endl;

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

    // generate Sensors

#include "generateSensors.H"

    // Create GEIM_(fieldName)_s_(SensorsVariance) folder

#include "createFolderResultsandResetTime.H"
        const wordList regionNames(selectRegionNames(args, runTime));
        if (regionNames.size() > 1)
        {
            Info << "Error: more than one region selected!" << endl;
            abort();
        }
        const word& region = regionNames[0];
#include "CreateMesh.H"

    // Magic Quantities declaration:

    PtrList<volScalarField> MagicFunctions;
    PtrList<volScalarField> MagicSensors;

    // definition of temp variables needed for the loop

    scalar max_L2norm=SMALL;
    scalar max_measure=SMALL;
    label tmp_sensID=0;
    label iterIndex=0;
    label generatingFunction=0;
    List<scalarField> InterpolationCoeffList;
    List<scalar> absoluteErrorList;
    List<scalar> relativeErrorList;
    scalarSquareMatrix Bmatrix;
    
    /***************************** GEIM First Iteration ***************************/
    
    //update index
    ++iterIndex;

    forAll (scalarSnapshotsList, snapI)
    {

        if(max_L2norm< MOR::L2norm(scalarSnapshotsList[snapI]))
        {
            max_L2norm = MOR::L2norm(scalarSnapshotsList[snapI]);
            generatingFunction=snapI;
        }

    }

    Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
    Info <<"maximum absolute interpolation error in L2 norm:     " << setprecision(14) << max_L2norm <<endl;
    Info <<"maximum relative interpolation error in L2 norm:     " << max_L2norm/max_L2norm  <<endl;

    //update errorLists
    absoluteErrorList.append(max_L2norm);
    relativeErrorList.append(max_L2norm/max_L2norm);

    // find magic functional

    forAll(sensors, sensI)
    {
        if (mag(max_measure) < mag(MOR::projection(scalarSnapshotsList[generatingFunction],sensors[sensI])) )
        {
            max_measure= MOR::projection(scalarSnapshotsList[generatingFunction],sensors[sensI]) ;
            tmp_sensID=sensI;
        }
    }

    // Update Magic Quantites
    MagicFunctions.append
    (
        GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                scalarSnapshotsList[generatingFunction].name(),
                scalarSnapshotsList[generatingFunction].time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            scalarSnapshotsList[generatingFunction]/max_measure,
            fvPatchField<scalar>::calculatedType()
        )
    );
    MagicSensors.append
    (
        sensors.set
        (
            tmp_sensID,
            NULL
        )
    );

    /******************************* GEIM MAIN LOOP ********************************/

    scalar max_relative_error=desidered_error+1;
    scalar max_absolute_error ;
    InterpolationCoeffList.setSize(scalarSnapshotsList.size());
    scalar snapshotsMeasure = SMALL;
    label coeffIndex=0;

    while ( (max_relative_error > desidered_error) && (iterIndex <= MaxSensorsNumber))
    {
        /* * * * update loop index* * * */

        ++iterIndex;

        /* * * * reset temparary parameters* * * */

        max_absolute_error=SMALL;
        max_relative_error=SMALL;
        max_measure=SMALL;
        volScalarField* tmp_residual_Ptr= new GeometricField<scalar, fvPatchField, volMesh>
        (
            scalarSnapshotsList[0]
        );

        coeffIndex=iterIndex-2;

        /* update interpolation Matrix*/

        updateBmatrix (Bmatrix, MagicFunctions, MagicSensors);

        /* * * * get new residual field* * * */
        forAll (scalarSnapshotsList, snapI)
        {

            /* find new GEIM coefficient for snapshot snapI */

            InterpolationCoeffList[snapI].setSize(coeffIndex+1);
            InterpolationCoeffList[snapI][coeffIndex] =MOR::projection(scalarSnapshotsList[snapI], MagicSensors[coeffIndex]);

            if ((coeffIndex)!=0)  // fist iteration alpha1 = b1/B11
            {
                for (label ii=0; ii< MagicSensors.size()-1; ii++)
                {
                    InterpolationCoeffList[snapI][coeffIndex]-= (InterpolationCoeffList[snapI][ii]*Bmatrix(coeffIndex, ii));

                }
            }


            GeometricField<scalar, fvPatchField, volMesh> interpolant
            (
                IOobject
                (
                    scalarSnapshotsList[snapI].name() + "GEIMInterpolation",
                    scalarSnapshotsList[snapI].time().timeName(),
                    scalarSnapshotsList[snapI].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                scalarSnapshotsList[snapI].mesh(),
                dimensioned<scalar>
                (
                    "zero",
                    scalarSnapshotsList[snapI].dimensions(),
                    pTraits<scalar>::zero
                )
            );


            forAll (InterpolationCoeffList[snapI], mfI)
            {

                interpolant+=
                    InterpolationCoeffList[snapI][mfI]*MagicFunctions[mfI];
            }

            GeometricField<scalar, fvPatchField, volMesh> residual
            (
                (scalarSnapshotsList[snapI] - interpolant)
            );


            scalar absolute_interpolation_error =MOR::L2norm(residual);
            snapshotsMeasure = MOR::L2norm(scalarSnapshotsList[snapI]);
            if (snapshotsMeasure < SMALL)
            {
                snapshotsMeasure=SMALL;
            }
            scalar relative_interpolation_error = absolute_interpolation_error/snapshotsMeasure;


            if (max_absolute_error< absolute_interpolation_error)
            {
                tmp_residual_Ptr-> ~GeometricField<scalar, fvPatchField, volMesh>();
                max_absolute_error=absolute_interpolation_error;
                tmp_residual_Ptr = new(tmp_residual_Ptr) GeometricField<scalar, fvPatchField, volMesh>
                (
                    residual
                );
            }

            if (max_relative_error< relative_interpolation_error)
            {
                max_relative_error=relative_interpolation_error;
            }

        }

        if ((max_relative_error > desidered_error) && (iterIndex <= MaxSensorsNumber))
        {
            /* * * *find magic functional * * * */

            forAll(sensors, sensI)
            {
                if (sensors(sensI))
                {
                    if (mag(max_measure) < mag(MOR::projection(*tmp_residual_Ptr,sensors[sensI])) )
                    {
                        max_measure= MOR::projection(*tmp_residual_Ptr,sensors[sensI]);
                        tmp_sensID=sensI;
                    }
                }
            }

            /* * * * Normalise MagicFUnction * * * */

            *tmp_residual_Ptr/=max_measure;

            /* * * * Update Magic Quantites * * * */

            MagicFunctions.append (tmp_residual_Ptr);
            MagicSensors.append
            (
                sensors.set
                (
                    tmp_sensID,
                    NULL
                )
            );

        }

        if (iterIndex > MaxSensorsNumber)
        {
            delete tmp_residual_Ptr;
        }

        /* * * * update errorLists* * * */
        absoluteErrorList.append(max_absolute_error);
        relativeErrorList.append(max_relative_error);

        /* * * * print info * * * */

        Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
        Info<<"maximum absolute interpolation error in L2 norm:     " <<  max_absolute_error << endl;
        Info<<"maximum relative interpolation error in L2 norm:     " <<  max_relative_error <<endl;
    }
    
    /*************************** Print Basis, Sensors & files******************************/

    /*reset IO and wite MagicFunction*/
    Info <<"\nWriting MagicFunctions in folder: "<<currentFolder << endl;
    for(label mfI=0; mfI<MagicFunctions.size(); mfI++)
    {
        MagicFunctions.set
        (   mfI,
            new volScalarField
            (   IOobject
                (
                    fieldName+"GEIMMagicFunction"+name(mfI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                MagicFunctions[mfI],
                fvPatchField<scalar>::calculatedType()
            )
        );


        MagicFunctions(mfI)->write();
    }
    Info <<"\nWriting MagicSensors in folder: "<<currentFolder << endl;
    for(label sensI=0; sensI<MagicSensors.size(); sensI++)
    {

        MagicSensors.set
        (
            sensI,
            new volScalarField
            (   IOobject
                (
                    fieldName+"MagicSensor"+name(sensI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                MagicSensors[sensI],
                fvPatchField<scalar>::calculatedType()
            )
        );

        
        MagicSensors[sensI].write();

    }

volScalarField MagicSensorToPlot
(   IOobject
    (
        fieldName+"MagicSensorToPlot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    MagicSensors[0],
    fvPatchField<scalar>::calculatedType()
);
for(label sensI = 1; sensI<5; sensI++)
{
    MagicSensorToPlot += MagicSensors[sensI];
}
Info <<"\nWriting "<< MagicSensorToPlot.name()<<" in folder: "<<runTime.timeName() << endl;
MagicSensorToPlot.write();


Info<<"\nComputing Lebesgue Constant"<<endl;
#include "computeLebesgueConstant.H"

#include "printAllFiles.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    auto stop_s=std::clock();

    Info<< nl << "ExecutionTime = " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
