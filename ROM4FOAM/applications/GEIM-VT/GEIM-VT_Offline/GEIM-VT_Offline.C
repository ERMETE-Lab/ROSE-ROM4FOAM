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
    GEIM-VT_Offline

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
#include "simpleMatrix.H"
// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //
#include "ReadGEIM-VTsolverDict.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- update B matrix (interpolation Matrix)
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

void assembleSourceTerm
(
    scalarField& source,
    const GeometricField<scalar, fvPatchField, volMesh>& field_to_interpolate,
    const PtrList<GeometricField<scalar, fvPatchField, volMesh>>& MagicSensors
)
{
    label ms_size = MagicSensors.size();
    source.setSize(ms_size);
    for (label zz=0; zz< ms_size; ++zz)
    {
        source[zz]= MOR::projection(field_to_interpolate, MagicSensors[zz]) ;
    }
}

int main(int argc, char *argv[])
{

    auto start_s=std::clock();

#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"


    //get GEIM_VT_Offline parameters GEIM-VTsolverDict

    GEIM_VTparameters GEIM_VT_parameters = getGEIMOfflineParameters(args);

    List<fileName> foldersList = GEIM_VT_parameters.folders_list ;

    pointField sensorsPositions = GEIM_VT_parameters.sensorsPositions;

    label MaxSensorsNumber = GEIM_VT_parameters.MaxSensorsNumber ;

    scalar SensorsVariance = GEIM_VT_parameters.SensorsVariance;

    scalar desidered_error = GEIM_VT_parameters.error;

    // Initialise all Lists necessary to store the snapshots

    PtrList<volScalarField> T_snapshots;
    PtrList<volScalarField> p_snapshots;
    PtrList<volVectorField> U_snapshots;

    // get Snapshots

    for (label folderI=0; folderI < foldersList.size(); ++ folderI)
    {
        chDir(args.rootPath()/foldersList[folderI]);

        Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), foldersList[folderI]);

#include "CreateMesh.H"

        Info << "reading snapshots in "<< foldersList[folderI]<<"\n"<< endl;

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

            T_snapshots.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "T",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );

            p_snapshots.append
            (
                new volScalarField
                (
                    IOobject
                    (
                        "p_rgh",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );

            U_snapshots.append
            (
                new volVectorField
                (
                    IOobject
                    (
                        "U",
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


    // Create GEIM-VT_s_(SensorsVariance) folder

#include "createFolderResultsandResetTime.H"
#include "CreateMesh.H"

    // Magic Quantities declaration:

    PtrList<volScalarField> T_MagicFunctions;
    PtrList<volScalarField> p_MagicFunctions;
    PtrList<volVectorField> U_MagicFunctions;
    PtrList<volScalarField> MagicSensors;

    // definition of temp varibles needed for the loop

    scalar max_Norm=SMALL;
    scalar max_measure=SMALL;
    label tmp_sensID=0;
    label iterIndex=0;
    label generatingFunction=0;
    List<scalar> relativeErrorList;
    scalarSquareMatrix Bmatrix;
    List<word> whichFieldisCorrected;
    List<scalarField> InterpolationCoeffList;

    /***************************** GEIM-VT First Iteration ***************************/

    //update index
    ++iterIndex;

    whichFieldisCorrected.setSize(iterIndex);

    forAll (T_snapshots, snapI)
    {
        if(max_Norm<MOR::L2norm(T_snapshots[snapI]))
        {
            max_Norm = MOR::L2norm(T_snapshots[snapI]);
            generatingFunction=snapI;
            whichFieldisCorrected[iterIndex-1]=T_snapshots[snapI].name();
        }
    }

    forAll (p_snapshots, snapI)
    {
        if(max_Norm<MOR::L2norm(p_snapshots[snapI]))
        {
            max_Norm = MOR::L2norm(p_snapshots[snapI]);
            generatingFunction=snapI;
            whichFieldisCorrected[iterIndex-1]=p_snapshots[snapI].name();
        }
    }

    forAll (U_snapshots, snapI)
    {
        if(max_Norm<MOR::L2norm(U_snapshots[snapI]))
        {
            max_Norm =MOR::L2norm(U_snapshots[snapI]);
            generatingFunction=snapI;
            whichFieldisCorrected[iterIndex-1]=U_snapshots[snapI].name();
        }
    }


    Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
    Info <<"maximum relative interpolation error for field "<< whichFieldisCorrected[iterIndex-1]  << " :     " << max_Norm/max_Norm  << endl;

    //update errorLists
    relativeErrorList.append(max_Norm/max_Norm);

    // find magic functional

    forAll(sensors, sensI)
    {
        if (mag(max_measure) < mag(MOR::projection(T_snapshots[generatingFunction],sensors[sensI])) )
        {
            max_measure= MOR::projection(T_snapshots[generatingFunction],sensors[sensI]) ;
            tmp_sensID=sensI;
        }
    }


    // Update Magic Quantites

    T_MagicFunctions.append
    (
        GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                T_snapshots[generatingFunction].name(),
                T_snapshots[generatingFunction].time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            T_snapshots[generatingFunction],
            fvPatchField<scalar>::calculatedType()
        )
    );


    p_MagicFunctions.append
    (
        GeometricField<scalar, fvPatchField, volMesh>
        (
            IOobject
            (
                p_snapshots[generatingFunction].name(),
                p_snapshots[generatingFunction].time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            p_snapshots[generatingFunction],
            fvPatchField<scalar>::calculatedType()
        )
    );
    U_MagicFunctions.append
    (
        GeometricField<vector, fvPatchField, volMesh>
        (
            IOobject
            (
                U_snapshots[generatingFunction].name(),
                U_snapshots[generatingFunction].time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_snapshots[generatingFunction],
            fvPatchField<vector>::calculatedType()
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


    /******************************* GEIM-VT MAIN LOOP ********************************/

    scalar error=desidered_error+1;
    InterpolationCoeffList.setSize(T_snapshots.size());
    label coeffIndex=0;


    while ( (error>desidered_error) && (iterIndex<=MaxSensorsNumber))
    {
        /* * * * update loop index* * * */

        ++iterIndex;

        whichFieldisCorrected.setSize(iterIndex);

        /* * * * reset temparary parameters* * * */

        error=SMALL;
        max_measure=SMALL;
        coeffIndex=iterIndex-2;
        // declare pointer which will contain the magic function
        volScalarField* T_tmp_residual_Ptr=new GeometricField<scalar, fvPatchField, volMesh>
        (
            T_snapshots[0]
        );
        volScalarField* p_tmp_residual_Ptr=new GeometricField<scalar, fvPatchField, volMesh>
        (
            p_snapshots[0]
        );
        volVectorField* U_tmp_residual_Ptr=new GeometricField<vector, fvPatchField, volMesh>
        (
            U_snapshots[0]
        );

        /* update interpolation Matrix*/

        updateBmatrix (Bmatrix, T_MagicFunctions, MagicSensors);

        forAll (T_snapshots, snapI)
        {

            /* find new GEIM coefficient for snapshot snapI */

            InterpolationCoeffList[snapI].setSize(coeffIndex+1);
            InterpolationCoeffList[snapI][coeffIndex] =MOR::projection(T_snapshots[snapI], MagicSensors[coeffIndex]);

            if ((coeffIndex)!=0)  // fist iteration alpha1 = b1/B11
            {
                for (label ii=0; ii< MagicSensors.size()-1; ii++)
                {
                    InterpolationCoeffList[snapI][coeffIndex]-= (InterpolationCoeffList[snapI][ii]*Bmatrix(coeffIndex, ii));

                }
            }

            InterpolationCoeffList[snapI][coeffIndex]/=Bmatrix(coeffIndex,coeffIndex);


            GeometricField<scalar, fvPatchField, volMesh> T_interpolant
            (
                IOobject
                (
                    T_snapshots[snapI].name() + "GEIMInterpolation",
                    T_snapshots[snapI].time().timeName(),
                    T_snapshots[snapI].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                T_snapshots[snapI].mesh(),
                dimensioned<scalar>
                (
                    "zero",
                    T_snapshots[snapI].dimensions(),
                    pTraits<scalar>::zero
                )
            );

            GeometricField<scalar, fvPatchField, volMesh> p_interpolant
            (
                IOobject
                (
                    p_snapshots[snapI].name() + "GEIMInterpolation",
                    p_snapshots[snapI].time().timeName(),
                    p_snapshots[snapI].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                p_snapshots[snapI].mesh(),
                dimensioned<scalar>
                (
                    "zero",
                    p_snapshots[snapI].dimensions(),
                    pTraits<scalar>::zero
                )
            );

            GeometricField<vector, fvPatchField, volMesh> U_interpolant
            (
                IOobject
                (
                    U_snapshots[snapI].name() + "GEIMInterpolation",
                    U_snapshots[snapI].time().timeName(),
                    U_snapshots[snapI].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U_snapshots[snapI].mesh(),
                dimensioned<vector>
                (
                    "zero",
                    U_snapshots[snapI].dimensions(),
                    pTraits<vector>::zero
                )
            );

            forAll (InterpolationCoeffList[snapI], mfI)
            {

                T_interpolant+=
                    InterpolationCoeffList[snapI][mfI]*T_MagicFunctions[mfI];
                p_interpolant+=
                    InterpolationCoeffList[snapI][mfI]*p_MagicFunctions[mfI];
                U_interpolant+=
                    InterpolationCoeffList[snapI][mfI]*U_MagicFunctions[mfI];
            }

            GeometricField<scalar, fvPatchField, volMesh> T_residual
            (
                (T_snapshots[snapI] - T_interpolant)
            );

            GeometricField<scalar, fvPatchField, volMesh> p_residual
            (
                (p_snapshots[snapI] - p_interpolant)
            );

            GeometricField<vector, fvPatchField, volMesh> U_residual
            (
                (U_snapshots[snapI] - U_interpolant)
            );

            scalar T_snapshots_measure = MOR::L2norm(T_snapshots[snapI]);
            if (T_snapshots_measure < SMALL)
            {
                T_snapshots_measure=SMALL;
            }
            scalar p_snapshots_measure = MOR::L2norm(p_snapshots[snapI]);
            if (p_snapshots_measure < SMALL)
            {
                p_snapshots_measure=SMALL;
            }
            scalar U_snapshots_measure = MOR::L2norm(U_snapshots[snapI]);
            if (U_snapshots_measure < SMALL)
            {
                U_snapshots_measure=SMALL;
            }

            scalar relative_T_interpolation_error =MOR::L2norm(T_residual)/T_snapshots_measure;
            scalar relative_p_interpolation_error =MOR::L2norm(p_residual)/p_snapshots_measure;
            scalar relative_U_interpolation_error =MOR::L2norm(U_residual)/U_snapshots_measure;

            if ((error< relative_T_interpolation_error) || (error< relative_U_interpolation_error) || (error< relative_p_interpolation_error))
            {
                generatingFunction=snapI;
                T_tmp_residual_Ptr-> ~GeometricField<scalar, fvPatchField, volMesh>();
                T_tmp_residual_Ptr = new(T_tmp_residual_Ptr) GeometricField<scalar, fvPatchField, volMesh>
                (
                    T_residual
                );
                p_tmp_residual_Ptr-> ~GeometricField<scalar, fvPatchField, volMesh>();
                p_tmp_residual_Ptr = new(p_tmp_residual_Ptr) GeometricField<scalar, fvPatchField, volMesh>
                (
                    p_residual
                );
                U_tmp_residual_Ptr-> ~GeometricField<vector, fvPatchField, volMesh>();
                U_tmp_residual_Ptr = new(U_tmp_residual_Ptr) GeometricField<vector, fvPatchField, volMesh>
                (
                    U_residual
                );

                if (error< relative_T_interpolation_error)
                {
                    error=relative_T_interpolation_error;
                    whichFieldisCorrected[iterIndex-1]=T_snapshots[snapI].name();
                }

                if (error< relative_p_interpolation_error)
                {
                    error=relative_p_interpolation_error;
                    whichFieldisCorrected[iterIndex-1]=p_snapshots[snapI].name();
                }

                if(error< relative_U_interpolation_error)
                {
                    error=relative_U_interpolation_error;
                    whichFieldisCorrected[iterIndex-1]=U_snapshots[snapI].name();
                }
            }

        }

        /* * * *find magic functional * * * */
        if ((error > desidered_error) && (iterIndex <= MaxSensorsNumber))
        {
            forAll(sensors, sensI)
            {
                if (sensors(sensI))
                {
                    if (mag(max_measure) < mag(MOR::projection(*T_tmp_residual_Ptr,sensors[sensI])) )
                    {
                        max_measure=MOR::projection(*T_tmp_residual_Ptr,sensors[sensI]);
                        tmp_sensID=sensI;
                    }
                }
            }


            // Update Magic Quantites

            T_MagicFunctions.append(T_tmp_residual_Ptr);
            p_MagicFunctions.append(p_tmp_residual_Ptr);
            U_MagicFunctions.append(U_tmp_residual_Ptr);

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
            delete T_tmp_residual_Ptr;
            delete p_tmp_residual_Ptr;
            delete U_tmp_residual_Ptr;
        }

        //update errorLists

        relativeErrorList.append(error);

        /* * * * print info * * * */
        Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
        Info <<"maximum relative interpolation error for field "<< whichFieldisCorrected[iterIndex-1]  << " :     " << setprecision(14)<<  error  <<endl;

    }
    
    /*************************** Print Basis, Sensors & files******************************/

    /*reset IO and wite MagicFunction*/

    for(label mfI=0; mfI<T_MagicFunctions.size(); mfI++)
    {
        T_MagicFunctions.set
        (   mfI,
            new volScalarField
            (   IOobject
                (
                    "TGEIM-VTMagicFunction"+name(mfI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                T_MagicFunctions[mfI],
                fvPatchField<scalar>::calculatedType()
            )
        );

        Info <<"\nWriting "<< T_MagicFunctions[mfI].name()<<" in folder: "<<runTime.timeName() << endl;
        T_MagicFunctions(mfI)->write();
    }

    for(label mfI=0; mfI<p_MagicFunctions.size(); mfI++)
    {
        p_MagicFunctions.set
        (   mfI,
            new volScalarField
            (   IOobject
                (
                    "p_rghGEIM-VTMagicFunction"+name(mfI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                p_MagicFunctions[mfI],
                fvPatchField<scalar>::calculatedType()
            )
        );

        Info <<"\nWriting "<< p_MagicFunctions[mfI].name()<<" in folder: "<<runTime.timeName() << endl;
        p_MagicFunctions(mfI)->write();
    }

    for(label mfI=0; mfI<U_MagicFunctions.size(); mfI++)
    {
        U_MagicFunctions.set
        (   mfI,
            new volVectorField
            (   IOobject
                (
                    "UGEIM-VTMagicFunction"+name(mfI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                U_MagicFunctions[mfI],
                fvPatchField<scalar>::calculatedType()
            )
        );

        Info <<"\nWriting "<< U_MagicFunctions[mfI].name()<<" in folder: "<<runTime.timeName() << endl;
        U_MagicFunctions(mfI)->write();
    }

    for(label sensI=0; sensI<MagicSensors.size(); sensI++)
    {

        MagicSensors.set
        (
            sensI,
            new volScalarField
            (   IOobject
                (
                    "TMagicSensor"+name(sensI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                MagicSensors[sensI],
                fvPatchField<scalar>::calculatedType()
            )
        );

        Info <<"\nWriting "<< MagicSensors[sensI].name()<<" in folder: "<<runTime.timeName() << endl;
        MagicSensors[sensI].write();

    }

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
