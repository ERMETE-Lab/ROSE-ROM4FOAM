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
    ScalarEIM_Offline

Author
	Stefano Riva and Carolina Introini

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"
#include "regionProperties.H"
#include "turbulentFluidThermoModel.H" // necessary for the multi-region BC
// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //
#include "ReadEIMsolverDict.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- update B matrix (Interpolation Matrix)
void updateBmatrix
(
    scalarSquareMatrix& BMatrix,
    const PtrList<GeometricField<scalar, fvPatchField, volMesh>>& magic_functions,
    const labelList& magic_cellsID
)
{
    label mf_size = magic_functions.size();
    label mc_size = magic_cellsID.size();
    if ( mc_size!=mf_size)
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
        BMatrix(mf_size-1,jj)=(magic_functions[jj])[(magic_cellsID[mf_size-1])];

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


    //get EIM_Offline parameters from EIMsolverDict.

    EIMparameters EIM_parameters = getEIMparameters(args);

    word fieldName = EIM_parameters.fieldName ;

    List<fileName> foldersList = EIM_parameters.folders_list ;

    label maxBasis = EIM_parameters.maxBasis;

    scalar desidered_error = EIM_parameters.error;

    // Initialize List necessary to store the snapshots

    PtrList<volScalarField> scalarSnapshotsList;

    // get ScalarSnapshots

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
                        << "polyMesh has changed. EIM can be performed only on unchanged polyMesh"
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


    // Create EIM_(fieldName) folder

#include "createFolderResultsandResetTime.H"
#include "CreateMesh.H"

    // Magic Quantities declaration:

    PtrList<volScalarField> MagicFunctions;
    labelList MagicCellsID;

    // definition of the variables needed for the EIM Offline loop

    scalar max_L_infinity_norm=SMALL;
    scalar valueAtMagicPoint=SMALL;
    label  CellID=0;
    label  iterIndex=0;
    label  generatingFunction=0;
    List<scalarField> InterpolationCoeffList;
    List<scalar> absoluteErrorList;
    List<scalar> relativeErrorList;
    List<scalar> idxGeneratingFunction;
    scalarSquareMatrix Bmatrix;

    /***************************** EIM First Iteration ****************************/

    ++iterIndex;

    forAll (scalarSnapshotsList, snapI)
    {

        if(max_L_infinity_norm< MOR::L_infinity_norm(scalarSnapshotsList[snapI]))
        {
            max_L_infinity_norm = MOR::L_infinity_norm(scalarSnapshotsList[snapI]);
            generatingFunction=snapI;
        }

    }
    idxGeneratingFunction.append(generatingFunction);

    Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
    Info <<"maximum absolute interpolation error in Linfinty norm:     " << setprecision(14) << max_L_infinity_norm <<endl;
    Info <<"maximum relative interpolation error in Linfinty norm:     " << max_L_infinity_norm/max_L_infinity_norm <<endl;

    /* * * * update errorLists * * * */

    absoluteErrorList.append(max_L_infinity_norm);
    relativeErrorList.append(max_L_infinity_norm/max_L_infinity_norm);

    /* * * * find first magic point * * * */

    forAll(scalarSnapshotsList[generatingFunction], CellI)
    {
        if (mag(valueAtMagicPoint)< mag((scalarSnapshotsList[generatingFunction])[CellI]))
        {
            valueAtMagicPoint=(scalarSnapshotsList[generatingFunction])[CellI];
            CellID=CellI;
        }
    }

    /* * * * Update Magic Quantities * * * */

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
            scalarSnapshotsList[generatingFunction]/valueAtMagicPoint,
            fvPatchField<scalar>::calculatedType()
        )
    );
    MagicCellsID.append (CellID) ;

    /******************************* EIM MAIN LOOP ********************************/

    scalar max_relative_error=desidered_error+1;
    scalar max_absolute_error ;
    InterpolationCoeffList.setSize(scalarSnapshotsList.size());
    scalar snapshotsMeasure = SMALL;
    label coeffIndex=0;


    while ( (max_relative_error > desidered_error) && (iterIndex <= maxBasis))
    {

        /* * * * update loop index* * * */

        ++iterIndex;

        /* * * * reset temp parameters* * * */
        
        max_absolute_error=SMALL;
        max_relative_error=SMALL;
        valueAtMagicPoint =SMALL;

        volScalarField* tmp_residual_Ptr =new GeometricField<scalar, fvPatchField, volMesh>
        (
            scalarSnapshotsList[0]
        );

        coeffIndex=iterIndex-2;


        /* * * * update interpolation Matrix* * * */

        updateBmatrix (Bmatrix, MagicFunctions, MagicCellsID);

        /* * * * get new residual field* * * */

        forAll (scalarSnapshotsList, snapI)
        {
            /* * * * find new EIM coefficient for snapshot snapI * * * */

            InterpolationCoeffList[snapI].setSize(coeffIndex+1);
            InterpolationCoeffList[snapI][coeffIndex] = scalarSnapshotsList[snapI][(MagicCellsID[coeffIndex])];


            if ((coeffIndex)!=0)  // fist iteration alpha1 = b1/B11 B11=1
            {
                for (label ii=0; ii< MagicCellsID.size()-1; ii++)
                {
                    InterpolationCoeffList[snapI][coeffIndex]-= (InterpolationCoeffList[snapI][ii]*Bmatrix(coeffIndex, ii));

                }
            }


            GeometricField<scalar, fvPatchField, volMesh> interpolant
            (
                IOobject
                (
                    scalarSnapshotsList[snapI].name() + "EIMInterpolation",
                    scalarSnapshotsList[snapI].time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<scalar>
                (
                    "zero",
                    scalarSnapshotsList[snapI].dimensions(),
                    pTraits<scalar>::zero
                )
            );

            for (label mfI=0; mfI< InterpolationCoeffList[snapI].size(); mfI++)
            {

                interpolant+=
                    InterpolationCoeffList[snapI][mfI]*MagicFunctions[mfI];
            }

            GeometricField<scalar, fvPatchField, volMesh> residual
            (
                (scalarSnapshotsList[snapI]-interpolant)
            );


            scalar absolute_interpolation_error =MOR::L_infinity_norm(residual);
            snapshotsMeasure = MOR::L_infinity_norm(scalarSnapshotsList[snapI]);
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
		generatingFunction = snapI;
            }

            if (max_relative_error< relative_interpolation_error)
            {
                max_relative_error=relative_interpolation_error;
            }

        }

        idxGeneratingFunction.append(generatingFunction);

        if ((max_relative_error > desidered_error) && (iterIndex <= maxBasis))
        {
            /* * * *find magic point * * * */

            forAll(*tmp_residual_Ptr, CellI)
            {
                if (mag(valueAtMagicPoint)< mag((*tmp_residual_Ptr)[CellI]))
                {
                    valueAtMagicPoint=(*tmp_residual_Ptr)[CellI];
                    CellID=CellI;
                }
            }


            /* * * * normalize MagicFunction * * * */

            *tmp_residual_Ptr/=valueAtMagicPoint;

            /* * * * update Magic Functions & Magic Points lists* * * */

            MagicFunctions.append (tmp_residual_Ptr);
            MagicCellsID.append (CellID) ;

        }

        if (iterIndex > maxBasis)
        {
            delete tmp_residual_Ptr;
        }

        /* * * * update errorLists* * * */

        absoluteErrorList.append(max_absolute_error);
        relativeErrorList.append(max_relative_error);

        /* * * * print info * * * */

        Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
        Info<<"maximum absolute interpolation error in Linfinty norm:     " <<  max_absolute_error << endl;
        Info<<"maximum relative interpolation error in Linfinty norm:     " <<  max_relative_error <<endl;
    }

    /*************************** Print Basis & files******************************/

    /*reset IO and write MagicFunction*/

    for(label mfI=0; mfI<MagicFunctions.size(); mfI++)
    {
        //runTime.setTime(mfI, mfI);

        MagicFunctions.set
        (   mfI,
            new volScalarField
            (   IOobject
                (
                    fieldName+"MagicFunction"+name(mfI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                *MagicFunctions(mfI),
                fvPatchField<scalar>::calculatedType()
            )
        );

        Info <<"\nWriting "<< MagicFunctions[mfI].name()<<" in folder: "<<runTime.timeName() << endl;
        MagicFunctions(mfI)->write();
    }

    /* create MagicPointField and write in "constant" folder */

    pointField MagicPointsField (MagicCellsID.size());

    forAll(MagicCellsID, CellI)
    {
        MagicPointsField[CellI]= mesh.C()[MagicCellsID[CellI]];
    }

    pointIOField MagicPointsIO
    (
        IOobject
        (
            fieldName+"MagicPoints",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        MagicPointsField
    );


    MagicPointsIO.writeObject(IOstream::ASCII, IOstream::currentVersion, IOstream::UNCOMPRESSED, true);

    Info <<"\nWriting "<< fieldName+"MagicPoints"<<" in folder: "<<runTime.constant()<< endl;

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
