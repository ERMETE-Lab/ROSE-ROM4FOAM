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
    ScalarPBDW_Offline

Author
	Stefano Riva

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
#include "QRMatrix.H"
#include "SortableList.H"

#include "regionProperties.H"
#include "turbulentFluidThermoModel.H" // necessary for the BC used in fluid-solid cases

// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //
#include "ReadPBDWsolverDict.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//- update B matrix (interpolation matrix)
void updateBmatrix
(
    scalarSquareMatrix& BMatrix,
    const PtrList<GeometricField<scalar, fvPatchField, volMesh>>& Basis_functions,
    const PtrList<GeometricField<scalar, fvPatchField, volMesh>>& MagicSensors
)
{
    label mf_size = Basis_functions.size();
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
        BMatrix(mf_size-1,jj)=MOR::projection(Basis_functions[jj], MagicSensors[mf_size-1]);

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

argList::addOption
    (
        "algoRS",
        "name",
        "The algorithm to build the Reduced Space can be chosen (WeakGreedy, GEIM). If GEIM is selected the sensors will be chosen accordingly"
    );

#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"


    //get PBDW_Offline parameters from PBDWsolverDict

    PBDWparameters PBDW_parameters = getPBDWOfflineParameters(args);

    word fieldName = PBDW_parameters.fieldName ;

    List<fileName> foldersList = PBDW_parameters.folders_list ;

    scalar MaxSensorsNumber = PBDW_parameters.MaxSensorsNumber ;

    pointField sensorsPositions = PBDW_parameters.sensorsPositions;

    scalar SensorsVariance = PBDW_parameters.SensorsVariance;

    scalar BasisNumber = PBDW_parameters.BasisNumber;

    // Initialise List necessary to store the snapshots 

    PtrList<volScalarField> scalarSnapshotsList;

    // get ScalarSnapshots

    for (label folderI=0; folderI < foldersList.size(); ++ folderI)
    {
        chDir(args.rootPath()/foldersList[folderI]);

        Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), foldersList[folderI]);
        
#include "CreateMesh.H"

        Info << "\nReading snapshots in "<< foldersList[folderI]<< endl;

        // Get times list
        instantList timeDirs = timeSelector::select0(runTime, args);

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            Info<< "Time = " << runTime.timeName() << endl;

            if (mesh.readUpdate() != polyMesh::UNCHANGED)
            {

                FatalErrorInFunction
                        << "polyMesh has changed. PBDW can be performed only on unchanged polyMesh"
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
    
// initialize GEIM Sensors
PtrList<volScalarField> sensors ;

// Selection of the Algorithm

word algorithmRS = "WeakGreedy";

algorithmRS = args.optionLookupOrDefault("algoRS", algorithmRS);

Info<<"\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"<<endl;

if ( algorithmRS != "WeakGreedy" && algorithmRS != "GEIM")
{
    Info<<"Warning: the available algorithms to build the reduced space are WeakGreedy or GEIM only\n"
        <<"Default algorithm (WeakGreedy) has been selected"<<endl;
    algorithmRS = "WeakGreedy";
}
else
{
    Info<<"The Reduced Space will be built with algorithm = "<<algorithmRS<<endl;
}
// Create PBDW_(fieldName)_algorithmRS folder

#include "createFolderResultsandResetTime.H"
#include "CreateMesh.H"

if (algorithmRS == "GEIM")
{
    #include "generateSensors.H"
}

Info<<"\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<endl;

//+++++++++++++++++++++++++++++++++++ Construction of the reduced Space Z_N +++++++++++++++++++++++++++++++++++

    PtrList<volScalarField> BasisFunctions;
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
    
    Info<<"Starting the generation of the Reduced Space with "<<algorithmRS<<endl;
    if ( algorithmRS == "GEIM" )
    {
        if ( BasisNumber != MaxSensorsNumber )
        {
            BasisNumber = MaxSensorsNumber ;
            Info<<"Warning: the BasisNumber cannot be different from the MaxSensorsNumber with GEIM"<<endl;
        }
        #include "GEIM_Loop.H"
    }
    else if ( algorithmRS == "WeakGreedy" )
    {
        #include "WeakGreedy_Loop.H"
    }

Info<<"\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<endl;

//++++++++++++++++++ Selection of the sensors using SGREEDY +++++++++++++++++++++++++++++++++++

PtrList<volScalarField> BasisSensors;
labelList SensorsCellsID;
scalarField infSupConstant;

if ( algorithmRS != "GEIM" )
{
    bool SensorSelection = false;
    Info<<"The sensors are selected using the SGREEDY algorithm"<<endl;
    while ( SensorSelection == false )
    {
        #include "S_GREEDY.H"
        SensorSelection = true;
    }
}
else
{
    Info<<"The sensors are taken from GEIM algorithm"<<endl;
}

Info<<"\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<endl;

//+++++++++++++++++++++++++++++++++++ Print Basis, Sensors & files +++++++++++++++++++++++++++++++++++

    
scalarSquareMatrix checkMatrix ( BasisFunctions.size() );
for (label ii = 0; ii < BasisFunctions.size(); ++ii )
{
    for (label jj = 0; jj < BasisFunctions.size(); ++jj )
    {
        checkMatrix[ii][jj] = MOR::projection(BasisFunctions[ii], BasisFunctions[jj] );
    }
}

List<scalar> LebesgueConstantList;
if ( algorithmRS == "GEIM" )
{
    Info<<"\nComputing Lebesgue Constant"<<endl;
    //#include "computeLebesgueConstant.H"
}
    


Info <<"\nWriting BasisFunctions in folder: "<<runTime.timeName() << endl;
for(label mfI=0; mfI<BasisFunctions.size(); mfI++)
{
    BasisFunctions.set
    (   mfI,
        new volScalarField
        (   IOobject
            (
                fieldName+"PBDW_BasisFunction"+name(mfI),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            BasisFunctions[mfI],
            fvPatchField<scalar>::calculatedType()
        )
    );

    BasisFunctions(mfI)->write();
}

Info <<"\nWriting BasisSensors in folder: "<<runTime.timeName() << endl;

if ( algorithmRS == "GEIM" )
{
    forAll(MagicSensors, msI)
    {
        
        MagicSensors.set
        (   msI,
            new volScalarField
            (   IOobject
                (
                    fieldName+"PBDW_BasisSensor"+name(msI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                MagicSensors[msI],
                fvPatchField<scalar>::calculatedType()
            )
        );

        MagicSensors(msI)->write();
        
    }
    
}
else
{
    
    forAll(BasisSensors, msI)
    {
        BasisSensors(msI)->write();
    }
    
}

#include "printAllFiles.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    auto stop_s=std::clock();

    Info<< nl << "ExecutionTime = " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
