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
    ScalarPBDW_Online

Author
	Stefano Riva

\*---------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"

#include "regionProperties.H"
#include "turbulentFluidThermoModel.H" // necessary for the BC used in fluid-solid cases
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "getPBDWOnlineParameters.H"

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
        source[zz]= MOR::projection(field_to_interpolate, MagicSensors[zz]);
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
	    "reg",
	    "value",
	    "Specify if the regularization should be performed, with the proper value for the paramter xi" 
	);
argList::addBoolOption
    (
        "writeEstimate",
        "If activate the estimation will be written in the folder"
    );
argList::addOption
    (
        "noise",
        "value", 
        "If activate random noise is introduced with std deviation = <value>"
    );
#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"


    //get PBDW_Online parameters from PBDWsolverDict

    PBDWOnlineParameters PBDW_parameters = getPBDWOnlineParameters(args);

    word fieldName = PBDW_parameters.fieldName ;

    List<fileName> foldersList = PBDW_parameters.folders_list ;

    scalar BasisNumber= PBDW_parameters.BasisNumber;

    scalar MaxSensors = PBDW_parameters.MaxSensors;
    
    fileName sensorsFolder = PBDW_parameters.sensorsFolder;

    scalar noise_std;

const bool noisyData = args.optionFound("noise");
    if (noisyData)
    {
        args.optionLookup("noise")() >> noise_std; 
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> NormalGenerator{0,noise_std};

    /***********************************************************************************/
    if ( BasisNumber > MaxSensors )
    {
        Info<<"WARNING: "<<endl;
        Info<<"The dimension of the reduced space N should not be larger than the number of observables M!"<<endl;
        Info<<"The solution may be unstable, and the problem is not well posed!"<<endl;
        BasisNumber = MaxSensors - 1;
        Info<<"N is set to be equal to M-1 ("<<BasisNumber<<")"<<endl;
    }
    PtrList<volScalarField> BasisFunctions (BasisNumber);
    PtrList<volScalarField> BasisSensors (MaxSensors);

    scalar regularizingParameter = 0.0; // xi in the PBDW
    scalarSquareMatrix A(MaxSensors);
    RectangularMatrix<doubleScalar> K(MaxSensors, BasisNumber);
    
    bool BasisFunctionLoaded = false;

    while (BasisFunctionLoaded == false)
    {

        Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"
        Info<<"\nLoading the basis functions"<<endl;
        // get Magic Functions
        for(label mfI=0; mfI< BasisNumber; mfI ++)
        {
            BasisFunctions.set
            (
                mfI,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName+"PBDW_BasisFunction"+name(mfI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }
        BasisFunctionLoaded = true;
    }
    // get Magic Sensors
    Info<<"Loading the magic sensors"<<endl;
    bool MagicQuantitiesLoaded = false;

    while (MagicQuantitiesLoaded == false)
    {
        chDir(args.rootPath()/sensorsFolder);
    
    Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), sensorsFolder);
#include "CreateMesh.H"
        for(label msI=0; msI< MaxSensors; msI ++)
        {
            BasisSensors.set
            (
                msI,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName+"PBDW_BasisSensor"+name(msI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        Info<<"Generating the matrices"<<endl;

        for (label ii = 0; ii < MaxSensors; ++ii)
        {
            for (label jj = 0; jj < MaxSensors; ++jj)
            {
                A[ii][jj] = MOR::projection(BasisSensors[ii], BasisSensors[jj]);
            }

            for (label jj = 0; jj < BasisNumber; ++jj)
            {
                K[ii][jj] = MOR::projection(BasisSensors[ii], BasisFunctions[jj]);
            }
        }
    
        MagicQuantitiesLoaded = true;
    }

const bool writeEstimate = args.optionFound("writeEstimate");
bool writingPhase = false; // if true the interpolant will be written
if (writeEstimate)
    {
       writingPhase  = true;
       Info<<"The state estimations will be written in the folders"<<endl;
    }

    if (noisyData)
    {
        Info<<"\nThe data are polluted by random noise with std deviation = "<<noise_std<<endl;
        if ( args.optionFound("reg") )
        {
            args.optionLookup("reg")() >> regularizingParameter;
        }
        else
        {
            regularizingParameter = 0.0;
        }
        if (regularizingParameter != 0.0)
        {
            Info<<"The regularization is introduced with xi = "<<regularizingParameter<<endl;
        }
    }
// +++++++++++++++++++++++++ Assembling and solving the linear system ++++++++++++++++++//
/*
                        | xi*M*I+A     K |       | eta_M |       | y_obs |       
                        |                |   *   |       |   =   |       |
                        | K^T          0 |       |  z_N  |       |   0   |
where M is the observables used, I is the identity matrix.
*/

    List<scalar> maxAbsL2ErrorList  (MaxSensors - BasisNumber + 1 );
    List<scalar> averageAbsL2ErrorList ( MaxSensors - BasisNumber + 1 );
    List<scalar> maxRelL2ErrorList  ( MaxSensors - BasisNumber + 1 );
    List<scalar> averageRelL2ErrorList ( MaxSensors - BasisNumber + 1 );

    scalarRectangularMatrix relL2ErrorList;
    scalarRectangularMatrix absL2ErrorList;
    scalar rows = 0;
    
    dimensionedScalar dummyVariable 
    ( 
        "dummy", 
        dimensionSet(0,3,0,0,0,0,0), 
        scalar(1.0) 
    ); 
    
forAll (foldersList, folderI)
{
    chDir(args.rootPath()/foldersList[folderI]);
    Info <<"\nReconstructing fields in folder " <<foldersList[folderI]<< endl;
    Foam::Time runTime(args.rootPath(), foldersList[folderI]);

#include "CreateMesh.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "  Time = " << runTime.timeName() << endl;
        ++rows;
        relL2ErrorList.setSize (rows, MaxSensors - BasisNumber + 1);
        absL2ErrorList.setSize (rows, MaxSensors - BasisNumber + 1);

        for(label kk=BasisNumber; kk<=MaxSensors; ++kk)
        {
            volScalarField FieldToEstimate
                (
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );

            PtrList<volScalarField> tmp_BasisSensors (kk);;
            for(label ll=0; ll<kk; ll++)
            {
                tmp_BasisSensors.set
                (
                    ll,
                    new volScalarField
                    (
                        BasisSensors[ll]
                    )
                );
            }

            scalarField source;
            assembleSourceTerm(source, FieldToEstimate, tmp_BasisSensors);

            scalarSquareMatrix systemMatrix( kk + BasisNumber );
            for ( label ii = 0; ii < kk+BasisNumber; ++ii )
            {
                for ( label jj = 0; jj < kk+BasisNumber; ++jj )
                {
                    if ( ii < kk )
                    {
                        if ( jj < kk) // top left
                        {
                            systemMatrix[ii][jj] = A[ii][jj];
                            if (ii == jj)
                            {
                                systemMatrix[ii][jj] += (regularizingParameter * kk);
                            }
                        }
                        else // top right
                        {
                            label tmpLabel = jj - kk;
                            systemMatrix[ii][jj] = K[ii][tmpLabel];
                        }
                    }
                    else 
                    {
                        if ( jj < kk )
                        {
                            label tmpLabel = ii - kk;
                            systemMatrix[ii][jj] = K.T()[tmpLabel][jj];
                        }
                        else 
                        {
                            systemMatrix[ii][jj] = 0.0;
                        }
                    }
                }
            }

            scalarField rhs   ( kk+BasisNumber );

            for ( label ii = 0; ii < kk+BasisNumber; ++ii )
            {
                if ( ii < kk )
                {
                    rhs[ii] = source[ii]; // y_obs
                    if (noisyData)
                    {
                        rhs[ii] += NormalGenerator(gen);
                    }
                }
                else
                {
                    rhs[ii] = 0.0;
                }
            }

            simpleMatrix<scalar> LinearSystem
            (
                systemMatrix,
                rhs
            );
    
            scalarField coeffs = LinearSystem.solve();

            volScalarField  stateEstimation
            (
                IOobject
                (
                    fieldName + "PBDW_estimation",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                FieldToEstimate.mesh(),
                dimensioned<scalar>
                (
                    "zero",
                    FieldToEstimate.dimensions(),
                    pTraits<scalar>::zero
                )
            );
            dimensionedScalar dummyVariableField
            ( 
                "dummyField", 
                FieldToEstimate.dimensions(),
                scalar(1.0) 
            ); 
            forAll (coeffs, msI)
            {
                if (msI < kk)
                {
                    stateEstimation += coeffs[msI] * BasisSensors[msI] * (dummyVariableField * dummyVariable);
                }
                else
                {
                    label tmpLabel = msI-kk;
                    stateEstimation += coeffs[msI] * BasisFunctions[tmpLabel];
                }
            }

            volScalarField residual
            (
                IOobject
                (
                    fieldName + "PBDW_residual",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                FieldToEstimate-stateEstimation
            );

            scalar sumFieldError =
                MOR::L2norm
                (
                    residual
                );
            
            scalar measure =
                MOR::L2norm(FieldToEstimate);

            if (measure< SMALL)
            {
                measure = SMALL;
            }

            label tmpLabel = kk - BasisNumber;
            absL2ErrorList[rows-1][tmpLabel] = sumFieldError ;
            relL2ErrorList[rows-1][tmpLabel] = sumFieldError/(measure) ;

            if ( (writingPhase == true) & (kk == BasisNumber) )
            {
                stateEstimation.write();
                residual.write();
            }

            }

        }
        
    }

Info<<"\n\n"<<endl;

for(label kk=0; kk <= MaxSensors - BasisNumber ; ++kk)
{
    scalarField tmpAbs ( rows );
    scalarField tmpRel ( rows );
    for (label ii = 0; ii < rows; ++ii)
    {
        tmpRel[ii] = relL2ErrorList[ii][kk];
        tmpAbs[ii] = absL2ErrorList[ii][kk];
    }
    maxRelL2ErrorList[kk]=max(tmpRel);
    averageRelL2ErrorList[kk]=average(tmpRel);
    maxAbsL2ErrorList[kk]=max(tmpAbs);
    averageAbsL2ErrorList[kk]=average(tmpAbs);
    Info<<"\nSensors =  "<<kk+BasisNumber<<endl;
    Info << "Ave L2 relative reconstruction error with : " <<average(tmpRel)<<endl;
    Info << "Ave L2 absolute reconstruction error with : " <<average(tmpAbs)<<endl;
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
