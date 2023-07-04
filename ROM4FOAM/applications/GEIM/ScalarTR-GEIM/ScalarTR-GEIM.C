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
    ScalarTR-GEIM

Author
	Stefano Riva and Carolina Introini
	
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <random>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "simpleMatrix.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"
#include "EigenSolver.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "ReadGEIMsolverDict.H"

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

#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"


    //get GEIM_Online parameters

    GEIMOnlineParameters GEIM_parameters = getGEIMOnlineParameters(args);

    word fieldName = GEIM_parameters.fieldName ;

    List<fileName> foldersList = GEIM_parameters.folders_list ;

    label msNumber= GEIM_parameters.msNumber;

    scalar noise_std= GEIM_parameters.noise_std;

    label NumberRepeatedExperiments= GEIM_parameters.NumberRepeatedExperiments;

    std::random_device rd;

    std::mt19937 gen(rd());

    std::normal_distribution<> NormalGenerator{0,noise_std};

    List<scalar> averageL2ErrorList (msNumber);

    /*********************** load magic functions, sensors, interpolation matrix etc.. ***************************/

#include "LoadVariablesNeeded.H"

    // Noise Generator


    // repeat field reconstructions for each additional basis; kk number of basis 
    
    for(label kk=1; kk<=(msNumber); ++kk)
    {

        PtrList<volScalarField> tmp_MagicSensors (kk);;

        for(label ll=0; ll<kk; ll++)
        {
            tmp_MagicSensors.set
            (
                ll,
                new volScalarField
                (
                    MagicSensors[ll]
                )
            );
        }

        scalarSquareMatrix tmp_BMatrix (kk);

        for(label ll=0; ll<(kk); ll++)
        {
            for (label ii=0; ii<kk; ii++)
            {
                tmp_BMatrix[ll][ii]=BMatrix[ll][ii];
            }

        }

        scalarSquareMatrix BmatrixTBmatrix (tmp_BMatrix.T()*tmp_BMatrix);
        
        // Generate Tikhonov regularization matrix

        scalarSquareMatrix TR_Matrix (BmatrixTBmatrix.n(),zero());


        for(label ll=0; ll<(kk); ll++)
        {
            TR_Matrix(ll,ll)= 1.0/(std::pow(stdVector[ll],2));
        }


        List<scalar> L2ErrorList;

        forAll (foldersList, folderI)
        {

            chDir(args.rootPath()/foldersList[folderI]);

            Foam::Time runTime(args.rootPath(), foldersList[folderI]);

#include "CreateMesh.H"

            // Get times list
            instantList timeDirs = timeSelector::select0(runTime, args);

			
            for (label noiseI=0; noiseI<NumberRepeatedExperiments; noiseI++) // repeat Experiment
            {

                forAll(timeDirs, timeI)
                {
                    runTime.setTime(timeDirs[timeI], timeI);

                    volScalarField FieldToInterpolate
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

                    scalarField sourceWithNoise;
                    assembleSourceTerm(sourceWithNoise, FieldToInterpolate, tmp_MagicSensors);

                    forAll (sourceWithNoise, sensI)
                    {
                        sourceWithNoise[sensI]+= NormalGenerator(gen);
                    }

                    scalarField StabilisedSource(kk,zero());

                    for(label ii=0; ii<StabilisedSource.size() ; ii++)
                    {
                        for(label jj=0; jj<sourceWithNoise.size() ; jj++)
                        {
                            StabilisedSource[ii]+= tmp_BMatrix.T()(ii,jj)*sourceWithNoise[jj]+(std::pow(noise_std,2)*TR_Matrix(ii,jj))*meanVector[ii];
                        }
                    }


                    simpleMatrix<scalar> minimisation_problem
                    (
                        (BmatrixTBmatrix+(std::pow(noise_std,2)*TR_Matrix)),
                        StabilisedSource
                    );

                    scalarField coeffs = minimisation_problem.solve();

                    volScalarField  interpolant
                    (
                        IOobject
                        (
                            fieldName + "_TR-GEIMinterpolant",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        FieldToInterpolate.mesh(),
                        dimensioned<scalar>
                        (
                            "zero",
                            FieldToInterpolate.dimensions(),
                            pTraits<scalar>::zero
                        )
                    );

                    for (label ii=0; ii< coeffs.size(); ii++)
                    {
                        interpolant +=
                            coeffs[ii]*(MagicFunctions[ii]);
                    }

                    volScalarField residual
                    (
                        IOobject
                        (
                            fieldName + "_TR-GEIMresidual",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        FieldToInterpolate-interpolant
                    );

                    scalar interpolationdError =
                        MOR::L2norm
                        (
                            residual
                        );

                    scalar measure =
                        MOR::L2norm(FieldToInterpolate);

                    if (measure< SMALL)
                    {
                        measure = SMALL;
                    }

                    L2ErrorList.append(interpolationdError/(measure));

                    if ((kk==msNumber) & (noiseI==NumberRepeatedExperiments-1))
                    {
                        interpolant.write();
                        residual.write();
                    }

                }

            }

        }


        Info << "average L2 relative reconstruction error with " << kk << " basis: " << average(L2ErrorList)<< endl;

        averageL2ErrorList[kk-1]=average(L2ErrorList);


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
