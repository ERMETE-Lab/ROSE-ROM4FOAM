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
    ScalarGEIM_Online

Author
	Stefano Riva and Carolina Introini

\*---------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <random>
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

#include "regionProperties.H"
#include "turbulentFluidThermoModel.H" // necessary for the BC used in fluid-solid cases
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "getGEIMOnlineParameters.H"

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
argList::addBoolOption
    (
        "writeInterpolant",
        "If activate the intepolant will be written in the folder"
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


    //get GEIM_Online parameters from GEIMsolverDict

    GEIMOnlineParameters GEIM_parameters = getGEIMOnlineParameters(args);

    word fieldName = GEIM_parameters.fieldName ;

    List<fileName> foldersList = GEIM_parameters.folders_list ;

    label msNumber= GEIM_parameters.msNumber;

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

    PtrList<volScalarField> MagicFunctions (msNumber);
    PtrList<volScalarField> MagicSensors (msNumber);
    scalarSquareMatrix BMatrix (msNumber);
    bool MagicQuantitiesLoaded = false;

    while (MagicQuantitiesLoaded == false)
    {
        Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"

        // get Magic Functions
        for(label mfI=0; mfI< msNumber; mfI ++)
        {
            MagicFunctions.set
            (
                mfI,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName+"GEIMMagicFunction"+name(mfI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        // get Magic Sensors
        for(label msI=0; msI< msNumber; msI ++)
        {
            MagicSensors.set
            (
                msI,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName+"MagicSensor"+name(msI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }


        for (label ii=0; ii< MagicFunctions.size(); ++ii)
        {
            for(label jj=0; jj< MagicFunctions.size(); ++jj)
            {
                if (jj<=ii)
                {
                    BMatrix(ii,jj)=MOR::projection(MagicFunctions[jj], MagicSensors[ii]);
                }
                else
                {
                    BMatrix(ii,jj)=0;
                }
            }
        }

        MagicQuantitiesLoaded = true;
    }

const bool writeInterpolant = args.optionFound("writeInterpolant");
bool writingPhase = false; // if true the interpolant will be written
if (writeInterpolant)
    {
       writingPhase  = true;
       Info<<"The interpolant will be written in the folders"<<endl;
    }
    
    /***********************************************************************************/

    List<scalar> maxAbsL2ErrorList  ( msNumber );
    List<scalar> averageAbsL2ErrorList ( msNumber );
    List<scalar> maxRelL2ErrorList  ( msNumber );
    List<scalar> averageRelL2ErrorList ( msNumber );

    scalarRectangularMatrix relL2ErrorList;
    scalarRectangularMatrix absL2ErrorList;
    scalar rows = 0;

    if (noisyData)
    {
        Info<<"\nThe data are polluted by random noise with std deviation = "<<noise_std<<endl;
    }

    // repeat field reconstructions for each additional basis; kk number of basis 
forAll (foldersList, folderI)
{
    chDir(args.rootPath()/foldersList[folderI]);

    Foam::Time runTime(args.rootPath(), foldersList[folderI]);
    Info <<"\nReconstructing fields in folder " <<foldersList[folderI]<< endl;
#include "CreateMesh.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "  Time = " << runTime.timeName() << endl;

        ++rows;
        relL2ErrorList.setSize (rows, msNumber);
        absL2ErrorList.setSize (rows, msNumber);

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
        for(label kk=1; kk<=msNumber; ++kk)
        {

            PtrList<volScalarField> tmp_MagicSensors (kk);

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

            for(label ll=0; ll<kk; ll++)
            {
                for (label ii=0; ii<kk; ii++)
                {
                    tmp_BMatrix[ll][ii]=BMatrix[ll][ii];
                }

            }

                scalarField source;
                assembleSourceTerm(source, FieldToInterpolate, tmp_MagicSensors);

                if (noisyData)
                {
                    forAll (source, idxSource)
                    {
                        source[idxSource] +=  NormalGenerator(gen);
                    }
                }
                
                simpleMatrix<scalar> Interpolation_problem
                (
                    tmp_BMatrix,
                    source
                );

                scalarField coeffs = Interpolation_problem.solve();

                volScalarField  interpolant
                (
                    IOobject
                    (
                        fieldName + "GEIMinterpolant",
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


                forAll (coeffs, msI)
                {
                    interpolant +=
                        coeffs[msI]*(*MagicFunctions(msI));

                }


                volScalarField residual
                (
                    IOobject
                    (
                        fieldName + "GEIMresidual",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    FieldToInterpolate-interpolant
                );

                scalar interpolationError =
                    MOR::L2norm
                    (
                        residual
                    );

                scalar measure =
                    MOR::L2norm(FieldToInterpolate);

                if (measure < SMALL)
                {
                    measure = SMALL;
                }
                absL2ErrorList[rows-1][kk-1] = interpolationError ;
                relL2ErrorList[rows-1][kk-1] = interpolationError/(measure) ;
                if ( (writingPhase == true) & (kk == msNumber) )
                {
                    interpolant.write();
                    residual.write();
                }

        }

    }

}


Info<<"\n\n"<<endl;

for(label kk=0; kk < msNumber; ++kk)
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
    Info<<"\nSensors =  "<<kk<<endl;
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
