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
    ScalarPOD_Online

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
#include "simpleMatrix.H"
#include "IOmanip.H"
#include "IOstream.H"

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

    //get POD_Online parameters

    PODOnlineParameters POD_parameters = getPODOnlineParameters(args);

    word fieldName = POD_parameters.fieldName ;

    List<fileName> foldersList = POD_parameters.folders_list ;

    label BasisNumber= POD_parameters.BasisNumber;


    /***********************************************************************************/

    PtrList<volScalarField>  scalarPODOrthoNormalBasis (BasisNumber);

    bool PODBasisLoaded = false;
    Info<<"\nLoading POD basis"<<endl;
    while (PODBasisLoaded == false)
    {
        Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"

        // load Magic Functions
        runTime.setTime(0, 0);
        for(label mfI=0; mfI< BasisNumber; mfI ++)
        {
            //runTime.setTime(mfI, mfI);
            scalarPODOrthoNormalBasis.set
            (
                mfI,
                new volScalarField
                (
                    IOobject
                    (
                        fieldName+"POD"+name(mfI),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        PODBasisLoaded = true;
    }


    /***********************************************************************************/
    
    List<scalar> maxAbsL2ErrorList  (BasisNumber);
    List<scalar> averageAbsL2ErrorList (BasisNumber);
    List<scalar> maxRelL2ErrorList  (BasisNumber);
    List<scalar> averageRelL2ErrorList (BasisNumber);

    scalarRectangularMatrix relL2ErrorList;
    scalarRectangularMatrix absL2ErrorList;
    label snapIdx=0;

Info<<"\nBeginning Reconstruction"<<endl;
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
                
                ++snapIdx;
                relL2ErrorList.setSize(snapIdx, BasisNumber);
                absL2ErrorList.setSize(snapIdx, BasisNumber);

                volScalarField FieldToReconstruct
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
            List<scalar> coeffs;
            for(label kk=1; kk<=BasisNumber; ++kk)
            {
                        coeffs.setSize(kk);
                        coeffs[kk-1] = MOR::projection
                        (
                            FieldToReconstruct,
                            scalarPODOrthoNormalBasis[kk-1]
                        );

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
                        FieldToReconstruct.dimensions(),
                        pTraits<scalar>::zero
                    )
                );

                for (label baseI = 0; baseI < kk; baseI++)
                {
                    reconstruct +=
                        coeffs[baseI]*scalarPODOrthoNormalBasis[baseI];
                }

                volScalarField residual
                (
                    IOobject
                    (
                        fieldName + "PODresidual",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    FieldToReconstruct-reconstruct
                );

		// Assembling L2 error
                scalar sumFieldError =
                    MOR::L2norm
                    (
                        residual
                    );

                scalar measure =
                    MOR::L2norm(FieldToReconstruct);

                if (measure < SMALL)
                {
                    measure = SMALL;
                }

                absL2ErrorList[snapIdx-1][kk-1] = sumFieldError ;
                 relL2ErrorList[snapIdx-1][kk-1] = sumFieldError/(measure) ;

                if (kk==BasisNumber)
                {
                    reconstruct.write();
                    residual.write();
                }

            }
        }

    }


Info<<"\n\n"<<endl;

for(label kk=1; kk<=BasisNumber; ++kk)
{
    scalarField tmpAbs ( snapIdx );
    scalarField tmpRel ( snapIdx );
    for (label ii = 0; ii < snapIdx; ++ii)
    {
        tmpRel[ii] = relL2ErrorList[ii][kk-1];
        tmpAbs[ii] = absL2ErrorList[ii][kk-1];
    }
    maxRelL2ErrorList[kk-1]=max(tmpRel);
    averageRelL2ErrorList[kk-1]=average(tmpRel);
    maxAbsL2ErrorList[kk-1]=max(tmpAbs);
    averageAbsL2ErrorList[kk-1]=average(tmpAbs);
    Info<<"\nBasisNumber =  #"<<kk<<endl;
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
