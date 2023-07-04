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
    VectorialPOD_Online

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

    PtrList<volVectorField>  vectorPODOrthoNormalBasis (BasisNumber);

    bool PODBasisLoaded = false;
    Info<<"\nLoading POD Basis"<<endl;
    while (PODBasisLoaded == false)
    {
        Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"

        // load Magic Functions
        runTime.setTime(0, 0);
        for(label mfI=0; mfI< BasisNumber; mfI ++)
        {
            //runTime.setTime(mfI, mfI);
            vectorPODOrthoNormalBasis.set
            (
                mfI,
                new volVectorField
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

    List<scalar> maxL2ErrorList  (BasisNumber);
    List<scalar> averageL2ErrorList (BasisNumber);
    List<scalar> maxL2ErrorList_Energy  (BasisNumber);
    List<scalar> averageL2ErrorList_Energy (BasisNumber);
    scalarRectangularMatrix L2ErrorMatrix;
    scalarRectangularMatrix L2ErrorMatrix_Energy;
    label snapIdx=0;
    Info<<"Beginning Reconstruction"<<endl;
    // repeat field reconstructions for each additional basis; kk number of basis 

       forAll (foldersList, folderI)
        {
            chDir(args.rootPath()/foldersList[folderI]);

            Foam::Time runTime(args.rootPath(), foldersList[folderI]);

            //Info <<"Reconstructing fields in folder " <<foldersList[folderI]<< endl;

#include "CreateMesh.H"


            // Get times list
            instantList timeDirs = timeSelector::select0(runTime, args);

            forAll(timeDirs, timeI)
            {
                ++snapIdx;
                L2ErrorMatrix.setSize(snapIdx, BasisNumber);
                L2ErrorMatrix_Energy.setSize(snapIdx, BasisNumber);
                runTime.setTime(timeDirs[timeI], timeI);

                volVectorField FieldToReconstruct
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
                            vectorPODOrthoNormalBasis[kk-1]
                        );

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
                        FieldToReconstruct.dimensions(),
                        pTraits<vector>::zero
                    )
                );

                for (label baseI = 0; baseI < kk; baseI++)
                {
                    reconstruct +=
                        coeffs[baseI]*vectorPODOrthoNormalBasis[baseI];
                }


                volVectorField residual
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

                if (measure< SMALL)
                {
                    measure = SMALL;
                }

                L2ErrorMatrix[snapIdx-1][kk-1] = sumFieldError/(measure);
               
                volScalarField FieldToReconstructEnergy
                (
                    IOobject
                    (
                        fieldName+"Energy",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    FieldToReconstruct&FieldToReconstruct
                );
                
                volScalarField reconstructEnergy
                (
                    IOobject
                    (
                        fieldName + "PODreconstructEnergy",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    reconstruct&reconstruct
                );

                volScalarField residualEnergy
                (
                    IOobject
                    (
                        fieldName + "PODresidualEnergy",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    ( FieldToReconstructEnergy - reconstructEnergy )
                );
                
                scalar energyDifference;
                energyDifference = MOR::L2norm(
                            0.5 * residualEnergy
                );

                scalar EnergyMeasure = MOR::L2norm( 0.5 * FieldToReconstructEnergy );
                L2ErrorMatrix_Energy[snapIdx-1][kk-1] = energyDifference/(EnergyMeasure);

                if (kk==BasisNumber)
                {
                    reconstruct.write();
                    residual.write();
                }

            }
            Info<<"Completed folder: "<<foldersList[folderI]<<" and Time = "<< runTime.timeName() <<endl;
        }
    }
    
Info<<"\nRecostrunction Errors: "<<endl;
List<scalar> maxIndex;
 for(label kk=1; kk<=BasisNumber; ++kk)
    {
        List<scalar> tmp;
        List<scalar> tmpEnergy;
        for (label snapI = 0; snapI < L2ErrorMatrix.size() / BasisNumber; ++snapI)
        {
            tmp.append(L2ErrorMatrix[snapI][kk-1]);
            tmpEnergy.append(L2ErrorMatrix_Energy[snapI][kk-1]);
            
        }

Info << "max L2 relative reconstruction error with " << kk << " basis: " <<max(tmp)<<endl;

        maxL2ErrorList[kk-1]=max(tmp);
        averageL2ErrorList[kk-1]=average(tmp);
        maxL2ErrorList_Energy[kk-1]=max(tmpEnergy);
        averageL2ErrorList_Energy[kk-1]=average(tmpEnergy);


        // Selection of the snapshot that maximizes the error
        scalar maximum = SMALL;
        maxIndex.setSize(kk);
        for(label snapI = 0; snapI < tmp.size(); ++snapI)
        {
            if (tmp[snapI] > maximum)
            {
                maxIndex[kk-1] = snapI;
                maximum = tmp[snapI];
            }
        }
        // *************************************************** //

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
