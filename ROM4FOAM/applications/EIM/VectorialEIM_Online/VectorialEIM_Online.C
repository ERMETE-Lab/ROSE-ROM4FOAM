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
    VectorialEIM_Online

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
#include "turbulentFluidThermoModel.H" // necessary for the multi-region BC
// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //
#include "ReadEIMsolverDict.H"


void assembleSourceTerm
(
    scalarField& source,
    const PtrList<GeometricField<vector, fvPatchField, volMesh>> & MagicFunctions,
    const GeometricField<vector, fvPatchField, volMesh>& field_to_interpolate,
    const labelList& magic_cellsID
)
{
    label mc_size = magic_cellsID.size();
    source.setSize(mc_size);
    for (label zz=0; zz< mc_size; ++zz)
    {
        source[zz]= (field_to_interpolate[magic_cellsID[zz]])&(MagicFunctions[zz][(magic_cellsID[zz])]);
    }
}

int main(int argc, char *argv[])
{
	argList::addOption
		(
			"region",
			"name",
			"Specify the mesh region"
		);
#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"

    //get EIM_Online parameters

    EIMOnlineParameters EIM_parameters = getEIMOnlineParameters(args);

    word fieldName = EIM_parameters.fieldName ;

    List<fileName> foldersList = EIM_parameters.folders_list ;

    label mfNumber= EIM_parameters.mfNumber;


    /***********************************************************************************/

    PtrList<volVectorField> MagicFunctions (mfNumber);
    labelList MagicCellsID (mfNumber);
    scalarSquareMatrix BMatrix (mfNumber);

    bool MagicQuantitiesLoaded = false;

    while (MagicQuantitiesLoaded == false)
    {
        Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"

        // load Magic Functions
        runTime.setTime(0, 0);
        for(label mfI=0; mfI< mfNumber; mfI ++)
        {
            MagicFunctions.set
            (
                mfI,
                new volVectorField
                (
                    IOobject
                    (
                        fieldName+"MagicFunction"+name(mfI),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        // load Magic Points
        pointIOField MagicPoints(
            IOobject
            (
                fieldName+"MagicPoints",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ
            )
        );

        forAll (MagicCellsID, mpI)
        {
            MagicCellsID[mpI]=mesh.findCell(MagicPoints[mpI]);
        }


        for (label ii=0; ii< MagicFunctions.size(); ++ii)
        {
            for(label jj=0; jj< MagicFunctions.size(); ++jj)
            {
                if (jj<=ii)
                {
                    BMatrix(ii,jj)=(MagicFunctions[jj][(MagicCellsID[ii])])&(MagicFunctions[ii][(MagicCellsID[ii])]);
                }
                else
                {
                    BMatrix(ii,jj)=0;
                }
            }
        }

        MagicQuantitiesLoaded = true;
    }


    /***********************************************************************************/
    List<scalar> maxL2ErrorList  (mfNumber);
    List<scalar> averageL2ErrorList (mfNumber);

    // repeat field reconstructions for each additional basis; kk number of basis 

    Info <<"Reconstructing fields in folders: " << nl
         <<foldersList<< endl;

    for(label kk=1; kk<=mfNumber; ++kk)
    {

        labelList tmp_MagicCellsID (kk);

        for(label ll=0; ll<kk; ll++)
        {
            tmp_MagicCellsID[ll]=MagicCellsID[ll];
        }

        scalarSquareMatrix tmp_BMatrix (kk);

        for(label ll=0; ll<kk; ll++)
        {
            for (label ii=0; ii<kk; ii++)
            {
                tmp_BMatrix[ll][ii]=BMatrix[ll][ii];
            }

        }

        List<scalar> L2ErrorList;
        scalar error = 0.0;

        forAll (foldersList, folderI)
        {
            chDir(args.rootPath()/foldersList[folderI]);

            Foam::Time runTime(args.rootPath(), foldersList[folderI]);


#include "CreateMesh.H"


            // Get times list
            instantList timeDirs = timeSelector::select0(runTime, args);

            forAll(timeDirs, timeI)
            {
                runTime.setTime(timeDirs[timeI], timeI);

                volVectorField FieldToInterpolate
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

                scalarField source;
                assembleSourceTerm(source,MagicFunctions,FieldToInterpolate,tmp_MagicCellsID);
                simpleMatrix<scalar> Interpolation_problem
                (
                    tmp_BMatrix,
                    source
                );

                scalarField coeffs = Interpolation_problem.solve();

                volVectorField  interpolant
                (
                    IOobject
                    (
                        fieldName + "EIMinterpolant",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    FieldToInterpolate.mesh(),
                    dimensioned<vector>
                    (
                        "zero",
                        FieldToInterpolate.dimensions(),
                        pTraits<vector>::zero
                    )
                );


                forAll (coeffs, mfI)
                {
                    interpolant +=
                        coeffs[mfI]*(*MagicFunctions(mfI));
                }


                volVectorField residual
                (
                    IOobject
                    (
                        fieldName + "EIMresidual",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    FieldToInterpolate-interpolant
                );

                scalar sumFieldError =
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

                L2ErrorList.append(sumFieldError/(measure));
                if (error < sumFieldError/(measure) )
                {
                    error=sumFieldError/(measure);

                }
                
                if (kk==mfNumber)
                {
                    interpolant.write();
                    residual.write();
                }

            }

        }

        Info << "max L2 relative reconstruction error with " << kk << " basis: " <<max(L2ErrorList)<<endl;

        maxL2ErrorList[kk-1]=max(L2ErrorList);
        averageL2ErrorList[kk-1]=average(L2ErrorList);

    }

#include "printAllFiles.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //
