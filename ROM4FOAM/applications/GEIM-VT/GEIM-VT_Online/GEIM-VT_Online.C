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
    GEIM-VT_Online

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
#include "simpleMatrix.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "ReadGEIM-VTsolverDict.H"

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

    GEIM_VTOnlineParameters GEIM_VT_parameters = getGEIM_VTOnlineParameters(args);

    List<fileName> foldersList = GEIM_VT_parameters.folders_list ;

    label msNumber= GEIM_VT_parameters.msNumber;


    /***********************************************************************************/

    PtrList<volScalarField> T_MagicFunctions (msNumber);
    PtrList<volScalarField> p_MagicFunctions (msNumber);
    PtrList<volVectorField> U_MagicFunctions (msNumber);
    PtrList<volScalarField> MagicSensors (msNumber);
    scalarField VolumeField;
    scalarSquareMatrix BMatrix (msNumber);
    bool MagicQuantitiesLoaded = false;

    while (MagicQuantitiesLoaded == false)
    {
        Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"

        // get Magic Functions
        for(label mfI=0; mfI< msNumber; mfI ++)
        {
            T_MagicFunctions.set
            (
                mfI,
                new volScalarField
                (
                    IOobject
                    (
                        "TGEIM-VTMagicFunction"+name(mfI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        for(label mfI=0; mfI< msNumber; mfI ++)
        {
            p_MagicFunctions.set
            (
                mfI,
                new volScalarField
                (
                    IOobject
                    (
                        "p_rghGEIM-VTMagicFunction"+name(mfI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        for(label mfI=0; mfI< msNumber; mfI ++)
        {
            U_MagicFunctions.set
            (
                mfI,
                new volVectorField
                (
                    IOobject
                    (
                        "UGEIM-VTMagicFunction"+name(mfI),
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
                        "TMagicSensor"+name(msI),
                        "0",
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

        for (label ii=0; ii< T_MagicFunctions.size(); ++ii)
        {
            for(label jj=0; jj< T_MagicFunctions.size(); ++jj)
            {
                if (jj<=ii)
                {
                    BMatrix(ii,jj)=MOR::projection(T_MagicFunctions[jj], MagicSensors[ii]);
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

    List<scalar> T_maxL2ErrorList  (msNumber);
    List<scalar> T_averageL2ErrorList (msNumber);

    List<scalar> p_maxL2ErrorList  (msNumber);
    List<scalar> p_averageL2ErrorList (msNumber);

    List<scalar> U_maxL2ErrorList  (msNumber);
    List<scalar> U_averageL2ErrorList (msNumber);

    // repeat field reconstructions for each additional basis; kk number of basis 

    for(label kk=1; kk<=msNumber; ++kk)
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

        for(label ll=0; ll<kk; ll++)
        {
            for (label ii=0; ii<kk; ii++)
            {
                tmp_BMatrix[ll][ii]=BMatrix[ll][ii];
            }

        }

        List<scalar> TL2ErrorList;
        List<scalar> p_rghL2ErrorList;
        List<scalar> UL2ErrorList;

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

                volScalarField T_ToInterpolate
                (
                    IOobject
                    (
                        "T",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );

                volScalarField p_ToInterpolate
                (
                    IOobject
                    (
                        "p_rgh",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );

                volVectorField U_ToInterpolate
                (
                    IOobject
                    (
                        "U",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );

                scalarField source;
                assembleSourceTerm(source, T_ToInterpolate, tmp_MagicSensors);
                simpleMatrix<scalar> Interpolation_problem
                (
                    tmp_BMatrix,
                    source
                );

                scalarField coeffs = Interpolation_problem.solve();

                volScalarField  T_interpolant
                (
                    IOobject
                    (
                        "T_GEIM-VTinterpolant",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    T_ToInterpolate.mesh(),
                    dimensioned<scalar>
                    (
                        "zero",
                        T_ToInterpolate.dimensions(),
                        pTraits<scalar>::zero
                    )
                );

                volScalarField  p_interpolant
                (
                    IOobject
                    (
                        "p_rgh_GEIM-VTinterpolant",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    p_ToInterpolate.mesh(),
                    dimensioned<scalar>
                    (
                        "zero",
                        p_ToInterpolate.dimensions(),
                        pTraits<scalar>::zero
                    )
                );

                volVectorField  U_interpolant
                (
                    IOobject
                    (
                        "U_GEIM-VTinterpolant",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    U_ToInterpolate.mesh(),
                    dimensioned<vector>
                    (
                        "zero",
                        U_ToInterpolate.dimensions(),
                        pTraits<vector>::zero
                    )
                );

                forAll (coeffs, msI)
                {
                    T_interpolant +=
                        coeffs[msI]*(T_MagicFunctions[msI]);
                    p_interpolant +=
                        coeffs[msI]*(p_MagicFunctions[msI]);
                    U_interpolant +=
                        coeffs[msI]*(U_MagicFunctions[msI]);
                }

                volScalarField T_residual
                (
                    IOobject
                    (
                        "T_GEIM-VTresidual",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    T_ToInterpolate-T_interpolant
                );

                volScalarField p_residual
                (
                    IOobject
                    (
                        "p_rgh_GEIM-VTresidual",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    p_ToInterpolate-p_interpolant
                );

                volVectorField U_residual
                (
                    IOobject
                    (
                        "U_GEIM-VTresidual",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    U_ToInterpolate-U_interpolant
                );

                scalar T_interpolationdError =
                    MOR::L2norm
                    (
                        T_residual
                    );

                scalar p_interpolationdError =
                    MOR::L2norm
                    (
                        p_residual
                    );

                scalar U_interpolationdError =
                    MOR::L2norm
                    (
                        U_residual
                    );

                scalar T_measure =
                    MOR::L2norm(T_ToInterpolate);
                if (T_measure< SMALL)
                {
                    T_measure = SMALL;
                }
                scalar p_measure =
                    MOR::L2norm(p_ToInterpolate);
                if (p_measure< SMALL)
                {
                    p_measure = SMALL;
                }
                scalar U_measure =
                    MOR::L2norm(U_ToInterpolate);
                if (U_measure< SMALL)
                {
                    U_measure = SMALL;
                }

                TL2ErrorList.append(T_interpolationdError/(T_measure));
                p_rghL2ErrorList.append(p_interpolationdError/(p_measure));
                UL2ErrorList.append(U_interpolationdError/(U_measure));

                if (kk==msNumber)
                {
                    T_interpolant.write();
                    U_interpolant.write();
                    p_interpolant.write();
                    T_residual.write();
                    U_residual.write();
                    p_residual.write();
                }

            }

        }
        Info << " Basis Number: " << kk <<endl;
        Info << "max T L2 reconstruction error with " << kk << " basis: " <<max(TL2ErrorList)<<endl;
        Info << "max p_rgh L2 reconstruction error with " << kk << " basis: " <<max(p_rghL2ErrorList)<<endl;
        Info << "max U L2 reconstruction error with " << kk << " basis: " <<max(UL2ErrorList)<<endl;

        T_maxL2ErrorList[kk-1]=max(TL2ErrorList);
        T_averageL2ErrorList[kk-1]=average(TL2ErrorList);

        p_maxL2ErrorList[kk-1]=max(p_rghL2ErrorList);
        p_averageL2ErrorList[kk-1]=average(p_rghL2ErrorList);

        U_maxL2ErrorList[kk-1]=max(UL2ErrorList);
        U_averageL2ErrorList[kk-1]=average(UL2ErrorList);
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
