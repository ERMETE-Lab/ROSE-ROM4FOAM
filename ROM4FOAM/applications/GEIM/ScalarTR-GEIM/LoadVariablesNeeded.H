
PtrList<volScalarField> MagicFunctions (msNumber);
PtrList<volScalarField> MagicSensors (msNumber);
scalarSquareMatrix BMatrix (msNumber);
scalarField meanVector;
scalarField stdVector;

bool VariablesLoaded = false;


while (VariablesLoaded == false)
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

    /********************** load mean and std values *************************/

    scalarField tmp_meanVector(IFstream(args.rootPath()/args.caseName()/fieldName+"_GEIM_Offline_files/GEIM_Offline_Coeffs_average_values.txt")());
    scalarField temp_stdVector(IFstream(args.rootPath()/args.caseName()/fieldName+"_GEIM_Offline_files/GEIM_Offline_Coeffs_std.txt")());

    meanVector=tmp_meanVector;
    stdVector=temp_stdVector;

    /*********************************************************************************/

    VariablesLoaded = true;
}
