PtrList<volScalarField> sensors (sensorsPositions.size());
bool EndSensorsGeneration= false;

while (EndSensorsGeneration == false)
{
    Foam::Time runTime(args.rootPath(), args.caseName());
#include "CreateMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    tensor varianceTensor (1/(2*SensorsVariance), 0, 0, 0,1/(2*SensorsVariance),0,0,0,1/(2*SensorsVariance));

    Info << "Generating Sensors" << endl;

    //

    forAll (sensorsPositions, sensI)
    {
        sensors.set
        (
            sensI,
            new volScalarField
            (
                IOobject
                (
                    "sensor"+name(sensI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensioned<scalar>
                (
                    "zero",
                    dimless,
                    pTraits<scalar>::zero
                )
            )
        );

        forAll (sensors[sensI], CellI)
        {
            if ((Foam::exp(-(sensors[sensI].mesh().C()[CellI]-sensorsPositions[sensI])&(varianceTensor& (sensors[sensI].mesh().C()[CellI]-sensorsPositions[sensI]))) > SMALL))
            {
                sensors[sensI][CellI]= Foam::exp(-(sensors[sensI].mesh().C()[CellI]-sensorsPositions[sensI])&(varianceTensor& (sensors[sensI].mesh().C()[CellI]-sensorsPositions[sensI])));
            }
            else
            {
                sensors[sensI][CellI]=0;
            }
        }

        if ((fvc::domainIntegrate(sensors[sensI].internalField())).value() < SMALL)
		{
			Info << "Error: location of sensor " << sensI + 1 << " is likely outside the geometrical domain" << endl;
        	abort();
		}
		
		sensors[sensI]/= (fvc::domainIntegrate(sensors[sensI].internalField()));
    }

    EndSensorsGeneration = true;
    Info<< "End of sensors generation" << endl;
}
