bool EndSensorsGeneration= false;
Info<<"Generating the sensors for GEIM"<<endl;
while (EndSensorsGeneration == false)
{
    Foam::Time runTime(args.rootPath(), args.caseName());

// #include "CreateMesh.H"
    scalar maxSize = 12e3;
if ( sensorsPositions.size() > maxSize )
{
    Info<<"Warning: the mesh might too large (too much RAM occupied)"<<endl;
    scalar newSize = std::floor( sensorsPositions.size() / maxSize ) ;
    Info<<"Max available mesh size is "<<maxSize<<", hence the position will be sampled every "<<newSize<<endl;
    sensorsPositions.setSize (maxSize);
    pointField TMPsensorsPositions (mesh.C().field());
    label idxSens = 0;
    forAll(sensorsPositions, sensI)
    {
        sensorsPositions[sensI] = TMPsensorsPositions[idxSens];
        idxSens = idxSens + newSize;
    }
}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
sensors.setSize( sensorsPositions.size() );
tensor varianceTensor (1/(2*SensorsVariance), 0, 0, 0,1/(2*SensorsVariance),0,0,0,1/(2*SensorsVariance));

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
    // Info<<"The "<<sensI<<" sensor has been generated"<<endl; // this is used to check if there is any problem for the RAM
}

EndSensorsGeneration = true;
Info<< "End of sensors generation" << endl;
}
