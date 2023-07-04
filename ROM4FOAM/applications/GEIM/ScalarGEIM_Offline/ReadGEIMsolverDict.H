struct GEIMparameters
{
    List<fileName> folders_list;
    word fieldName;
    label  MaxSensorsNumber;
    scalar SensorsVariance;
    pointField sensorsPositions;
    scalar error;
};

void readSensorPositions
(
    const IOdictionary& GEIMsolverDict,
    pointField& sensorsPositions
)
{
    pointField tmp_sensorsPositions (
        GEIMsolverDict.subDict("Offline_parameters").lookup("SensorsPositions")
    );
    sensorsPositions = tmp_sensorsPositions;
};

// Main function
GEIMparameters getGEIMOfflineParameters ( Foam::argList args )
{
    Foam::Time runTime(Foam::Time::controlDictName, args);
    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );


    IOdictionary GEIMsolverDict
    (
        IOobject
        (
            "GEIMsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read entries from GEIM dictionary
    const entry* existField = GEIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = GEIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existsensorsNumber= GEIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("MaxSensorsNumber", false, true);
    const entry* existsensorsPositions= GEIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("SensorsPositions", false, true);
    const entry* existsensorsVariance= GEIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("SensorsVariance", false, true);
    const entry* existerror = GEIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("error", false, true);

    // Initiate variable from GEIMsolverDict
    if (!(existsensorsNumber) | !(existsensorsVariance) )
    {
        Info << "Error: you must define 'MaxSensorsNumber' and 'SensorsVariance'" << endl;
        abort();
    }

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' " << endl;
        abort();
    }

    if (!(existerror))
    {
        Info << "Error: 'error' entry not found. Must define 'error' " << endl;
        abort();
    }

    if (!existFoldersLists)
    {
        Info << "Error: 'foldersList' entry not found. Must define the list of folders containing the snapshots" << endl;
        abort();
    }

    word fieldName
    (
        GEIMsolverDict.subDict("Offline_parameters").lookup("field")
    );


    pointField sensorsPositions (mesh.C().field());

    if (!(existsensorsPositions))
    {
        Info << "No sensors location selected. As many sensors as Cells will be created  " << endl;
    }
    else
    {
        readSensorPositions(GEIMsolverDict,sensorsPositions);
        Info << "Reading sensor locations\n" << sensorsPositions << endl;
        for(label ii=0; ii< sensorsPositions.size(); ii ++)
        {
            for(label jj=ii+1; jj< sensorsPositions.size(); jj ++)
            {
                if(sensorsPositions[ii]==sensorsPositions[jj])
                {
                    Info << "Error: cannot have two or more sensors at the same location: " << sensorsPositions[ii] << endl;
                    abort();
                }
            }
        }
    }

    label MaxSensorsNumber =
        readInt
        (
            GEIMsolverDict.subDict("Offline_parameters").lookup("MaxSensorsNumber")
        );
    Info << "max number of sensors: " << MaxSensorsNumber << endl;

    if (MaxSensorsNumber > sensorsPositions.size())
    {
        Info << "ERROR! MaxSensorsNumber cannot be greater than the number of sensors available " << nl
             << "new MaxSensorsNumber is set = " << sensorsPositions.size() << endl;
        MaxSensorsNumber=sensorsPositions.size();
    }

    scalar SensorsVariance =
        readScalar
        (
            GEIMsolverDict.subDict("Offline_parameters").lookup("SensorsVariance")
        );
    Info << "sensors variance: " << SensorsVariance << endl;

    scalar error =
        readScalar
        (
            GEIMsolverDict.subDict("Offline_parameters").lookup("error")
        );

    Info << "Seeking error: " << error << endl;

    List<fileName> folders_list
    (
        GEIMsolverDict.subDict("Offline_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;


    GEIMparameters results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.MaxSensorsNumber= MaxSensorsNumber;
    results.sensorsPositions=sensorsPositions;
    results.SensorsVariance=SensorsVariance;
    results.error=error;

    return results;
}

