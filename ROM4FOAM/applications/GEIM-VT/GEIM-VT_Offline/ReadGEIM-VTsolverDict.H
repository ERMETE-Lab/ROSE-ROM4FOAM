struct GEIM_VTparameters
{
    List<fileName> folders_list;
    label  MaxSensorsNumber;
    scalar SensorsVariance;
    pointField sensorsPositions;
    scalar error;
};

void readSensorPositions
(
    const IOdictionary& GEIM_VTsolverDict,
    pointField& sensorsPositions
)
{
    pointField tmp_sensorsPositions (
        GEIM_VTsolverDict.subDict("Offline_parameters").lookup("SensorsPositions")
    );
    sensorsPositions = tmp_sensorsPositions;
};

// Main function
GEIM_VTparameters getGEIMOfflineParameters ( Foam::argList args )
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


    IOdictionary GEIM_VTsolverDict
    (
        IOobject
        (
            "GEIM-VTsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read entries from  GEIM_VT dictionary
    const entry* existFoldersLists = GEIM_VTsolverDict.subDict("Offline_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existsensorsNumber= GEIM_VTsolverDict.subDict("Offline_parameters").lookupEntryPtr("MaxSensorsNumber", false, true);
    const entry* existsensorsPositions= GEIM_VTsolverDict.subDict("Offline_parameters").lookupEntryPtr("SensorsPositions", false, true);
    const entry* existsensorsVariance= GEIM_VTsolverDict.subDict("Offline_parameters").lookupEntryPtr("SensorsVariance", false, true);
    const entry* existerror = GEIM_VTsolverDict.subDict("Offline_parameters").lookupEntryPtr("error", false, true);

    // Initiate variable from GEIM_VTsolverDict
    if (!(existsensorsNumber) | !(existsensorsVariance) )
    {
        Info << "Error: you must define 'MaxSensorsNumber' and 'SensorsVariance'" << endl;
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

    pointField sensorsPositions (mesh.C().field());

    if (!(existsensorsPositions))
    {
        Info << "No sensors location selected. As many sensors as Cells will be created  " << endl;
    }
    else
    {
        readSensorPositions(GEIM_VTsolverDict,sensorsPositions);
        Info << "Reading sensor locations\n" << sensorsPositions << endl;
        for(label ii=0; ii< sensorsPositions.size(); ii ++)
        {
            for(label jj=ii+1; jj< sensorsPositions.size(); jj ++)
            {
                if(sensorsPositions[ii]==sensorsPositions[jj])
                {
                    Info << "Error: cannot have two or more sensors at same location: " << sensorsPositions[ii] << endl;
                    abort();
                }
            }
        }
    }

    label MaxSensorsNumber =
        readInt
        (
            GEIM_VTsolverDict.subDict("Offline_parameters").lookup("MaxSensorsNumber")
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
            GEIM_VTsolverDict.subDict("Offline_parameters").lookup("SensorsVariance")
        );
    Info << "sensors variance: " << SensorsVariance << endl;

    scalar error =
        readScalar
        (
            GEIM_VTsolverDict.subDict("Offline_parameters").lookup("error")
        );

    Info << "Seeking error: " << error << endl;

    List<fileName> folders_list
    (
        GEIM_VTsolverDict.subDict("Offline_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;


    GEIM_VTparameters results;
    results.folders_list=folders_list;
    results.MaxSensorsNumber= MaxSensorsNumber;
    results.MaxSensorsNumber= MaxSensorsNumber;
    results.SensorsVariance=SensorsVariance;
    results.error=error;
    results.sensorsPositions=sensorsPositions;

    return results;
}

