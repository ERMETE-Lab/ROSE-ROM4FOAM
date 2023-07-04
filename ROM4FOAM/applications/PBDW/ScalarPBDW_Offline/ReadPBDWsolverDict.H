struct PBDWparameters
{
    List<fileName> folders_list;
    word fieldName;
    scalar BasisNumber;
    scalar  MaxSensorsNumber;
    scalar SensorsVariance;
    pointField sensorsPositions;
};

void readSensorPositions
(
    const IOdictionary& PBDWsolverDict,
    pointField& sensorsPositions
)
{
    pointField tmp_sensorsPositions (
        PBDWsolverDict.subDict("Offline_parameters").lookup("SensorsPositions")
    );
    sensorsPositions = tmp_sensorsPositions;
};

// Main function
PBDWparameters getPBDWOfflineParameters ( Foam::argList args )
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


    IOdictionary PBDWsolverDict
    (
        IOobject
        (
            "PBDWsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read entries from PBDW dictionary
    const entry* existField = PBDWsolverDict.subDict("Offline_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = PBDWsolverDict.subDict("Offline_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existsensorsNumber= PBDWsolverDict.subDict("Offline_parameters").lookupEntryPtr("MaxSensorsNumber", false, true);
    const entry* existsensorsPositions= PBDWsolverDict.subDict("Offline_parameters").lookupEntryPtr("SensorsPositions", false, true);
    const entry* existsensorsVariance= PBDWsolverDict.subDict("Offline_parameters").lookupEntryPtr("SensorsVariance", false, true);
    const entry* existBasisNumber = PBDWsolverDict.subDict("Offline_parameters").lookupEntryPtr("BasisNumber", false, true);

    // Initiate variable from PBDWsolverDict

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' " << endl;
        abort();
    }

    word fieldName
    (
        PBDWsolverDict.subDict("Offline_parameters").lookup("field")
    );

    Info<<"FieldName: "<<fieldName<<endl;
    if (!(existBasisNumber))
    {
        Info << "Error: 'BasisNumber' entry not found. Must define 'BasisNumber' " << endl;
        abort();
    }

    scalar BasisNumber =
        readScalar
        (
            PBDWsolverDict.subDict("Offline_parameters").lookup("BasisNumber")
        );

    Info << "BasisNumber: " << BasisNumber << endl;

    if (!existFoldersLists)
    {
        Info << "Error: 'foldersList' entry not found. Must define the list of folders containing the snapshots" << endl;
        abort();
    }

    List<fileName> folders_list
    (
        PBDWsolverDict.subDict("Offline_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;


    Info<<"\nParameters specific for sensors selection"<<endl;

    if (!(existsensorsNumber) | !(existsensorsVariance) )
    {
        Info << "Error: you must define 'MaxSensorsNumber' and 'SensorsVariance'" << endl;
        abort();
    }

    scalar MaxSensorsNumber =
        readScalar
        (
            PBDWsolverDict.subDict("Offline_parameters").lookup("MaxSensorsNumber")
        );

    Info << "MaxSensorsNumber: " << MaxSensorsNumber << endl;

    scalar SensorsVariance =
        readScalar
        (
            PBDWsolverDict.subDict("Offline_parameters").lookup("SensorsVariance")
        );
    Info << "Sensors variance: " << SensorsVariance << endl;

    pointField sensorsPositions (mesh.C().field());

    if (!(existsensorsPositions))
    {
        Info << "No sensors location selected. As many sensors as Cells will be created  " << endl;
    }
    else
    {
        readSensorPositions(PBDWsolverDict,sensorsPositions);
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

    if (MaxSensorsNumber > sensorsPositions.size())
    {
        Info << "ERROR! MaxSensorsNumber cannot be greater than the number of sensors available " << nl
             << "new MaxSensorsNumber is set = " << sensorsPositions.size() << endl;
        MaxSensorsNumber=sensorsPositions.size();
    }

    Info<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<endl;
    PBDWparameters results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.MaxSensorsNumber= MaxSensorsNumber;
    results.sensorsPositions=sensorsPositions;
    results.SensorsVariance=SensorsVariance;
    results.BasisNumber=BasisNumber;

    return results;
}

