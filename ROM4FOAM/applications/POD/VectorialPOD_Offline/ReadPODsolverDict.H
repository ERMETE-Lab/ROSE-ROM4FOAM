struct PODparameters
{
    List<fileName> folders_list;
    word fieldName;
    label  maxBasis;
    scalar accuracy;
};

PODparameters getPODparameters (Foam::argList args)
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


    IOdictionary PODsolverDict
    (
        IOobject
        (
            "PODsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read Initial and last time from the POD dictionary
    const entry* existField = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existmaxBasis = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("maxBasis", false, true);
    const entry* existAccuracy = PODsolverDict.subDict("Offline_parameters").lookupEntryPtr("accuracy", false, true);

    // Initiate variable from PODsolverDict
    if (!(existmaxBasis) | !(existAccuracy))
    {
        Info << "Error: you must define accuracy and max number of Basis" << endl;
        abort();
    }

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' " << endl;
        abort();
    }
    
    if (!existFoldersLists)
    {
        Info << "Error: 'foldersList' entry not found. Must define the list of folders containing the snapshots" << endl;
        abort();
    }

    word fieldName
    (
        PODsolverDict.subDict("Offline_parameters").lookup("field")
    );

    label maxBasis =
        readInt
        (
            PODsolverDict.subDict("Offline_parameters").lookup("maxBasis")
        );
    Info << "Seeking number of Basis: " << maxBasis << endl;

    scalar accuracy =
        readScalar
        (
            PODsolverDict.subDict("Offline_parameters").lookup("accuracy")
        );

    Info << "Seeking accuracy: " << accuracy << endl;

    List<fileName> folders_list
    (
        PODsolverDict.subDict("Offline_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;

    PODparameters results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.maxBasis=maxBasis;
    results.accuracy=accuracy;

    return results;
}

