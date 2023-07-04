struct EIMparameters
{
    List<fileName> folders_list;
    word fieldName;
    label  maxBasis;
    scalar error;
};

EIMparameters getEIMparameters ( Foam::argList args )
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

    IOdictionary EIMsolverDict
    (
        IOobject
        (
            "EIMsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read EIM dictionary entries
    const entry* existField = EIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = EIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existmaxBasis = EIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("maxBasis", false, true);
    const entry* existerror = EIMsolverDict.subDict("Offline_parameters").lookupEntryPtr("error", false, true);

    // Initiate variable from EIMsolverDict
    if (!(existmaxBasis) | !(existerror))
    {
        Info << "Error: must define 'error' and 'maxBasis'" << endl;
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
        EIMsolverDict.subDict("Offline_parameters").lookup("field")
    );

    label maxBasis =
        readInt
        (
            EIMsolverDict.subDict("Offline_parameters").lookup("maxBasis")
        );
    Info << "Seeking number of Basis: " << maxBasis << endl;

    scalar error =
        readScalar
        (
            EIMsolverDict.subDict("Offline_parameters").lookup("error")
        );

    Info << "Seeking error: " << error << endl;

    List<fileName> folders_list
    (
        EIMsolverDict.subDict("Offline_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;


    EIMparameters results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.maxBasis=maxBasis;
    results.error=error;

    return results;
}

