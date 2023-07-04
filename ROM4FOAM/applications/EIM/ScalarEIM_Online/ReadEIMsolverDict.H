struct EIMOnlineParameters
{
    List<fileName> folders_list;
    word fieldName;
    label  mfNumber;
};


EIMOnlineParameters getEIMOnlineParameters ( Foam::argList args )
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
    const entry* existField = EIMsolverDict.subDict("Online_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = EIMsolverDict.subDict("Online_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existMagicFunctionsNumber = EIMsolverDict.subDict("Online_parameters").lookupEntryPtr("mfNumber", false, true);

    // Initiate variable from EIMsolverDict

    if (!existMagicFunctionsNumber)
    {
        Info << "Error: must define 'mfNumber' number of magic functions to use" << endl;
        abort();
    }

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' to reconstruct " << endl;
        abort();
    }

    if (!(existFoldersLists))
    {
        Info << "Error: 'foldersList' entry not found. Must define the list of folders containing the snapshots to reconstruct" << endl;
        abort();
    }

    List<fileName> folders_list
    (
        EIMsolverDict.subDict("Online_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;

    word fieldName
    (
        EIMsolverDict.subDict("Online_parameters").lookup("field")
    );
    Info << "Field to reconstruct: " << fieldName << endl;

    label MagicFunctionsNumber =
        readInt
        (
            EIMsolverDict.subDict("Online_parameters").lookup("mfNumber")
        );
    Info << "Number of Magic Functions: " << MagicFunctionsNumber << endl;



    EIMOnlineParameters  results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.mfNumber=MagicFunctionsNumber;

    return results;
}

