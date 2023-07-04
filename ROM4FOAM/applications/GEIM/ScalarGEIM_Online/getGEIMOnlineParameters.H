struct GEIMOnlineParameters
{
    List<fileName> folders_list;
    word fieldName;
    label  msNumber;
};

GEIMOnlineParameters getGEIMOnlineParameters ( Foam::argList args )
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

    // Read entries from the GEIM dictionary
    const entry* existField = GEIMsolverDict.subDict("Online_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = GEIMsolverDict.subDict("Online_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existMagicFunctionsNumber = GEIMsolverDict.subDict("Online_parameters").lookupEntryPtr("msNumber", false, true);

    // Initiate variable from GEIMsolverDict
    if (!existMagicFunctionsNumber)
    {
        Info << "Error: must define number of magic functions and magic sensors to use" << endl;
        abort();
    }

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' to reconstruct " << endl;
        abort();
    }

    if (!(existFoldersLists))
    {
        Info << "Error: 'foldersList' entry not found. Must define list of folders inside which there is the field to reconstruct " << endl;
        abort();
    }


    List<fileName> folders_list
    (
        GEIMsolverDict.subDict("Online_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;

    word fieldName
    (
        GEIMsolverDict.subDict("Online_parameters").lookup("field")
    );
    Info << "Field to reconstruct: " << fieldName << endl;

    label MagicFunctionsNumber =
        readInt
        (
            GEIMsolverDict.subDict("Online_parameters").lookup("msNumber")
        );
    Info << "Number of Magic Functions: " << MagicFunctionsNumber << endl;



    GEIMOnlineParameters  results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.msNumber=MagicFunctionsNumber;

    return results;
}

