struct GEIM_VTOnlineParameters
{
    List<fileName> folders_list;
    label  msNumber;
};

GEIM_VTOnlineParameters getGEIM_VTOnlineParameters ( Foam::argList args )
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

    // Read entries from the GEIM_VT dictionary
    const entry* existFoldersLists = GEIM_VTsolverDict.subDict("Online_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existMagicFunctionsNumber = GEIM_VTsolverDict.subDict("Online_parameters").lookupEntryPtr("msNumber", false, true);

    // Initiate variable from GEIM_VTsolverDict
    if (!existMagicFunctionsNumber)
    {
        Info << "Error: must define number of magic functions and magic sensors to use" << endl;
        abort();
    }

    if (!(existFoldersLists))
    {
        Info << "Error: 'foldersList' entry not found. Must define list of folders inside which there is the field to reconstruct " << endl;
        abort();
    }


    List<fileName> folders_list
    (
        GEIM_VTsolverDict.subDict("Online_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;

    label MagicFunctionsNumber =
        readInt
        (
            GEIM_VTsolverDict.subDict("Online_parameters").lookup("msNumber")
        );
    Info << "Number of Magic Functions: " << MagicFunctionsNumber << endl;



    GEIM_VTOnlineParameters  results;
    results.folders_list=folders_list;
    results.msNumber=MagicFunctionsNumber;

    return results;
}

