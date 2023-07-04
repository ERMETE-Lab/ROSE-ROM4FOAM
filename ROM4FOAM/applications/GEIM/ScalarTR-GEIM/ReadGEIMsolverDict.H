struct GEIMOnlineParameters
{
    List<fileName> folders_list;
    word fieldName;
    label  msNumber;
    scalar noise_std;
    label NumberRepeatedExperiments;

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
    const entry* existNoise = GEIMsolverDict.subDict("Online_parameters").lookupEntryPtr("noise_std", false, true);
    const entry* existNRepeatedExp= GEIMsolverDict.subDict("Online_parameters").lookupEntryPtr("N_Repeated_Experiments", false, true);



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

    if (!(existNoise))
    {
        Info << "Error: 'noise_std' entry not found. Must define noise standard deviation " << endl;
        abort();
    }

    if (!(existNRepeatedExp))
    {
        Info << "Error: 'N_Repeated_Experiments' entry not found. Must define number of Repeated experiments " << endl;
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

    scalar noise_std=
        readScalar
        (
            GEIMsolverDict.subDict("Online_parameters").lookup("noise_std")
        );
    Info << "Gaussian noise standart deviation: " << noise_std << endl;

    label NumberRepeatedExperiments =
        readInt
        (
            GEIMsolverDict.subDict("Online_parameters").lookup("N_Repeated_Experiments")
        );
    Info << "Number of Repeated experiments: " << NumberRepeatedExperiments << endl;



    GEIMOnlineParameters  results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.msNumber=MagicFunctionsNumber;
    results.noise_std=noise_std;
    results.NumberRepeatedExperiments=NumberRepeatedExperiments;


    return results;
}

