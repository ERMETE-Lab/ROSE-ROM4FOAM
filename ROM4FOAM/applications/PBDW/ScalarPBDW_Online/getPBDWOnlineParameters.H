struct PBDWOnlineParameters
{
    List<fileName> folders_list;
    fileName sensorsFolder;
    word fieldName;
    scalar BasisNumber;
    scalar MaxSensors;
};

PBDWOnlineParameters getPBDWOnlineParameters ( Foam::argList args )
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

    // Read entries from the PBDW dictionary
    const entry* existField = PBDWsolverDict.subDict("Online_parameters").lookupEntryPtr("field", false, true);
    const entry* existSensorsFolder = PBDWsolverDict.subDict("Online_parameters").lookupEntryPtr("sensorsFolder", false, true);
    const entry* existFoldersLists = PBDWsolverDict.subDict("Online_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existBasisNumber = PBDWsolverDict.subDict("Online_parameters").lookupEntryPtr("BasisNumber", false, true);
    const entry* existMaxSensors = PBDWsolverDict.subDict("Online_parameters").lookupEntryPtr("MaxSensors", false, true);
    
    // Initiate variable from PBDWsolverDict
    if (!existBasisNumber)
    {
        Info << "Error: must define number of basis functions to use" << endl;
        abort();
    }

    if (!existMaxSensors)
    {
        Info << "Error: must define number of magic sensors to use" << endl;
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
    if (!(existSensorsFolder))
    {
        Info << "Error: 'sensorsFolder' entry not found. Must define the folder of the sensors from which they can be loaded" << endl;
        abort();
    }


    List<fileName> folders_list
    (
        PBDWsolverDict.subDict("Online_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders " << folders_list << endl;

    word fieldName
    (
        PBDWsolverDict.subDict("Online_parameters").lookup("field")
    );
    Info << "Field to reconstruct: " << fieldName << endl;

     scalar MaxSensors =
        readScalar
        (
            PBDWsolverDict.subDict("Online_parameters").lookup("MaxSensors")
        );
    Info << "Max sensors number: " << MaxSensors << endl;

        scalar BasisNumber =
        readScalar
        (
            PBDWsolverDict.subDict("Online_parameters").lookup("BasisNumber")
        );
    Info << "Number of Basis Functions: " << BasisNumber << endl;

    fileName sensorsFolder
    (
        PBDWsolverDict.subDict("Online_parameters").lookup("sensorsFolder")
    );
    Info << "Reading sensors Folder " << sensorsFolder << endl;

    Info<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<endl;


    PBDWOnlineParameters  results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.sensorsFolder = sensorsFolder;
    results.BasisNumber=BasisNumber;
    results.MaxSensors = MaxSensors;

    return results;
}

