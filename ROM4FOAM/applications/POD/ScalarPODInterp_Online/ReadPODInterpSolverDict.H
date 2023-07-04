struct PODOnlineParameters
{
    List<fileName> folders_list;
    word fieldName;
    label  BasisNumber;

};

PODOnlineParameters getPODOnlineParameters (Foam::argList args)
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


    IOdictionary PODInterpSolverDict
    (
        IOobject
        (
            "PODInterpSolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read entries from the POD dictionary
    const entry* existField = PODInterpSolverDict.subDict("Online_parameters").lookupEntryPtr("field", false, true);
    const entry* existFoldersLists = PODInterpSolverDict.subDict("Online_parameters").lookupEntryPtr("foldersList", false, true);
    const entry* existNumberBasis = PODInterpSolverDict.subDict("Online_parameters").lookupEntryPtr("BasisNumber", false, true);

    // Initiate variable from PODsolverDict
    if (!(existNumberBasis) )
    {
        Info << "Error: you must define number of Basis: 'BasisNumber'" << endl;
        abort();
    }

    if (!(existField))
    {
        Info << "Error: 'field' entry not found. Must define 'field' " << endl;
        abort();
    }
    
    if (!existFoldersLists)
    {
        Info << "Error: 'foldersList' entry not found. Must define the list of folders containing the snapshots to reconstruct" << endl;
        abort();
    }


    word fieldName
    (
        PODInterpSolverDict.subDict("Online_parameters").lookup("field")
    );

    label BasisNumber =
        readInt
        (
            PODInterpSolverDict.subDict("Online_parameters").lookup("BasisNumber")
        );
    Info << "Number of Basis: " << BasisNumber << endl;
    
    List<fileName> folders_list
    (
        PODInterpSolverDict.subDict("Online_parameters").lookup("foldersList")
    );
    Info << "Reading list of folders\n" << folders_list << endl;


    PODOnlineParameters results;
    results.folders_list=folders_list;
    results.fieldName=fieldName;
    results.BasisNumber=BasisNumber;

    return results;
}

