fileName currentFolder= args.caseName ();

    if(isDir(args.rootPath()/"POD_"+fieldName))
    {
        Info<<"\n'POD_"<<fieldName<<"' folder already exists" << endl;
        Info<<"\n'deleting 'POD_"<<fieldName<<"' folder" << endl;
        rmDir(args.rootPath()/"POD_"+fieldName);
    }

    Info<<"\nCreating 'POD_"<<fieldName<<"' folder \n" << endl;
    mkDir(args.rootPath()/"POD_"+fieldName);
    mkDir(args.rootPath()/"POD_"+fieldName/"0");
    cp(args.rootPath()/args.caseName()/"system", args.rootPath()/"POD_"+fieldName/"system");
    cp(args.rootPath()/args.caseName()/"constant", args.rootPath()/"POD_"+fieldName/"constant");


    chDir(args.rootPath()/"POD_"+fieldName);
    currentFolder="POD_"+fieldName;


Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), currentFolder);
runTime.setTime(0,0);
