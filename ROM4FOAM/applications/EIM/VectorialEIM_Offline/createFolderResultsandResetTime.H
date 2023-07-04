fileName currentFolder= args.caseName();

if(isDir(args.rootPath()/"EIM_"+fieldName))
{
    Info<<"\n'EIM_"<<fieldName<<"' folder already exists" << endl;
    Info<<"\n'deleting 'EIM_"<<fieldName<<"' folder \n" << endl;
    rmDir(args.rootPath()/"EIM_"+fieldName);
}

Info<<"\ncreating 'EIM_"<<fieldName<<"' folder \n" << endl;
mkDir(args.rootPath()/"EIM_"+fieldName);
mkDir(args.rootPath()/"EIM_"+fieldName/"0");
cp(args.rootPath()/args.caseName()/"system", args.rootPath()/"EIM_"+fieldName/"system");
cp(args.rootPath()/args.caseName()/"constant", args.rootPath()/"EIM_"+fieldName/"constant");


chDir(args.rootPath()/"EIM_"+fieldName);
currentFolder="EIM_"+fieldName;


Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), currentFolder);
runTime.setTime(0,0);

