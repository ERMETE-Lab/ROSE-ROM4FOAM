fileName currentFolder= args.caseName();


if(isDir(args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance)))
{
    Info<<"\n'GEIM_"<<fieldName<<"_s_"<<SensorsVariance<<"' folder already exists" << endl;
    Info<<"\n'deleting 'GEIM_"<<fieldName<<"' folder \n" << endl;
    rmDir(args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance));
}

Info<<"\ncreating 'GEIM_"<<fieldName<<"_s_"<<SensorsVariance<<"' folder \n" << endl;
mkDir(args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance));
mkDir(args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance)/"0");
cp(args.rootPath()/args.caseName()/"system", args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance)/"system");
cp(args.rootPath()/args.caseName()/"constant", args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance)/"constant");


chDir(args.rootPath()/"GEIM_"+fieldName+"_s_"+name(SensorsVariance));
currentFolder="GEIM_"+fieldName+"_s_"+name(SensorsVariance);


Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), currentFolder);
runTime.setTime(0,0);

