fileName currentFolder= args.caseName();


if( isDir(args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance) ) )
{
    Info<<"\n'PBDW_"<<fieldName<<"_"<<algorithmRS<<"' folder already exists --> Deleting 'PBDW_"<<fieldName<<"_"<<algorithmRS<<"_s_"<<SensorsVariance<<"' folder" << endl;
    rmDir( args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance) );
}

Info<<"\nCreating 'PBDW_"<<fieldName<<"_"<<algorithmRS<<"_s_"<<SensorsVariance<<"' folder \n" << endl;
mkDir(args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance) );
mkDir(args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance)/"0");
cp(args.rootPath()/args.caseName()/"system", args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance)/"system");
cp(args.rootPath()/args.caseName()/"constant", args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance)/"constant");


chDir(args.rootPath()/"PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance));
currentFolder="PBDW_"+fieldName+"_"+algorithmRS+"_s_"+name(SensorsVariance);


Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), currentFolder);
runTime.setTime(0,0);

