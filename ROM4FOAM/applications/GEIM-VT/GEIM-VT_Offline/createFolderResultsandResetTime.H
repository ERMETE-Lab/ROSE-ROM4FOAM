fileName currentFolder= args.caseName();


if(isDir(args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance)))
{
    Info<<"\n'GEIM-VT_s_"<<SensorsVariance<<"' folder already exists" << endl;
    Info<< "deleting 'GEIM-VT_s_"<<SensorsVariance<<"' folder"<< endl;
    rmDir (args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance));
}

Info<<"\ncreating 'GEIM-VT_s_"<<SensorsVariance<<"' folder \n" << endl;
mkDir(args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance));
mkDir(args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance)/"0");
cp(args.rootPath()/args.caseName()/"system", args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance)/"system");
cp(args.rootPath()/args.caseName()/"constant", args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance)/"constant");


chDir(args.rootPath()/"GEIM-VT_s_"+name(SensorsVariance));
currentFolder="GEIM-VT_s_"+name(SensorsVariance);



Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), currentFolder);
runTime.setTime(0,0);

