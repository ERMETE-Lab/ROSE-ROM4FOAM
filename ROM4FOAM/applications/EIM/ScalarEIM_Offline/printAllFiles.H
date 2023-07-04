// print files.txt

Info << "\nprinting files"<< endl;

if (! isDir(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files"))
{
    mkDir(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files");
}

chDir(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files");


if ( isFile(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/absolute_L_infinity_error.txt"))
{
    rm(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/absolute_L_infinity_error.txt");
}

if ( isFile(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/relative_L_infinity_error.txt"))
{
    rm(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/relative_L_infinity_error.txt");
}

if ( isFile(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/Lebesgue_constant.txt"))
{
    rm(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/Lebesgue_constant.txt");
}

if ( isFile(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/MagicParameter.txt"))
{
    rm(args.rootPath()/currentFolder/fieldName+"_EIM_Offline_files/MagicParameter.txt");
}

ofstream absolute_error ("max_absolute_L_infinity_error.txt", std::ios::out | std::ios::app);
if (absolute_error.is_open())
{
    forAll (absoluteErrorList, IndexI)
    {
        absolute_error << absoluteErrorList[IndexI];
        absolute_error << "\n";
    }
    absolute_error.close();
}

ofstream relative_error ("max_relative_L_infinity_error.txt", std::ios::out | std::ios::app);
if (relative_error.is_open())
{
    forAll (relativeErrorList, IndexI)
    {
        relative_error << relativeErrorList[IndexI];
        relative_error << "\n";
    }
    relative_error.close();
}

ofstream LebesgueConstant ("Lebesgue_constant.txt", std::ios::out | std::ios::app);
if (LebesgueConstant.is_open())
{
    forAll (LebesgueConstantList, IndexI)
    {
        LebesgueConstant << LebesgueConstantList[IndexI];
        LebesgueConstant << "\n";
    }
    LebesgueConstant.close();
}

// the index of the generating function is saved into a .txt file
ofstream Magic_Parameter ("MagicParameter.txt", std::ios::out | std::ios::app);
if (Magic_Parameter.is_open())
{
    forAll (idxGeneratingFunction, IndexI)
    {
        Magic_Parameter << idxGeneratingFunction[IndexI];
        Magic_Parameter << "\n";
    }
    Magic_Parameter.close();
}
