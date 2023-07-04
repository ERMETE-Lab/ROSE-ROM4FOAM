// print files.txt

Info << "\nprinting files in "<< args.caseName() << endl;

chDir(args.rootPath()/args.caseName());


if (! isDir(fieldName+"_EIM_Online_files"))
{
    mkDir(fieldName+"_EIM_Online_files");
}

if ( isFile(fieldName+"_EIM_Online_files/average_L2_relative_error.txt"))
{
    rm(fieldName+"_EIM_Online_files/average_L2_relative_error.txt");
}

if ( isFile(fieldName+"_EIM_Online_files/maximum_L2_relative_error.txt"))
{
    rm(fieldName+"_EIM_Online_files/maximum_L2_relative_error.txt");
}

chDir(fieldName+"_EIM_Online_files");

ofstream average_relative_error ("average_L2_relative_error.txt", std::ios::out | std::ios::app);
if (average_relative_error.is_open())
{
    average_relative_error << std::setprecision(14)<< "1";
    average_relative_error << "\n";

    forAll (averageL2ErrorList, IndexI)
    {
        average_relative_error << std::setprecision(14)<< averageL2ErrorList[IndexI];
        average_relative_error << "\n";
    }
    average_relative_error.close();
}

ofstream maximum_relative_error ("maximum_L2_relative_error.txt", std::ios::out | std::ios::app);
if (maximum_relative_error.is_open())
{
    maximum_relative_error  << std::setprecision(14)<< "1";
    maximum_relative_error  << "\n";
    forAll (maxL2ErrorList, IndexI)
    {
        maximum_relative_error << std::setprecision(14)<< maxL2ErrorList[IndexI];
        maximum_relative_error << "\n";
    }
    maximum_relative_error.close();
}

