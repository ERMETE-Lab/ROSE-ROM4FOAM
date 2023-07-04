// print files.txt

Info << "\nprinting files in "<< args.caseName() << endl;

chDir(args.rootPath()/args.caseName());
if (isDir(fieldName+"_reconstruction_POD_files"))
{
    rmDir((fieldName+"_reconstruction_POD_files"));
}
mkDir(fieldName+"_reconstruction_POD_files");
chDir(fieldName+"_reconstruction_POD_files");


ofstream average_L2_relative_error ("average_L2_relative_error.txt", std::ios::out | std::ios::app);
if (average_L2_relative_error.is_open())
{
    average_L2_relative_error << std::setprecision(14)<< "1";
    average_L2_relative_error << "\n";

    forAll (averageL2ErrorList, IndexI)
    {
        average_L2_relative_error << std::setprecision(14)<< averageL2ErrorList[IndexI];
        average_L2_relative_error << "\n";
    }
    average_L2_relative_error.close();
}

ofstream maximum_L2_relative_error ("maximum_L2_relative_error.txt", std::ios::out | std::ios::app);
if (maximum_L2_relative_error.is_open())
{
    maximum_L2_relative_error  << std::setprecision(14)<< "1";
    maximum_L2_relative_error  << "\n";
    forAll (maxL2ErrorList, IndexI)
    {
        maximum_L2_relative_error << std::setprecision(14)<< maxL2ErrorList[IndexI];
        maximum_L2_relative_error << "\n";
    }
    maximum_L2_relative_error.close();
}

ofstream average_L2_relative_error_Energy ("average_L2_relative_error_Energy.txt", std::ios::out | std::ios::app);
if (average_L2_relative_error_Energy.is_open())
{
    average_L2_relative_error_Energy << std::setprecision(14)<< "1";
    average_L2_relative_error_Energy << "\n";

    forAll (averageL2ErrorList_Energy, IndexI)
    {
        average_L2_relative_error_Energy << std::setprecision(14)<< averageL2ErrorList_Energy[IndexI];
        average_L2_relative_error_Energy << "\n";
    }
    average_L2_relative_error_Energy.close();
}

ofstream maximum_L2_relative_error_Energy ("maximum_L2_relative_error_Energy.txt", std::ios::out | std::ios::app);
if (maximum_L2_relative_error_Energy.is_open())
{
    forAll (maxL2ErrorList, IndexI)
    {
        maximum_L2_relative_error_Energy << std::setprecision(14)<< maxL2ErrorList_Energy[IndexI];
        maximum_L2_relative_error_Energy << "\n";
    }
    maximum_L2_relative_error_Energy.close();
}





