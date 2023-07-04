// print files.txt

Info << "\nprinting files in "<< args.caseName() << endl;

chDir(args.rootPath()/args.caseName());
if (isDir(fieldName+"_GEIM_Online_files"))
{
    rmDir((fieldName+"_GEIM_Online_files"));
}
mkDir(fieldName+"_GEIM_Online_files");
chDir(fieldName+"_GEIM_Online_files");

// L2 error
ofstream average_relative_error ("average_L2_relative_error.txt", std::ios::out | std::ios::app);
if (average_relative_error.is_open())
{

    forAll (averageRelL2ErrorList, IndexI)
    {
        average_relative_error << std::setprecision(14)<< averageRelL2ErrorList[IndexI];
        average_relative_error << "\n";
    }
    average_relative_error.close();
}

ofstream maximum_relative_error ("maximum_L2_relative_error.txt", std::ios::out | std::ios::app);
if (maximum_relative_error.is_open())
{

    forAll (maxRelL2ErrorList, IndexI)
    {
        maximum_relative_error << std::setprecision(14)<< maxRelL2ErrorList[IndexI];
        maximum_relative_error << "\n";
    }
    maximum_relative_error.close();
}

ofstream average_absolute_error ("average_L2_absolute_error.txt", std::ios::out | std::ios::app);
if (average_absolute_error.is_open())
{

    forAll (averageAbsL2ErrorList, IndexI)
    {
        average_absolute_error << std::setprecision(14)<< averageAbsL2ErrorList[IndexI];
        average_absolute_error << "\n";
    }
    average_absolute_error.close();
}

ofstream maximum_absolute_error ("maximum_L2_absolute_error.txt", std::ios::out | std::ios::app);
if (maximum_absolute_error.is_open())
{

    forAll (maxAbsL2ErrorList, IndexI)
    {
        maximum_absolute_error << std::setprecision(14)<< maxAbsL2ErrorList[IndexI];
        maximum_absolute_error << "\n";
    }
    maximum_absolute_error.close();
}
