// print files.txt

Info << "\nprinting files in "<< args.caseName() << endl;

chDir(args.rootPath()/args.caseName());
if (isDir("GEIM-VT_Online_files"))
{
    rmDir(("GEIM-VT_Online_files"));
}
mkDir("GEIM-VT_Online_files");
chDir("GEIM-VT_Online_files");

ofstream T_L2_average_relative_error ("T_average_L2_relative_error.txt", std::ios::out | std::ios::app);
if (T_L2_average_relative_error.is_open())
{
    T_L2_average_relative_error << std::setprecision(14)<< "1";
    T_L2_average_relative_error << "\n";

    forAll (T_averageL2ErrorList, IndexI)
    {
        T_L2_average_relative_error << std::setprecision(14)<< T_averageL2ErrorList[IndexI];
        T_L2_average_relative_error << "\n";
    }
    T_L2_average_relative_error.close();
}

ofstream T_L2_maximum_relative_error ("T_maximum_L2_relative_error.txt", std::ios::out | std::ios::app);
if (T_L2_maximum_relative_error.is_open())
{
    T_L2_maximum_relative_error  << std::setprecision(14)<< "1";
    T_L2_maximum_relative_error  << "\n";
    forAll (T_maxL2ErrorList, IndexI)
    {
        T_L2_maximum_relative_error << std::setprecision(14)<< T_maxL2ErrorList[IndexI];
        T_L2_maximum_relative_error << "\n";
    }
    T_L2_maximum_relative_error.close();
}

ofstream p_L2_average_relative_error ("p_rgh_average_L2_relative_error.txt", std::ios::out | std::ios::app);
if (p_L2_average_relative_error.is_open())
{
    p_L2_average_relative_error << std::setprecision(14)<< "1";
    p_L2_average_relative_error << "\n";

    forAll (p_averageL2ErrorList, IndexI)
    {
        p_L2_average_relative_error << std::setprecision(14)<< p_averageL2ErrorList[IndexI];
        p_L2_average_relative_error << "\n";
    }
    p_L2_average_relative_error.close();
}

ofstream p_L2_maximum_relative_error ("p_rgh_maximum_L2_relative_error.txt", std::ios::out | std::ios::app);
if (p_L2_maximum_relative_error.is_open())
{
    p_L2_maximum_relative_error  << std::setprecision(14)<< "1";
    p_L2_maximum_relative_error  << "\n";
    forAll (p_maxL2ErrorList, IndexI)
    {
        p_L2_maximum_relative_error << std::setprecision(14)<< p_maxL2ErrorList[IndexI];
        p_L2_maximum_relative_error << "\n";
    }
    p_L2_maximum_relative_error.close();
}

ofstream U_L2_average_relative_error ("U_average_L2_relative_error.txt", std::ios::out | std::ios::app);
if (U_L2_average_relative_error.is_open())
{
    U_L2_average_relative_error << std::setprecision(14)<< "1";
    U_L2_average_relative_error << "\n";

    forAll (U_averageL2ErrorList, IndexI)
    {
        U_L2_average_relative_error << std::setprecision(14)<< U_averageL2ErrorList[IndexI];
        U_L2_average_relative_error << "\n";
    }
    U_L2_average_relative_error.close();
}

ofstream U_L2_maximum_relative_error ("U_maximum_L2_relative_error.txt", std::ios::out | std::ios::app);
if (U_L2_maximum_relative_error.is_open())
{
    U_L2_maximum_relative_error  << std::setprecision(14)<< "1";
    U_L2_maximum_relative_error  << "\n";
    forAll (U_maxL2ErrorList, IndexI)
    {
        U_L2_maximum_relative_error << std::setprecision(14)<< U_maxL2ErrorList[IndexI];
        U_L2_maximum_relative_error << "\n";
    }
    U_L2_maximum_relative_error.close();
}


