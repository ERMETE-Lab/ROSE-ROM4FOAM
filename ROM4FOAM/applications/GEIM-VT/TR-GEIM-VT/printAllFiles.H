// print files.txt

Info << "\nprinting files in "<< args.caseName() << endl;

chDir(args.rootPath()/args.caseName());
if (! isDir("TR-GEIM-VT_files"))
{
    mkDir("TR-GEIM-VT_files");
}

chDir("TR-GEIM-VT_files");

if ( isFile("T_average_L2_relative_error_std_"+name(noise_std)+".txt"))
{
    rm ("T_average_L2_relative_error_std_"+name(noise_std)+".txt");
}

if ( isFile("U_average_L2_relative_error_std_"+name(noise_std)+".txt"))
{
    rm ("U_average_L2_relative_error_std_"+name(noise_std)+".txt");
}

if ( isFile("p_rgh_average_L2_relative_error_std_"+name(noise_std)+".txt"))
{
    rm ("p_rgh_average_L2_relative_error_std_"+name(noise_std)+".txt");
}


ofstream T_L2_average_relative_error ("T_average_L2_relative_error_std_"+name(noise_std)+".txt", std::ios::out | std::ios::app);
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

ofstream p_L2_average_relative_error ("p_rgh_average_L2_relative_error_std_"+name(noise_std)+".txt", std::ios::out | std::ios::app);
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


ofstream U_L2_average_relative_error ("U_average_L2_relative_error_std_"+name(noise_std)+".txt", std::ios::out | std::ios::app);
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


