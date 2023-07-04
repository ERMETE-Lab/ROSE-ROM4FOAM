// print files.txt

Info << "\nprinting files"<< endl;

mkDir(args.rootPath()/currentFolder/"GEIM-VT_Offline_files");
chDir(args.rootPath()/currentFolder/"GEIM-VT_Offline_files");


ofstream relative_error ("max_relative_error.txt", std::ios::out | std::ios::app);
if (relative_error.is_open())
{
    forAll (relativeErrorList, IndexI)
    {
        relative_error << std::setprecision(14)<< relativeErrorList[IndexI];
        relative_error << "\n";
    }
    relative_error.close();
}

ofstream Chosen_fields ("Chosen_fields.txt", std::ios::out | std::ios::app);
if (Chosen_fields.is_open())
{
    forAll (whichFieldisCorrected, IndexI)
    {
        Chosen_fields << whichFieldisCorrected[IndexI];
        Chosen_fields << "\n";
    }
    Chosen_fields.close();
}

ofstream GEIM_Coeffs_avarage_values ("GEIM-VT_T_Offline_Coeffs_avarage_values.txt", std::ios::out | std::ios::app);
if (GEIM_Coeffs_avarage_values.is_open())
{
    GEIM_Coeffs_avarage_values << std::setprecision(14)<< T_MagicFunctions.size() ;
    GEIM_Coeffs_avarage_values<< "\n";

    GEIM_Coeffs_avarage_values << "(" ;
    GEIM_Coeffs_avarage_values<< "\n";

    for (label ii=0; ii<T_MagicFunctions.size(); ii++)
    {

        List<scalar> temp_vector (T_snapshots.size());

        for (label jj=0; jj< T_snapshots.size(); jj++)
        {
            temp_vector[jj]=InterpolationCoeffList[jj][ii] ;
        }

        GEIM_Coeffs_avarage_values << std::setprecision(14)<< average(temp_vector) ;
        GEIM_Coeffs_avarage_values<< "\n";
    }

    GEIM_Coeffs_avarage_values << ")" ;

    GEIM_Coeffs_avarage_values.close();
}

ofstream GEIM_Coeffs_std ("GEIM-VT_T_Offline_Coeffs_std.txt", std::ios::out | std::ios::app);
if (GEIM_Coeffs_std.is_open())
{
    GEIM_Coeffs_std << std::setprecision(14)<< T_MagicFunctions.size() ;
    GEIM_Coeffs_std<< "\n";

    GEIM_Coeffs_std << "(" ;
    GEIM_Coeffs_std<< "\n";

    for (label ii=0; ii<T_MagicFunctions.size(); ii++)
    {

        List<scalar> temp_vector (T_snapshots.size());

        for (label jj=0; jj< T_snapshots.size(); jj++)
        {
            temp_vector[jj]=InterpolationCoeffList[jj][ii] ;
        }

        scalar mean = average(temp_vector);
        scalar variance =0;

        forAll (temp_vector, index)
        {
            variance+= std::pow(temp_vector[index]-mean,2);
        }

        variance/=(T_snapshots.size()-1);

        scalar standrd_deviation= std::sqrt(variance);

        GEIM_Coeffs_std << std::setprecision(14)<< standrd_deviation;
        GEIM_Coeffs_std<< "\n";
    }

    GEIM_Coeffs_std << ")" ;

    GEIM_Coeffs_std.close();
}

ofstream LebesgueConstant ("Lebesgue_constant.txt", std::ios::out | std::ios::app);
if (LebesgueConstant.is_open())
{
    forAll (LebesgueConstantList, IndexI)
    {
        LebesgueConstant << std::setprecision(14)<< LebesgueConstantList[IndexI];
        LebesgueConstant << "\n";
    }
    LebesgueConstant.close();
}

