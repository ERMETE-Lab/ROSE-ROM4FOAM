// print files.txt

Info << "\nprinting files"<< endl;

if (isDir(args.rootPath()/currentFolder/fieldName+"_GEIM_Offline_files"))
{
    rmDir(args.rootPath()/currentFolder/fieldName+"_GEIM_Offline_files");
}

mkDir(args.rootPath()/currentFolder/fieldName+"_GEIM_Offline_files");
chDir(args.rootPath()/currentFolder/fieldName+"_GEIM_Offline_files");

ofstream GEIM_Coeffs_avarage_values ("GEIM_Offline_Coeffs_average_values.txt", std::ios::out | std::ios::app);
if (GEIM_Coeffs_avarage_values.is_open())
{
    GEIM_Coeffs_avarage_values << std::setprecision(14)<< MagicFunctions.size() ;
    GEIM_Coeffs_avarage_values<< "\n";

    GEIM_Coeffs_avarage_values << "(" ;
    GEIM_Coeffs_avarage_values<< "\n";

    for (label ii=0; ii<MagicFunctions.size(); ii++)
    {

        List<scalar> temp_vector (scalarSnapshotsList.size());

        for (label jj=0; jj< scalarSnapshotsList.size(); jj++)
        {
            temp_vector[jj]=InterpolationCoeffList[jj][ii] ;
        }

        GEIM_Coeffs_avarage_values << std::setprecision(14)<< average(temp_vector) ;
        GEIM_Coeffs_avarage_values<< "\n";
    }

    GEIM_Coeffs_avarage_values << ")" ;

    GEIM_Coeffs_avarage_values.close();
}

ofstream GEIM_Coeffs_std ("GEIM_Offline_Coeffs_std.txt", std::ios::out | std::ios::app);
if (GEIM_Coeffs_std.is_open())
{
    GEIM_Coeffs_std << std::setprecision(14)<< MagicFunctions.size() ;
    GEIM_Coeffs_std<< "\n";

    GEIM_Coeffs_std << "(" ;
    GEIM_Coeffs_std<< "\n";

    for (label ii=0; ii<MagicFunctions.size(); ii++)
    {

        List<scalar> temp_vector (scalarSnapshotsList.size());

        for (label jj=0; jj< scalarSnapshotsList.size(); jj++)
        {
            temp_vector[jj]=InterpolationCoeffList[jj][ii] ;
        }

        scalar mean = average(temp_vector);
        scalar variance =0;

        forAll (temp_vector, index)
        {
            variance+= std::pow(temp_vector[index]-mean,2);
        }

        variance/=(scalarSnapshotsList.size()-1);

        scalar standrd_deviation= std::sqrt(variance);

        GEIM_Coeffs_std << std::setprecision(14)<< standrd_deviation;
        GEIM_Coeffs_std<< "\n";
    }

    GEIM_Coeffs_std << ")" ;

    GEIM_Coeffs_std.close();
}


ofstream absolute_error ("max_absolute_L2_error.txt", std::ios::out | std::ios::app);
if (absolute_error.is_open())
{
    forAll (absoluteErrorList, IndexI)
    {
        absolute_error << std::setprecision(14)<< absoluteErrorList[IndexI];
        absolute_error << "\n";
    }
    absolute_error.close();
}

ofstream relative_error ("max_relative_L2_error.txt", std::ios::out | std::ios::app);
if (relative_error.is_open())
{
    forAll (relativeErrorList, IndexI)
    {
        relative_error << std::setprecision(14)<< relativeErrorList[IndexI];
        relative_error << "\n";
    }
    relative_error.close();
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


ofstream betaCoefficients (fieldName+"betaCoefficients.txt", std::ios::out | std::ios::app);
if (betaCoefficients.is_open())
{
	forAll (scalarSnapshotsList, snapI)
	{
    		forAll (InterpolationCoeffList[snapI], ii)
    		{
        		betaCoefficients << std::setprecision(14)<< InterpolationCoeffList[snapI][ii]<<", ";
    		}
        	betaCoefficients << "\n";
	}
    	betaCoefficients.close();
}

ofstream matrixB (fieldName+"matrixB.txt", std::ios::out | std::ios::app);
if (matrixB.is_open())
{
    for (label ii=0; ii< MagicFunctions.size(); ++ii)
    {
        for(label jj=0; jj< MagicFunctions.size(); ++jj)
        {
            matrixB << std::setprecision(14)<< Bmatrix(ii,jj)<<", ";
        }
        matrixB << "\n";
    }
matrixB.close();
}