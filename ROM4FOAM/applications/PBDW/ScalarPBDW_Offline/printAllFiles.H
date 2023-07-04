// print files.txt

Info << "\nprinting files"<< endl;

if (isDir(args.rootPath()/currentFolder/fieldName+"_PBDW_Offline_files"))
{
    rmDir(args.rootPath()/currentFolder/fieldName+"_PBDW_Offline_files");
}

mkDir(args.rootPath()/currentFolder/fieldName+"_PBDW_Offline_files");
chDir(args.rootPath()/currentFolder/fieldName+"_PBDW_Offline_files");


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

ofstream checkMatrixFile (fieldName+"orthog_matrix.txt", std::ios::out | std::ios::app);
    if (checkMatrixFile.is_open())
    {
        for (label ii=0; ii< BasisFunctions.size(); ++ii)
        {
            for(label jj=0; jj< BasisFunctions.size(); ++jj)
            {
                checkMatrixFile << std::setprecision(14)<< checkMatrix(ii,jj)<<", ";
            }
            checkMatrixFile << "\n";
        }
    checkMatrixFile.close();
    }

if ( algorithmRS == "GEIM" )
{
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
    ofstream matrixB (fieldName+"matrixB.txt", std::ios::out | std::ios::app);
    if (matrixB.is_open())
    {
        for (label ii=0; ii< BasisFunctions.size(); ++ii)
        {
            for(label jj=0; jj< BasisFunctions.size(); ++jj)
            {
                matrixB << std::setprecision(14)<< Bmatrix(ii,jj)<<", ";
            }
            matrixB << "\n";
        }
    matrixB.close();
    }
}
else
{
    ofstream InfSupConstant ("InfSup_constant.txt", std::ios::out | std::ios::app);
    if (InfSupConstant.is_open())
    {
        forAll (infSupConstant, IndexI)
        {
            InfSupConstant << std::setprecision(14)<< infSupConstant[IndexI];
            InfSupConstant << "\n";
        }
        InfSupConstant.close();
    }
}

ofstream Coefficients (fieldName+"_Coefficients.txt", std::ios::out | std::ios::app);
if (Coefficients.is_open())
{
	forAll (scalarSnapshotsList, snapI)
	{
    		forAll (InterpolationCoeffList[snapI], ii)
    		{
        		Coefficients << std::setprecision(14)<< InterpolationCoeffList[snapI][ii]<<", ";
    		}
        	Coefficients << "\n";
	}
    	Coefficients.close();
}