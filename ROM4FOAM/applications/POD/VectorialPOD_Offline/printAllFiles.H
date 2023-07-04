// print files.txt

Info << "\nprinting files"<< endl;

if (isDir(args.rootPath()/currentFolder/fieldName+"_POD_Offline_files"))
{
    rmDir(args.rootPath()/currentFolder/fieldName+"_POD_Offline_files");
}

mkDir(args.rootPath()/currentFolder/fieldName+"_POD_Offline_files");
chDir(args.rootPath()/currentFolder/fieldName+"_POD_Offline_files");


ofstream EigenValuesfile (fieldName+"EigenValues.txt", std::ios::out | std::ios::app );
if (EigenValuesfile.is_open())
{
    for (label jj=0; jj< EigenValues.size(); jj++)
    {
        EigenValuesfile << std::setprecision(14)<< EigenValues[jj] ;
        EigenValuesfile << "\n";

    }
    EigenValuesfile.close();
}

ofstream coefficientsFiles (fieldName+"alphaCoeffs.txt", std::ios::out | std::ios::app );
if (coefficientsFiles.is_open())
{
    forAll (vectorSnapshotsList, snapI)
    {
	    for (label baseI = 0; baseI < vectorPODOrthoNormalBasis.size(); baseI++)
	    {
		coefficientsFiles << std::setprecision(14)<< coeffs[snapI][baseI] << ", " ;
	    }
    coefficientsFiles << "\n";
    }
    coefficientsFiles.close();
}

ofstream absErrorFile (fieldName+"L2AbsError.txt", std::ios::out | std::ios::app );
if (absErrorFile.is_open())
{
    forAll (L2AbsErrorList, baseI)
    {
	    absErrorFile << std::setprecision(14)<< L2AbsErrorList[baseI];
    	absErrorFile << "\n";
    }
    absErrorFile.close();
}

ofstream relErrorFile (fieldName+"L2RelError.txt", std::ios::out | std::ios::app );
if (relErrorFile.is_open())
{
    forAll (L2RelErrorList, baseI)
    {
	    relErrorFile << std::setprecision(14)<< L2RelErrorList[baseI];
    	relErrorFile << "\n";
    }
    relErrorFile.close();
}
