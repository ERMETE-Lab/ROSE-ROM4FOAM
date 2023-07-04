// print files.txt

Info << "\nprinting files in "<< args.caseName() << endl;

chDir(args.rootPath()/args.caseName());
if (! isDir(fieldName+"_TR-GEIM_files"))
{
    mkDir(fieldName+"_TR-GEIM_files");
}

chDir(fieldName+"_TR-GEIM_files");

if ( isFile("average_L2_relative_error_noiseStd_"+name(noise_std)+".txt"))
{
    rm("average_L2_relative_error_noiseStd_"+name(noise_std)+".txt");
}


ofstream average_relative_error ("average_L2_relative_error_noiseStd_"+name(noise_std)+".txt", std::ios::out | std::ios::app);
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


