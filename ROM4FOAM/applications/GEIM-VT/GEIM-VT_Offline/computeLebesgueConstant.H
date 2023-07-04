

// compute Lagrangian Magic Functions

label mfsize= T_MagicFunctions.size();

List<scalar> LebesgueConstantList (mfsize);

PtrList<volScalarField> OrthoNormalT_MagicFunctions (mfsize);
PtrList<volScalarField> OrthonormalMagicSensors (mfsize);

for(label ll=0; ll<mfsize; ll++)
{
    OrthoNormalT_MagicFunctions.set
    (
        ll,
        new volScalarField
        (
            T_MagicFunctions[ll]
        )
    );

    OrthonormalMagicSensors.set
    (
        ll,
        new volScalarField
        (
            MagicSensors[ll]
        )
    );
}

// Gram-schmidt orthonormalization

for(label ii=0; ii<mfsize; ii++)
{
    for (label jj= 0; jj< ii; jj++)
    {

        OrthoNormalT_MagicFunctions[ii]-= ((MOR::projection(OrthoNormalT_MagicFunctions[ii],OrthoNormalT_MagicFunctions[jj]))/(MOR::projection(OrthoNormalT_MagicFunctions[jj],OrthoNormalT_MagicFunctions[jj]))*(OrthoNormalT_MagicFunctions[jj]));

        OrthonormalMagicSensors[ii]-= ((MOR::projection(OrthonormalMagicSensors[ii],OrthonormalMagicSensors[jj]))/(MOR::projection(OrthonormalMagicSensors[jj],OrthonormalMagicSensors[jj]))*(OrthonormalMagicSensors[jj]));

    }
}

	//normalize

for(label ii=0; ii<mfsize; ++ii)
{
    OrthoNormalT_MagicFunctions[ii]/=MOR::L2norm(OrthoNormalT_MagicFunctions[ii]);
    OrthonormalMagicSensors[ii]/=MOR::L2norm(OrthonormalMagicSensors[ii]);
}


// compute Lebesgue Constant

for(label kk=1; kk<=mfsize; ++kk)
{
    scalarSquareMatrix A (kk);

    for(label ii=0; ii<kk; ++ii)
    {
        for (label jj= 0; jj< kk; jj++)
        {
            A[ii][jj]=MOR::projection(OrthoNormalT_MagicFunctions[ii], OrthonormalMagicSensors[jj]);
        }
    }

    A = ((A.T())*A) ;

    // Calculate eigen-values

    EigenSolver<scalar> eigenSolver(A);

    // Sort and assemble

    SortableList<scalar> sortedList (A.n());

    forAll (sortedList, i)
    {
        sortedList[i] = eigenSolver.eigenValue(i);
    }

    // Sort will sort the values in ascending order and insert values
    sortedList.sort();

    //
    LebesgueConstantList[kk-1]= (1/std::sqrt(sortedList[0]));

}



