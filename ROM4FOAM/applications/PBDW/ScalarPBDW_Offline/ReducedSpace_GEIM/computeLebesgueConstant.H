
// compute Lagrangian Magic Functions

label mfsize= BasisFunctions.size();
LebesgueConstantList.setSize ( mfsize );
PtrList<volScalarField> OrthoNormalBasisFunctions (mfsize);
PtrList<volScalarField> OrthonormalMagicSensors (mfsize);

for(label ll=0; ll<mfsize; ll++)
{
    OrthoNormalBasisFunctions.set
    (
        ll,
        new volScalarField
        (
            BasisFunctions[ll]
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

        OrthoNormalBasisFunctions[ii]-= ((MOR::projection(OrthoNormalBasisFunctions[ii],OrthoNormalBasisFunctions[jj]))/(MOR::projection(OrthoNormalBasisFunctions[jj],OrthoNormalBasisFunctions[jj]))*(OrthoNormalBasisFunctions[jj]));

        OrthonormalMagicSensors[ii]-= ((MOR::projection(OrthonormalMagicSensors[ii],OrthonormalMagicSensors[jj]))/(MOR::projection(OrthonormalMagicSensors[jj],OrthonormalMagicSensors[jj]))*(OrthonormalMagicSensors[jj]));

    }
}

	//normalize

for(label ii=0; ii<mfsize; ++ii)
{
    OrthoNormalBasisFunctions[ii]/=MOR::L2norm(OrthoNormalBasisFunctions[ii]);
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
            A[ii][jj]=MOR::projection(OrthoNormalBasisFunctions[ii], OrthonormalMagicSensors[jj]);
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



