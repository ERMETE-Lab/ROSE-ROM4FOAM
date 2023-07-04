

// compute Lagrangian Magic Functions

label mfsize = MagicFunctions.size();

List<scalar> LebesgueConstantList (mfsize);

// repeat Lebesgue constant  calculation for each additional basis; kk number of basis 

for(label kk=1; kk<=mfsize; ++kk)
{
    PtrList<volVectorField> LagrangianFunctions (kk);


    for(label ll=0; ll<kk; ll++)
    {
        LagrangianFunctions.set
        (
            ll,
            new volVectorField
            (
                MagicFunctions[ll]
            )
        );
    }

    // computing lagrangian basis

    for(label ii=1; ii<=kk; ++ii)
    {
        for (label jj= (kk-ii); jj< kk-1; jj++)
        {

            LagrangianFunctions[kk-ii]-= Bmatrix(jj+1,kk-ii)*(LagrangianFunctions[jj+1]);
        }
    }

    // computing Lebesgue Constant

    volScalarField tmp_Field
    (
        IOobject
        (
            "tmp_field",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<scalar>
        (
            "zero",
            LagrangianFunctions[0].dimensions(),
            pTraits<scalar>::zero
        )
    );

    for (label jj= 0; jj<kk; ++jj)
    {
        tmp_Field+= mag(LagrangianFunctions[jj]);
    }

    LebesgueConstantList[kk-1]=max(tmp_Field.internalField()).value();

}



