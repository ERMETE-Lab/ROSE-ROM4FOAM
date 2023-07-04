

// compute Lagrangian Magic Functions

label mfsize= MagicFunctions.size();

List<scalar> LebesgueConstantList (mfsize);

// repeat Lebesgue constant calculation for each additional basis; kk number of basis 

for(label kk=1; kk<=mfsize; ++kk)
{
    PtrList<volScalarField> LagrangianFunctions (kk);

    //copia le kk magicFunctions in nuova PtrList

    for(label ll=0; ll<kk; ll++)
    {
        LagrangianFunctions.set
        (
            ll,
            new volScalarField
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
            LagrangianFunctions[0].name(),
            LagrangianFunctions[0].time().timeName(),
            LagrangianFunctions[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        LagrangianFunctions[0].mesh(),
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



