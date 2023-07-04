BasisFunctions.append
(
    GeometricField<scalar, fvPatchField, volMesh>
    (
        IOobject
        (
            scalarSnapshotsList[generatingFunction].name(),
            scalarSnapshotsList[generatingFunction].time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalarSnapshotsList[generatingFunction],
        fvPatchField<scalar>::calculatedType()
    )
);

for (label ii = 0; ii < iterIndex-1; ++ii)
{
    BasisFunctions[iterIndex - 1] -=   MOR::projection( scalarSnapshotsList[generatingFunction],BasisFunctions[ii]) / 
                                       MOR::projection( BasisFunctions[ii],                     BasisFunctions[ii]) *                                      BasisFunctions[ii];
}

BasisFunctions[iterIndex - 1] /= MOR::L2norm(BasisFunctions[iterIndex - 1]);
