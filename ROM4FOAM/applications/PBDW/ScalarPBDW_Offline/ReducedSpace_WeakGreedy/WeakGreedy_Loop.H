
 // ++++++++++++++++++++++++++++++++++++++ WeakGreedy First Iteration +++++++++++++++++++++++++++++++ 
    
    //update index
    ++iterIndex;

    forAll (scalarSnapshotsList, snapI)
    {

        if(max_L2norm< MOR::L2norm(scalarSnapshotsList[snapI]))
        {
            max_L2norm = MOR::L2norm(scalarSnapshotsList[snapI]);
            generatingFunction=snapI;
        }

    }

    Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
    Info <<"maximum absolute error in L2 norm:     " << setprecision(14) << max_L2norm <<endl;
    Info <<"maximum relative error in L2 norm:     " << max_L2norm/max_L2norm  <<endl;

    //update errorLists
    absoluteErrorList.append(max_L2norm);
    relativeErrorList.append(max_L2norm/max_L2norm);


    // Update Magic Quantites
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
            scalarSnapshotsList[generatingFunction]/max_L2norm,
            fvPatchField<scalar>::calculatedType()
        )
    );

 // ++++++++++++++++++++++++++++++++++++++ WeakGreedy Main Loop +++++++++++++++++++++++++++++++ 
    
    scalar max_relative_error;
    scalar max_absolute_error ;
    InterpolationCoeffList.setSize(scalarSnapshotsList.size());
    scalar snapshotsMeasure = SMALL;
    label coeffIndex=0;

    while ( iterIndex <= BasisNumber )
    {
        // + + + +  update loop index + + + + +  

        ++iterIndex;

       
        max_absolute_error=SMALL;
        max_relative_error=SMALL;
        coeffIndex=iterIndex-2;

        //  get new residual field 
        forAll (scalarSnapshotsList, snapI)
        {

            // find new WeakGreedy coefficient for snapshot snapI 

            InterpolationCoeffList[snapI].setSize(coeffIndex+1);
            InterpolationCoeffList[snapI][coeffIndex] = MOR::projection(scalarSnapshotsList[snapI], BasisFunctions[coeffIndex]);

            GeometricField<scalar, fvPatchField, volMesh> FieldEstimate
            (
                IOobject
                (
                    scalarSnapshotsList[snapI].name() + "_WeakGreedyEstimate",
                    scalarSnapshotsList[snapI].time().timeName(),
                    scalarSnapshotsList[snapI].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                scalarSnapshotsList[snapI].mesh(),
                dimensioned<scalar>
                (
                    "zero",
                    scalarSnapshotsList[snapI].dimensions(),
                    pTraits<scalar>::zero
                )
            );


            forAll (InterpolationCoeffList[snapI], mfI)
            {

                FieldEstimate+=
                    InterpolationCoeffList[snapI][mfI]*BasisFunctions[mfI];
            }

            GeometricField<scalar, fvPatchField, volMesh> residual
            (
                (scalarSnapshotsList[snapI] - FieldEstimate)
            );


            scalar absolute_interpolation_error = MOR::L2norm(residual);
            snapshotsMeasure = MOR::L2norm(scalarSnapshotsList[snapI]);
            if (snapshotsMeasure < SMALL)
            {
                snapshotsMeasure=SMALL;
            }

            if (max_absolute_error < absolute_interpolation_error)
            {
                max_absolute_error = absolute_interpolation_error;
                generatingFunction = snapI;
            }

            scalar relative_interpolation_error = absolute_interpolation_error/snapshotsMeasure;

            if (max_relative_error < relative_interpolation_error)
            {
                max_relative_error=relative_interpolation_error;
            }
        }

        // update errorLists 
        absoluteErrorList.append(max_absolute_error);
        relativeErrorList.append(max_relative_error);

        // print info 

        Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
        Info<<"maximum absolute error in L2 norm:     " <<  max_absolute_error << endl;
        Info<<"maximum relative error in L2 norm:     " <<  max_relative_error <<endl;
        
        #include "add_BasisFunction.H"
        
    }

    
