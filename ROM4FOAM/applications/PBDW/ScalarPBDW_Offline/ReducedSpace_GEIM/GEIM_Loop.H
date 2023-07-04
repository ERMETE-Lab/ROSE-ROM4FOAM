
    // ++++++++++++++++++++++++++++++++++++++ GEIM First Iteration +++++++++++++++++++++++++++++++ 
    
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
    Info <<"maximum absolute interpolation error in L2 norm:     " << setprecision(14) << max_L2norm <<endl;
    Info <<"maximum relative interpolation error in L2 norm:     " << max_L2norm/max_L2norm  <<endl;

    //update errorLists
    absoluteErrorList.append(max_L2norm);
    relativeErrorList.append(max_L2norm/max_L2norm);

    // find magic functional

    forAll(sensors, sensI)
    {
        if (mag(max_measure) < mag(MOR::projection(scalarSnapshotsList[generatingFunction],sensors[sensI])) )
        {
            max_measure= MOR::projection(scalarSnapshotsList[generatingFunction],sensors[sensI]) ;
            tmp_sensID=sensI;
        }
    }

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
            scalarSnapshotsList[generatingFunction]/max_measure,
            fvPatchField<scalar>::calculatedType()
        )
    );
    MagicSensors.append
    (
        sensors.set
        (
            tmp_sensID,
            NULL
        )
    );

 // ++++++++++++++++++++++++++++++++++++++ GEIM Main Loop +++++++++++++++++++++++++++++++ 
    
    scalar max_relative_error;
    scalar max_absolute_error ;
    InterpolationCoeffList.setSize(scalarSnapshotsList.size());
    scalar snapshotsMeasure = SMALL;
    label coeffIndex=0;

    while ( iterIndex <= MaxSensorsNumber )
    {
        // + + + +  update loop index + + + + +  

        ++iterIndex;

       
        max_absolute_error=SMALL;
        max_relative_error=SMALL;
        max_measure=SMALL;
        volScalarField* tmp_residual_Ptr= new GeometricField<scalar, fvPatchField, volMesh>
        (
            scalarSnapshotsList[0]
        );

        coeffIndex=iterIndex-2;

        // update interpolation Matrix

        updateBmatrix (Bmatrix, BasisFunctions, MagicSensors);

        //  get new residual field 
        forAll (scalarSnapshotsList, snapI)
        {

            // find new GEIM coefficient for snapshot snapI 

            InterpolationCoeffList[snapI].setSize(coeffIndex+1);
            InterpolationCoeffList[snapI][coeffIndex] =MOR::projection(scalarSnapshotsList[snapI], MagicSensors[coeffIndex]);

            if ((coeffIndex)!=0)  // fist iteration alpha1 = b1/B11
            {
                for (label ii=0; ii< MagicSensors.size()-1; ii++)
                {
                    InterpolationCoeffList[snapI][coeffIndex]-= (InterpolationCoeffList[snapI][ii]*Bmatrix(coeffIndex, ii));

                }
            }


            GeometricField<scalar, fvPatchField, volMesh> interpolant
            (
                IOobject
                (
                    scalarSnapshotsList[snapI].name() + "GEIMInterpolation",
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

                interpolant+=
                    InterpolationCoeffList[snapI][mfI]*BasisFunctions[mfI];
            }

            GeometricField<scalar, fvPatchField, volMesh> residual
            (
                (scalarSnapshotsList[snapI] - interpolant)
            );


            scalar absolute_interpolation_error =MOR::L2norm(residual);
            snapshotsMeasure = MOR::L2norm(scalarSnapshotsList[snapI]);
            if (snapshotsMeasure < SMALL)
            {
                snapshotsMeasure=SMALL;
            }
            scalar relative_interpolation_error = absolute_interpolation_error/snapshotsMeasure;


            if (max_absolute_error< absolute_interpolation_error)
            {
                tmp_residual_Ptr-> ~GeometricField<scalar, fvPatchField, volMesh>();
                max_absolute_error=absolute_interpolation_error;
                tmp_residual_Ptr = new(tmp_residual_Ptr) GeometricField<scalar, fvPatchField, volMesh>
                (
                    residual
                );
            }

            if (max_relative_error< relative_interpolation_error)
            {
                max_relative_error=relative_interpolation_error;
            }

        }

        if ( iterIndex <= MaxSensorsNumber )
        {
            // find magic functional 

            forAll(sensors, sensI)
            {
                if (sensors(sensI))
                {
                    if (mag(max_measure) < mag(MOR::projection(*tmp_residual_Ptr,sensors[sensI])) )
                    {
                        max_measure= MOR::projection(*tmp_residual_Ptr,sensors[sensI]);
                        tmp_sensID=sensI;
                    }
                }
            }

            // Normalise MagicFUnction 

            *tmp_residual_Ptr/=max_measure;

            // Update Magic Quantites 

            BasisFunctions.append (tmp_residual_Ptr);
            MagicSensors.append
            (
                sensors.set
                (
                    tmp_sensID,
                    NULL
                )
            );

        }

        if (iterIndex > MaxSensorsNumber)
        {
            delete tmp_residual_Ptr;
        }

        // update errorLists 
        absoluteErrorList.append(max_absolute_error);
        relativeErrorList.append(max_relative_error);

        // print info 

        Info<<"\nIteration # "<< iterIndex-1 << ": "<< nl;
        Info<<"maximum absolute interpolation error in L2 norm:     " <<  max_absolute_error << endl;
        Info<<"maximum relative interpolation error in L2 norm:     " <<  max_relative_error <<endl;
    }
