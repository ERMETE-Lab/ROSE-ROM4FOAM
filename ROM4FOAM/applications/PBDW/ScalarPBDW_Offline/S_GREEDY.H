tensor varianceTensor (1/(2*SensorsVariance), 0, 0, 0,1/(2*SensorsVariance),0,0,0,1/(2*SensorsVariance));

scalar M = 0;

// it is now imposed that two sensors cannot be too close
// the minimum distance is given by the inequality exp(-r^2/2s^2) >= 0.1 
scalar distanceTol = std::sqrt ( 2 * 1 * 2.3 * SensorsVariance ); 
Info<<"\n\nThe tolerance for the distanceChek is "<<distanceTol<<endl;

// First Iteration
scalar valueAtSensorPoint = SMALL;
label CellID = 0;

forAll(BasisFunctions[M], CellI)
    {
        if (mag(valueAtSensorPoint) < mag( BasisFunctions[M][CellI]) )
        {
            valueAtSensorPoint = ( BasisFunctions[M][CellI] );
            CellID=CellI;
        }
    }

SensorsCellsID.append( CellID );
Info<<"\nThe sensor #"<<M+1<<" is located in "<<CellID<<endl;

BasisSensors.setSize(M+1);
BasisSensors.set
(
    M,
    new volScalarField
    (
        IOobject
        (
            fieldName+"PBDW_BasisSensor"+name(M),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<scalar>
        (
            "zero",
            dimless,
            pTraits<scalar>::zero
        )
    )
);

forAll (BasisSensors[M], CellI)
{
    if ((Foam::exp(-(BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]])&(varianceTensor& (BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]] ))) > SMALL))
    {
        BasisSensors[M][CellI] = Foam::exp(-(BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]] )&(varianceTensor& (BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]] )));
    }
    else
    {
        BasisSensors[M][CellI]=0;
    }
    
}
        
if ((fvc::domainIntegrate(BasisSensors[M].internalField())).value() < SMALL)
{
    Info << "Error: location of sensor " << M + 1 << " is likely outside the geometrical domain" << endl;
    abort();
}

BasisSensors[M] /= (fvc::domainIntegrate(BasisSensors[M].internalField()));
// While Loop

bool eigenProblem = true;

while ( M < MaxSensorsNumber )
{   

    if ( !eigenProblem)
    {
        scalar tmpDistance = SMALL;
        label generatingCell = 0;
        forAll(BasisSensors[M], CellI)
        {
            scalar distance = 0;
            for (label m = 0; m <= M; ++m)
            {
                distance += std::sqrt (  (BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[m]]) & (BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[m]]) );
            }
            distance /= (M+1);
            
            if ( tmpDistance < distance )
            {
                tmpDistance = distance;
                generatingCell = CellI;
            }
        }
        
        Info<<generatingCell<<" whose point is "<<BasisSensors[M].mesh().C()[generatingCell]<<endl;
        M = M+1;
        SensorsCellsID.append (generatingCell);
    }
    else
    {
        M = M+1;
        infSupConstant.setSize(M+1);
        scalarRectangularMatrix K (M, BasisNumber);
        scalarRectangularMatrix A(M,M);

        for (label ii = 0; ii < M; ++ii)
        {
            for (label jj = 0; jj < M; ++jj)
            {
                A[ii][jj] = MOR::projection(BasisSensors[ii], BasisSensors[jj]);
            }

            for (label kk = 0; kk < BasisNumber; ++kk)
            {
                K[ii][kk] = MOR::projection(BasisSensors[ii], BasisFunctions[kk]);
            }
        }
        
        scalarRectangularMatrix A_inverse (M,M);
        A_inverse = SVDinv(A);
        scalarSquareMatrix SchurComplement (BasisNumber);
        scalarRectangularMatrix tmpSchurComplement (BasisNumber,BasisNumber);
        tmpSchurComplement = K.T() * ( A_inverse * K );
        for (label ii = 0; ii < BasisNumber; ++ii)
        {
            for (label jj = 0; jj < BasisNumber; ++jj)
            {
                SchurComplement[ii][jj] = tmpSchurComplement[ii][jj];
                if (jj<ii) // the assymetry is very low, this makes sure to have a symmetric matrix
                {
                    SchurComplement[ii][jj] = SchurComplement[jj][ii];
                }
            }
        }
        
        EigenSolver<scalar> eigenSolver(SchurComplement);
        scalar tmpMin = 1e9;
        label labelEigenVector = 0; 
        for (label eigenI = 0; eigenI < BasisNumber; ++eigenI)
        {
            if ( eigenSolver.eigenValue(eigenI) < tmpMin)
            {
                tmpMin = eigenSolver.eigenValue(eigenI);
                labelEigenVector = eigenI;
            }
        }
        infSupConstant[M-1] = tmpMin;
        if ( mag(infSupConstant[M-1]) < 1e-12 )
        {
            infSupConstant[M-1] = 0.0;
        }
        if ( infSupConstant[M-1] < 0 )
        {
            Info<<"ERROR: The inf-sup constant is negative!"<<endl;
            abort();
        }
        else 
        {
            infSupConstant[M-1] = std::sqrt(infSupConstant[M-1]);
        }
        Info<<"The inf-sup constant with "<<M<<" sensors is "<<infSupConstant[M-1]<<endl;
        scalarField generatingEigenVector (BasisNumber);
        generatingEigenVector = eigenSolver.eigenVector(labelEigenVector);
        
        volScalarField  z_min
        (
            IOobject
            (
                "z_min",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensioned<scalar>
            (
                "zero",
                dimless, // BasisFunctions[0].dimensions(),
                pTraits<scalar>::zero
            )
        );

        dimensionedScalar dummyVariableField
        (
            "dummyField",
            BasisFunctions[0].dimensions(),
            scalar(1.0)
        );
        for (label ii = 0; ii < BasisNumber; ++ii)
        {
            z_min += generatingEigenVector[ii] * BasisFunctions[ii] / dummyVariableField;
        }

        volScalarField  projected_z_min
        (
            IOobject
            (
                "z_min_Projected",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensioned<scalar>
            (
                "zero",
                dimless, // BasisFunctions[0].dimensions(),
                pTraits<scalar>::zero
            )
        );
        dimensionedScalar dummyVariable
        ( 
            "dummy", 
            dimensionSet(0,3,0,0,0,0,0), 
            scalar(1.0) 
        ); 
        for (label jj = 0; jj < M; ++jj)
        {
            scalar coeff = MOR::projection(z_min, BasisSensors[jj]);
            projected_z_min += dummyVariable * coeff * BasisSensors[jj];
        }

        valueAtSensorPoint = SMALL;
        CellID = 0;
        forAll ( z_min, CellI )
        {
            bool skipCell = false;
            scalar distanceCheck = 0.0;
            forAll (SensorsCellsID, checkI)
            {
                if (SensorsCellsID[checkI] == CellI)
                {
                    skipCell = true;
                }
                else
                {
                    distanceCheck = std::sqrt ( (z_min.mesh().C()[SensorsCellsID[checkI]] - z_min.mesh().C()[CellI] ) & (z_min.mesh().C()[SensorsCellsID[checkI]] - z_min.mesh().C()[CellI] ) );
                    if ( distanceCheck < distanceTol )
                    {
                        skipCell = true;
                    }
                }
                
            }
            if (skipCell == false)
            {
                if ( mag ( z_min[CellI] -projected_z_min[CellI] ) > valueAtSensorPoint )
                {
                    valueAtSensorPoint = mag ( z_min[CellI] -projected_z_min[CellI] );
                    CellID = CellI;
                }
            }
        }
        SensorsCellsID.append( CellID );
        Info<<"\nThe sensor #"<<M+1<<" is located in "<<CellID<<endl;
    
    }
    BasisSensors.setSize(M+1);
    BasisSensors.set
    (
        M,
        new volScalarField
        (
            IOobject
            (
                fieldName+"PBDW_BasisSensor"+name(M),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>
            (
                "zero",
                dimless,
                pTraits<scalar>::zero
            )
        )
    );

    forAll (BasisSensors[M], CellI)
    {
        if ((Foam::exp(-(BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]])&(varianceTensor& (BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]] ))) > SMALL))
        {
            BasisSensors[M][CellI] = Foam::exp(-(BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]] )&(varianceTensor& (BasisSensors[M].mesh().C()[CellI] - BasisSensors[M].mesh().C()[SensorsCellsID[M]] )));
        }
        else
        {
            BasisSensors[M][CellI]=0;
        }
        
    }
            
    if ((fvc::domainIntegrate(BasisSensors[M].internalField())).value() < SMALL)
    {
        Info << "Error: location of sensor " << M + 1 << " is likely outside the geometrical domain" << endl;
        abort();
    }

    BasisSensors[M] /= (fvc::domainIntegrate(BasisSensors[M].internalField()));
    
}
