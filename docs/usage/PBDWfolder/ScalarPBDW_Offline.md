# Scalar Field: Offline 

Offline phase of the Parametrized-Backward Data-Weak (PBDW) formulation applied to scalar fields.

- The reduced space can be built using two different algorihtm: WeakGreedy and GEIM.
- The selection of the sensors to be used for the measurements is hidden into the GEIM algorithm, if this algorithm is used to build the reduced space, otherwise the SGREEDY algorithm is used.

## Preparation

The structure of the case study folder is the following

```
>> ./Study_case
	>> /Folder_1
		>> /0
		>> ...
		>> /constant
		>> /system
			controlDict
			blockMeshDict
			...
			PBDWsolverDict  <--- Dictionary needed for the input parameters				
	>> /Folder_2
```

The *PBDWsolverDict* must be put inside *./Study_case/Folder_1/system/*

An example of *PBDWsolverDict* can be found in *application/PBDW/ScalarPBDW_Offline*, which requires the following entries:
```
Offline_parameters
{
	field      T; 			<---- ScalarField on which PBDW is performed
	BasisNumber      50;		<---- Size of the reduced space (should be lower than MaxSensorsNumber)
	MaxSensorsNumber 50;		<---- Max number of Sensors
	SensorsVariance 0.0001;		<---- variance of the Gaussian function describing the sensors
	foldersList (
			"Folder_1"
			"Folder_2") ;	<---- List of folder names containig the snapshots
	SensorsPositions ( 
			(0.1 0.5 0.3)
			(0.6 0.8 0.9)
			...
			(0.25 0.41 0.9) 
			 )		<---- List of admissible sensor locations		
} 
```

If 'SensorsPositions' is missing, as many sensors as cell centers will be created.
**However, for big mashes this option may be very RAM consuming !!!**
There has a limitation (based on the RAM capabilities) that samples every $n$ cell the center of the sensors, this is necessary imposed with big meshes.
	
## Usage

Inside *./Study_case/Folder_1* launch 
```bash
ScalarPBDW_Offline
```
To include folder "0" use 
```bash
ScalarPBDW_Offline -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarPBDW_Offline -region <regionName>
```
To select which algorithm to use launch (default WeakGreedy and SGREEDY will be used, alternative GEIM)
```bash
ScalarPBDW_Offline -algoRS <name>
````

## Results

The Basis Functions, Basis Sensors and all the other .txt files are saved in a separate folder called *PBDW_(fieldName)_(algorithmRS)_s_(SensorsVariance)*, which has the classical OpenFOAM structure.

```
>> ./Study_case
	>> /Folder_1
	>> /Folder_2
	>> /PBDW_(fieldName)_(algorithmRS)_s_(SensorsVariance)
		>> /0
			(fieldName)PBDW_BasisFunction0
			(fieldName)PBDW_BasisFunction1
			(fieldName)PBDW_BasisFunction2
			...
			(fieldName)_PBDW_BasisSensor0
			(fieldName)_PBDW_BasisSensor1
			(fieldName)_PBDW_BasisSensor2
			...	
		>> /constant
		>> /system
		>> /(fieldName)_PBDW_Offline_files
			Lebesgue_constant.txt		<---- Only for GEIM
			max_relative_L2_error.txt	<---- Relative reconstruction error of the Train Set
			max_absolute_L2_error.txt	<---- Absolute reconstruction error of the Train Set
			(fieldName)_Coefficients.txt	<---- Reduced basis coefficients
			(fieldName)matrixB.txt		<---- Only for GEIM
			InfSupConstant.txt		<---- Only for SGREEDY
			(fieldName)orthog_matrix.txt	<---- matrix with (zeta_i, zeta_j)_L2 to check the orthonormality of the basis functions
```

The absolute and relative error are computed as
```{math}
E_M = || T-\mathcal{I}_M[T]||_{L^2}\qquad 
\epsilon_M = \frac{|| T-\mathcal{I}_M[T]||_{L^\infty}}{||T||_{L^2}}
```