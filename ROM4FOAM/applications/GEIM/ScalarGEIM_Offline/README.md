# Description 

Offline phase of the Generalised Empirical Interpolation Method applied to scalar fields

# How to use it

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
			GEIMsolverDict  <--- Dictionary needed for the input parameters					
	>> /Folder_2
```

The *GEIMsolverDict* must be put inside *./Study_case/Folder_1/system/*

An example of *GEIMsolverDict* can be found in *application/GEIM/ScalarGEIM_Offline*, which requires the following entries:
```
{
	field      T;			<---- ScalarField on which GEIM is performed
	error      0.001;		<---- relative L_2 error desidered
	MaxSensorsNumber 20;		<---- Max number of Sensors
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
It has been implemented a way to tacke these problem, that samples every n cell (the maximum size can be changed according to the available capabilities).
	
## Usage

Inside *./Study_case/Folder_1* launch 
```bash
ScalarGEIM_Offline
```
To include folder "0" use 
```bash
ScalarGEIM_Offline -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarGEIM_Offline -region <regionName>
```

# Results

The Magic Functions, Magic Sensors and all the other .txt files are saved in a separate folder called *GEIM_(fieldName)_s_(SensorsVariance)*, which has the classical OpenFOAM structure.

```
>> ./Study_case
	>> /Folder_1  				 		
	>> /Folder_2		
	>> /GEIM_T_s_0.0001		
		>> /0
			TGEIMMagicFunction0
			TGEIMMagicFunction1
			TGEIMMagicFunction2
			...
			TMagicSensor0
			TMagicSensor1
			TMagicSensor2
			...						
		>> /system			
		>> /constant
		>> /T_GEIM_Offline_files
			Lebesgue_constant.txt
			max_relative_L2_error.txt
			max_absolute_L2_error.txt
			betaCoefficients.txt			<---- file needed in the indirect reconstruction (required for the mapping procedure)
			matrixB.txt				<---- file needed in the indirect reconstruction (required for the mapping procedure)
			GEIM_Offline_Coeffs_avarage_values.txt	<---- file needed for 'ScalarTR-GEIM'
			GEIM_Offline_Coeffs_std.txt		<---- file needed for 'ScalarTR-GEIM'
```

The absolute and relative error are computed as
```math
E_M = || T-\mathcal{I}_M[T]||_{L^2}\qquad 
\epsilon_M = \frac{|| T-\mathcal{I}_M[T]||_{L^\infty}}{||T||_{L^2}}
```
recalling that the norms are defined as
```math
|| T ||_{L^2(\Omega)}^2 =\int_\Omega T^2\, d\Omega
``` 