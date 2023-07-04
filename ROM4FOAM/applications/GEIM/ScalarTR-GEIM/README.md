# Description 

Online phase of the Tikhonov regularized Generalised Empirical Interpolation Method applied to scalar fields.
Tikhonov Regularizarion allows better reconstructions when signal noise is considered.

# How to use it

To be used only after *ScalarGEIM_Offline*

## Preparation

The structure of the case study folder is the following (in this example *Folder3* and *Folder4* are the test case folders)

```
>> ./Study_case
	>> /Folder_1  			
	>> /Folder_2
	>> /Folder_3  			
	>> /Folder_4		
	>> /GEIM_(fieldName)
		>> /0		  
		>> /constant
		>> /(fieldName)_GEIM_Offline_files      		
		>> /system		
			controlDict
			blockMeshDict
			...
			GEIMsolverDict  <--- Dictionary needed for the input parameters	
```

The *GEIMsolverDict* must be put inside *./Study_case/GEIM_(fieldName)_s\_(SensorsVariance)/system/*

An example of *GEIMsolverDict* can be found in *application/GEIM/ScalarGEIM_Online*, which requires the following entries:
```
Online_parameters
{
	field       T;                <---- ScalarField on which GEIM is performed 
	msNumber    20;               <---- number of GEIM magic sensors to use
	foldersList  ( "Folder_3", 
			"Folder_4") ; <---- List of folder names containig the
										snapshots to be reconstructed
	noise_std   0.001;	      <---- noise Gaussian standard deviation
	N_Repeated_Experiments 10;    <---- Number of repeated "experiments" needed in order to obtain statistically relevant average reconstruction errors
}
```

## Usage

Inside *./Study_case/GEIM_(fieldName)_s\_(SensorsVariance)* launch 
```bash
ScalarTR-GEIM
```
To include folder "0" use 
```bash
ScalarTR-GEIM -withZero
```

# Results

The interpolant and the residual field, defined as 
```math
r_M = \left| T-\mathcal{I}_M[T]\,\right|
````
are stored in the correspondent snapshot folders

```
>> ./Study_case
	>> /Folder_1  		  		
	>> /Folder_2
	>> /Folder_3
		>> /0
			T_TR-GEIMInterpolant  <---(fieldName) TR-GEIM interpolant obtained with mfNumber basis
			T_TR-GEIMresidual     <---(fieldName) TR-GEIM residual obtained with mfNumber basis
		>> /1	
			T_TR-GEIMInterpolant
			T_TR-GEIMresidual
		>>  ...			
				
	>> /Folder_4
		>> /0
			T_TR-GEIMInterpolant
			T_TR-GEIMresidual
		>> /1	
			T_TR-GEIMInterpolant
			T_TR-GEIMresidual
		>>  ...		
			
	>> /GEIM_T_s_0.0001		
		>> /0		        				
		>> /system			
		>> /constant
		>> /T_GEIM_Offline_files
		>> /T_GEIM_Online_files
		>> /T_TR-GEIM_files
			average_L2_relative_error_noiseStd_0.001.txt" <---- average L2 realtive error as a function of basis number with noise_std 0.001	
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