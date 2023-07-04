# Description 

Online phase of the Generalised Empirical Interpolation Method applied to scalar fields (with an option to activate synthetic random noise).

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
		>> /system
			controlDict
			blockMeshDict
			...
			GEIMsolverDict  <--- Dictionary needed for the input parameters	
		>> /(fieldName)_GEIM_Offline_files
```

The *GEIMsolverDict* must be put inside *./Study_case/GEIM_(fieldName)_s\_(SensorsVariance)/system/*

An example of *GEIMsolverDict* can be found in *application/GEIM/ScalarGEIM_Online*, which requires the following entries:
```
Online_parameters
{
	field        T;			<---- ScalarField on which GEIM is performed 
	msNumber     20;		<---- number of GEIM magic sensors to use
	foldersList  ( 
			"Folder_3"
			"Folder_4") ;	<---- List of folder names containig the snapshots to be reconstructed
}
```

## Usage

Inside *./Study_case/GEIM_(fieldName)_s\_(SensorsVariance)* launch 
```bash
ScalarGEIM_Online
```
To include folder "0" use 
```bash
ScalarGEIM_Online -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarGEIM_Online -region <regionName>
```
The interpolants using msNumber magic functions will be written in the folder only by activating the option as 
```bash
ScalarGEIM_Online -writeInterpolant
```
Synthetic random noise can be introduced to the data term as follows
```bash
ScalarGEIM_Online -noise <value>
```
where <value> is the std deviation of the noise, assumed zero-mean Gaussian.

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
			TGEIMInterpolant  <---(fieldName) GEIM interpolant obtained with mfNumber basis
			TGEIMresidual     <---(fieldName) GEIM residual obtained with mfNumber basis
		>> /1	
			TGEIMInterpolant
			TGEIMresidual
		>>  ...			
				
	>> /Folder_4
		>> /0
			TGEIMInterpolant
			TGEIMresidual
		>> /1	
			TGEIMInterpolant
			TGEIMresidual
		>> /...		
			
	>> /GEIM_T_s_0.0001		
		>> /0		        				
		>> /system			
		>> /constant
		>> /T_GEIM_Offline_files
		>> /T_GEIM_Online_files
			maximum_L2_absolute_error.txt <---- max L2 absolute error as a function of basis number
			average_L2_absolute_error.txt <---- max L2 absolute error as a function of basis number
			maximum_L2_relative_error.txt <---- max L2 relative error as a function of basis number
			average_L2_relative_error.txt <---- max L2 relative error as a function of basis number
					
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