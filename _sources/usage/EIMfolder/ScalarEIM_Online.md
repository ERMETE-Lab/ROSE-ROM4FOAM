# Scalar Field: Online 

Online phase of the Empirical Interpolation Method applied to scalar fields

## Preparation

The structure of the case study folder is the following (in this example *Folder3* and *Folder4* are the test case folders)

```
>> ./Study_case
	>> /Folder_1  			
	>> /Folder_2
	>> /Folder_3  			
	>> /Folder_4		
	>> /EIM_(fieldName)
		>> /0		
		>> /constant        		
		>> /system
			controlDict
			blockMeshDict
			...
			EIMsolverDict  <--- Dictionary needed for the input parameters	
		>> /(fieldName)_EIM_Offline_files
```

The *EIMsolverDict* must be put inside *./Study_case/EIM_(fieldName)/system/*

An example of *EIMsolverDict* can be found in *application/EIM/ScalarEIM_Online*, which requires the following entries:
```
Online_parameters
{
	field      T;			<---- ScalarField on which EIM is performed 
	mfNumber   20;			<---- number of EIM magic functions to use
	foldersList  ( 
			"Folder_3" 
			"Folder_4") ;	<---- List of folder names containig the snapshots to be reconstructed
}
```
## Usage

Inside *./Study_case/EIM_(fieldName)* launch 
```bash
ScalarEIM_Online
```
To include folder "0" use 
```bash
ScalarEIM_Online -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarEIM_Online -region <regionName>
```

# Results

The interpolant and the residual field, defined as 
```{math}
r_M = \left| \phi-\mathcal{I}_M[\phi]\,\right|
````
are stored in the correspondent snapshot folders

```
>> ./Study_case
	>> /Folder_1  		  		
	>> /Folder_2
	>> /Folder_3
		>> /0
			TEIMInterpolant  <---(fieldName) EIM interpolant obtained with mfNumber basis
			TEIMresidual     <---(fieldName) EIM residual obtained with mfNumber basis
		>> /1	
			TEIMInterpolant
			TEIMresidual
		>>  ...			
				
	>> /Folder_4
		>> /0
			TEIMInterpolant 
			TEIMresidual
		>> /1	
			TEIMInterpolant
			TEIMresidual
		>> /...		
			
	>> /EIM_T		
		>> /0		        				
		>> /system			
		>> /constant
		>> /T_EIM_Offline_files
		>> /T_EIM_Online_files
			maximum_L2_relative_error.txt <---- max L2 absolute error as a function of basis number
			average_L2_relative_error.txt <---- max L2 realtive error as a function of basis number
```

The absolute and relative error are computed as
```{math}
E_M = || \phi-\mathcal{I}_M[\phi]||_{L^2}\qquad 
\epsilon_M = \frac{|| \phi-\mathcal{I}_M[\phi]||_{L^\infty}}{||\phi||_{L^2}}
```
recalling that the norms are defined as
```{math}
|| \phi ||_{L^2(\Omega)}^2 =\int_\Omega \phi^2\, d\Omega
``` 