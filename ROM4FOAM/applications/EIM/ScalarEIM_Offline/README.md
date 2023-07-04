# Description 

Offline phase of the Empirical Interpolation Method applied to scalar fields

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
			EIMsolverDict  <--- Dictionary needed for the input parameters					
	>> /Folder_2
```

The *EIMsolverDict* must be put inside *./Study_case/Folder_1/system/*

An example of *EIMsolverDict* can be found in "application/EIM/ScalarEIM_Offline", which requires the following entries:
```
Offline_parameters
{
	field        T;               <---- ScalarField on which EIM is performed 
	error        0.001;           <---- relative L_infinity error desidered
	maxBasis     20;              <---- Max number of EIM basis functions
	foldersList  (
			"Folder_1" 
			"Folder_2") ; <---- List of folder names containig the snapshots
}
```
## Usage

Inside *./Study_case/Folder_1* launch 
```bash
ScalarEIM_Offline
```
To include folder "0" use 
```bash
ScalarEIM_Offline -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarEIM_Offline -region <regionName>
```

# Results

The Magic Functions, Magic Points and all the other .txt files are saved in a separate folder called *EIM_(fieldName)*, which has the classical OpenFOAM structure.

```
>> ./Study_case
	>> /Folder_1  				 		
	>> /Folder_2		
	>> /EIM_T		
		>> /0		        
				TMagicFunction0
				TMagicFunction1
				TMagicFunction2
				...
									
		>> /system			
		>> /constant
				TMagicPoints
			
		>> /T_EIM_Offline_files
				Lebesgue_constant.txt
				max_relative_L_infinity_error.txt
				max_absolute_L_infinity_error.txt
				MagicParameter.txt
```

The absolute and relative error are computed as
```math
E_M^{\infty} = || \phi-\mathcal{I}_M[\phi]||_{L^\infty}\qquad 
\epsilon_M^\infty = \frac{|| \phi-\mathcal{I}_M[\phi]||_{L^\infty}}{||\phi||_{L^\infty}}
```