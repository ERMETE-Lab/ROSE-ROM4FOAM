# Vectorial Field: Offline 

Offline phase of the Empirical Interpolation Method applied to vector fields

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

An example of *EIMsolverDict* can be found in *application/EIM/VectorialEIM_Offline*, which requires the following entries:
```
Offline_parameters
{
	field        U;               <---- VectorField on which EIM is performed 
	error        0.001;           <---- relative L_infinity error desidered
	maxBasis     20;              <---- Max number of EIM basis functions
	foldersList  ( 
			"Folder_1" 
			"Folder_2") ; <---- List of folder names containig the snapshots
}
```
## Usage

Inside "./Study_case/Folder_1" launch 
```bash
VectorialEIM_Offline
```
To include folder "0" use 
```bash
VectorialEIM_Offline -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
VectorialEIM_Offline -region <regionName>
```

# Results

The Magic Functions, Magic Points and all the other .txt files are saved in a separate folder called *EIM_(fieldName)*, which has the classical OpenFOAM structure.

```
>> ./Study_case
	>> /Folder_1  				 		
	>> /Folder_2		
	>> /EIM_U		
		>> /0
			UMagicFunction0
			UMagicFunction1
			UMagicFunction2
			...						
		>> /system			
		>> /constant
			TMagicPoints
		>> /U_EIM_Offline_files
			Lebesgue_constant.txt
			max_relative_L_infinity_error.txt
			max_absolute_L_infinity_error.txt
			MagicParameter.txt
```

The absolute and relative error are computed as
```{math}
E_M^{\infty} = || \mathbf{u}-\mathcal{I}_M[\mathbf{u}]||_{L^\infty}\qquad 
\epsilon_M^\infty = \frac{|| \mathbf{u}-\mathcal{I}_M[\mathbf{u}]||_{L^\infty}}{||\mathbf{u}||_{L^\infty}}
```