# Vectorial Field: Offline

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
			PODsolverDict  <--- Dictionary needed for the input parameters				
	>> /Folder_2
```

The *PODsolverDict* must be put inside *./Study_case/Folder_1/system/*

An example of *PODsolverDict* can be found in *application/POD/VectorialPOD_Offline*, which requires the following entries:
```
Offline_parameters
{
	field      U;			<---- VectorialField on which POD is performed 
	accuracy   0.999;		<---- relative energy retained by the POD modes
	maxBasis   20;			<---- Max number of POD modes
	foldersList  ( 
			"Folder_1" 
			"Folder_2") ;	<---- List of folder names containig the snapshots
}
```

## Usage

Inside *./Study_case/Folder_1* launch 
```bash
VectorialPOD_Offline
```
To include folder "0" use 
```bash
VectorialPOD_Offline -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
VectorialPOD_Offline -region <regionName>
```

## Results

POD basis all the other .txt files are saved in a separate folder called *POD_(fieldName)*, which has the classical OpenFOAM structure.

```
>> ./Study_case
	>> /Folder_1  				 		
	>> /Folder_2		
	>> /POD_U		
		>> /0		        
			UPOD0
			UPOD1
			UPOD2
			...
							
		>> /constant				
		>> /system	
		>> /U_POD_Offline_files
			UEigenValues.txt  <---- POD generalised Eigenvalued
			UalphaCoeffs.txt  <---- POD expansion coefficients to be used for mapping in the POD-I
			UL2AbsError.txt   <---- Absolute error in L2
			UL2RelError.txt   <---- Relative error in L2
```

The absolute and relative error are computed as
```{math}
E_N^{L^2} = || \mathbf{u}-\mathbf{u}_{N}^{POD}||_{L^2}\qquad 
\epsilon_N^{L^2} = \frac{|| \mathbf{u}-\mathbf{u}_{N}^{POD}||_{L^2}}{|| \mathbf{u} ||_{L^2}}
```