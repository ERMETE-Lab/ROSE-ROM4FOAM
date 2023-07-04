# Scalar Field: Offline

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

An example of *PODsolverDict* can be found in *application/POD/ScalarPOD_Offline*, which requires the following entries:
```
Offline_parameters
{
	field      T;			<---- ScalarField on which POD is performed 
	accuracy   0.999;		<---- relative energy retained by the POD modes
	maxBasis   20;			<---- Max number of POD modes
	foldersList  (
			"Folder_1" 
			"Folder_2") ;	<---- List of folder names containig the snapshots
}
```

## Usage

Inside *./Study_case/Folder_1*"* launch 
```bash
ScalarPOD_Offline
```
To include folder "0" use 
```bash
ScalarPOD_Offline -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarPOD_Offline -region <regionName>
```

## Results

POD basis all the other .txt files are saved in a separate folder called *POD_(fieldName)*, which has the classical OpenFOAM structure.

```
>> ./Study_case
	>> /Folder_1  				 		
	>> /Folder_2		
	>> /POD_T		
		>> /0		        
			TPOD0
			TPOD1
			TPOD2
			...
							
		>> /constant				
		>> /system	
		>> /T_POD_Offline_files
			TEigenValues.txt  <---- POD generalised Eigenvalued
			TalphaCoeffs.txt  <---- POD expansion coefficients to be used for mapping in the POD-I
			TL2AbsError.txt   <---- Absolute error in L2
			TL2RelError.txt   <---- Relative error in L2
```

The absolute and relative error are computed as
```math
E_N^{L^2} = || T-T_{N}^{POD}||_{L^2}\qquad 
\epsilon_N^{L^2} = \frac{|| T-T_{N}^{POD}||_{L^2}}{|| T ||_{L^2}}
```