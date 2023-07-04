# Vectorial Field: Online

To be used only after *VectorialPOD_Offline*: this solver implements a POD projection to vector fields.

## Preparation

The structure of the case study folder is the following (in this example *Folder3* and *Folder4* are the test case folders)

```
>> ./Study_case
	>> /Folder_1  			
	>> /Folder_2
	>> /Folder_3  			
	>> /Folder_4		
	>> /POD_(fieldName)
		>> /0		        		
		>> /system		
			controlDict
			blockMeshDict
			...
			PODsolverDict  <--- Dictionary needed for the input parameters	
		>> /constant
		>> /(fieldName)_POD_Offline_files
```

The *PODsolverDict* must be put inside *./Study_case/POD_(fieldName)/system/*

An example of *PODsolverDict* can be found in *application/POD/VectorialPOD_Online*, which requires the following entries:
```
Online_parameters
{
	field      U;			<---- ScalarField on which POD is performed 
	BasisNumber   20;		<---- number of POD modes to use
	foldersList  ( 
			"Folder_3" 
			"Folder_4" ) ;	<---- List of folder names containig the snapshots to be reconstructed
}
```

## Usage

Inside *./Study_case/POD_(fieldName)* launch 
```bash
VectorialPOD_Online
```
To include folder "0" use 
```bash
VectorialPOD_Online -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
VectorialPOD_Online -region <regionName>
```

## Results

The residual field is defined as the absolute difference between the test snapshot and the reconstruction and it is stored in the snapshot folder, as well.

```
>> ./Study_case
	>> /Folder_1  		  		
	>> /Folder_2
	>> /Folder_3
		>> /0
			UPODreconstruct  <--- U POD reconstruction obtained with BasisNumber basis
			UPODresidual     <--- U POD reconstruction obtained with BasisNumber basis
		>> /1	
			UPODreconstruct
			UPODresidual
		>>  ...			
				
	>> /Folder_4
		>> /0
			UPODreconstruct  <---(fieldName) POD reconstruction obtained with BasisNumber basis
			UPODresidual     <---(fieldName) POD reconstruction obtained with BasisNumber basis
		>> /1	
			UPODreconstruct
			UPODresidual
		>>  ...		
			
	>> /POD_U		
		>> /0		        				
		>> /system			
		>> /constant
		>> /U_POD_Offline_files
		>> /U_POD_Online_files
			maximum_L2_relative_error.txt <---- max L2 absolute error as a function of basis number
			average_L2_relative_error.txt <---- max L2 relative error as a function of basis number
```

The absolute and relative error are computed as
```math
E_N^{L^2} = || \mathbf{u}-\mathbf{u}_{N}^{POD}||_{L^2}\qquad 
\epsilon_N^{L^2} = \frac{|| \mathbf{u}-\mathbf{u}_{N}^{POD}||_{L^2}}{|| \mathbf{u} ||_{L^2}}
```