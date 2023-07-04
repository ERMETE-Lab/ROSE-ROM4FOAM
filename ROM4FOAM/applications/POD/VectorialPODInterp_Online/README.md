# Vectorial Field: Interpolation Online

To be used only after *ScalarPOD_Offline*: this solver implements a POD interpolation to vectorial fields.

 
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
			PODInterpSolverDict  <--- Dictionary needed for the input parameters	
		>> /constant
		>> /(fieldName)_POD_Offline_files
			(fieldName)Test_alphaCoeffs <--- matrix generated outside this environment. 
```

The $\alpha_n$ coefficients matrix is structered as follows: at fixed row the basis changes, at fixed column the folder of the Test set changes, thus it's $\mathbb{R}^{N\times N_s^{test}}$.
			

The *PODInterpSolverDict* must be put inside "./Study_case/POD_(fieldName)/system/"

An example of *PODInterpSolverDict* can be found in "application/POD/VectorialPODInterp_Online", which requires the following entries:

```
Online_parameters
{
	field      T;                <---- ScalarField on which POD is performed 
	BasisNumber   20;            <---- number of POD modes to use
	foldersList  ( "Folder_3" 
				"Folder_4") ; <---- List of folder names containig the snapshots to be reconstructed
}
```


## Usage

Inside "./Study_case/POD_(fieldName)" launch 
```bash
VectorialPODInterp_Online
```
To include folder "0" use 
```bash
VectorialPODInterp_Online -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
VectorialPODInterp_Online -region <regionName>
```

## Results

The residual field is defined as the absolute difference between the test snapshot and the reconstruction and it is stored in the snapshot folder, as well.

```
>> ./Study_case

	>> /Folder_1  		  		
	>> /Folder_2
	>> /Folder_3
		>> /0
			(fieldName)PODreconstruct  <---(fieldName) POD reconstruction obtained with mfNumber basis
			(fieldName)PODresidual     <---(fieldName) POD reconstruction obtained with mfNumber basis
		>> /1	
			(fieldName)PODreconstruct
			(fieldName)PODresidual
		>>  ...			
				
	>> /Folder_4
		>> /0
			(fieldName)PODreconstruct 
			(fieldName)PODresidual
		>> /1	
			(fieldName)PODreconstruct
			(fieldName)PODresidual
		>>  ...		
			
	>> /POD_T		
		>> /0		        				
		>> /system			
		>> /constant
		>> /T_POD_Offline_files
		>> /T_POD_Online_files
				maximum_L2_relative_error.txt <---- max L2 absolute error as a function of basis number
				average_L2_relative_error.txt <---- max L2 realtive error as a function of basis number
```

The absolute and relative error are computed as
```math
E_N^{L^2} = || T-T_{N}^{POD}||_{L^2}\qquad 
\epsilon_N^{L^2} = \frac{|| T-T_{N}^{POD}||_{L^2}}{|| T ||_{L^2}}
```