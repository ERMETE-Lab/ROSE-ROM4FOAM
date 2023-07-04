# Description 

Online phase of the Parametrized-Backward Data-Weak formulation applied to scalar fields.
There is an option to activate random noise, the regularization can be activating by using the suitable option.

# How to use it

To be used only after *ScalarPBDW_Offline*

## Preparation

The structure of the case study folder is the following (in this example *Folder3* and *Folder4* are the test case folders)

```
>> ./Study_case
	>> /Folder_1  	
	>> /Folder_2
	>> /Folder_3  	
	>> /Folder_4
	>> /PBDW_(fieldName)_(algorithmRS)_s_(sensorsVariance)
		>> /0
		>> /constant
		>> /T_PBDW_Online_files
		>> /system
			controlDict
			blockMeshDict
			...
			PBDWsolverDict  <--- Dictionary needed for the input parameters
```

The *PBDWsolverDict* must be put inside *./Study_case/PBDW_(fieldName)_(algorithmRS)_s_(SensorsVariance)/system/*

An example of *PBDWsolverDict* can be found in *application/PBDW/ScalarPBDW_Online*, which requires the following entries:
```
Online_parameters
{
	field T;				<---- ScalarField on which PBDW is performed 
	msNuMaxSensorsber 20;			<---- number of PBDW magic sensors to use
	BasisNumber 10;				<---- number of PBDW basis functions to use
	sensorsFolder  "sensorsFolderName";	<----- name of the folder in which the sensors are saved
	foldersList  ( 
			"Folder_3"
			"Folder_4");		<---- List of folder names containig the snapshots to be reconstructed
}
```

## Usage

Inside *./Study_case/PBDW_(fieldName)_(algorithmRS)_s_(SensorsVariance)* launch 
```bash
ScalarPBDW_Online
```
To include folder "0" use 
```bash
ScalarPBDW_Online -withZero
```
To perform on a specified region (for multi-region cases) use 
```bash
ScalarPBDW_Online -region <regionName>
```
Synthetic random noise can be introduced to the data term as follows
```bash
ScalarPBDW_Online -noise <value>
```
where <value> is the std deviation of the noise, assumed zero-mean Gaussian.
To include regularization launch
```bash
ScalarPBDW_Online -reg <value>
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
			TPBDW_estimation   <---(fieldName) PBDW estimate obtained with MaxSensors and BasisNumber
			TPBDW_residual     <---(fieldName) PBDW residual obtained with MaxSensors and BasisNumber
		>> /1
			TPBDW_estimation
			TPBDW_residual
		>>  /...
	>> /Folder_4
		>> /0
			TPBDW_estimation 
			TPBDW_residual
		>> /1
			TPBDW_estimation
			TPBDW_residual
		>>  /...
	>> /PBDW_T_algorithmRS_s_sensorsVariance
		>> /0
		>> /constant
		>> /system
		>> /T_PBDW_Offline_files
		>> /T_PBDW_Online_files
			maximum_L2_relative_error.txt	<---- max L2 relative error as a function of sensor number
			average_L2_relative_error.txt	<---- ave L2 relative error as a function of sensor number
			maximum_L2_absolute_error.txt	<---- max L2 absolute error as a function of sensor number
			average_L2_absolute_error.txt	<---- ave L2 absolute error as a function of sensor number
			TmatrixA.txt			<---- matrix A_ij = (g_i, g_j) of the PBDW
			TmatrixK.txt			<---- matrix K_ij = (g_i, zeta_j) of the PBDW			
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