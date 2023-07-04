# Description 

Online phase of the Generalised Empirical Interpolation Method Vectorial Treatment applied to the vector field $(T,\mathbf{u},p_{rgh})$ using temperature sensors only.

Tikhonov Regularizarion allows better reconstruction when signal noise is considered. 

# How to use it

To be used only after *GEIM-VT_Offline*.

## Preparation

The structure of the case study folder is the following

```
>> ./Study_case
	>> /Folder_1  		  		
	>> /Folder_2
	>> /Folder_3  		  		
	>> /Folder_4
	>> /GEIM-VT_s_0.0001
		>> /0		        		
		>> /constant
		>> /GEIM-VT_Offline_files
		>> /system
			controlDict
			blockMeshDict
			...
			GEIM-VTsolverDict  <--- Dictionary needed for the input parameters	
```

The *GEIM-VTsolverDict* must be put inside *./Study_case/GEIM-VT_s_(SensorsVariance)/system/*

An example of *GEIM-VTsolverDict* can be found in *application/GEIM-VT/GEIM-VT_Online*, which requires the following entries:
```
Online_parameters
{
	msNumber    20;               <---- number of GEIM magic sensors to use
	foldersList  ( 
			"Folder_3" 
			"Folder_4") ; <---- List of folder names containig the
										snapshots to be reconstructed
	noise_std   0.001;	      <---- noise Gaussian standard deviation
	N_Repeated_Experiments 10;    <---- Number of repeated "experiments" needed in order to obtain statistically relevant average reconstruction errors
}
```

## Usage

Inside *./Study_case/GEIM-VT_s_(SensorsVariance)* launch 
```bash
TR-GEIM-VT
```
To include folder "0" use 
```bash
TR-GEIM-VT -withZero
```

# Results

The interpolant and the residual field, defined as 
```math
r_M = \left| \phi-\mathcal{I}_M[\phi]\,\right|
````
are stored in the correspondent snapshot folders.

```
>> ./Study_case
	>> /Folder_1  		  		
	>> /Folder_2
	>> /Folder_3
		>> /0
			TGEIM-VTInterpolant	<--- T GEIM-VT interpolant obtained with msNumber basis
			TGEIM-VTresidual	<--- T GEIM-VT residual obtained with msNumber basis
			p_rghGEIM-VTInterpolant	<--- p_rgh GEIM-VT interpolant obtained with msNumber basis
			p_grhGEIM-VTresidual	<--- p_rgh GEIM-VT residual obtained with msNumber basis
			UGEIM-VTInterpolant	<--- U GEIM-VT interpolant obtained with msNumber basis
			UGEIM-VTresidual	<--- U GEIM-VT residual obtained with msNumber basis
		>> /1	
			TGEIM-VTInterpolant  
			TGEIM-VTresidual    
			p_rghGEIM-VTInterpolant 
			p_grhGEIM-VTresidual   
			UGEIM-VTInterpolant 
			UGEIM-VTresidual    
		>> /...			  		
	>> /Folder_4
		>> /0
			TGEIM-VTInterpolant  
			TGEIM-VTresidual  
			...
		>> /constant
		>> /system
	>> /GEIM-VT_s_0.0001		
		>> /0		   
		>> /constant
		>> /system	
		>> /GEIM-VT_Offline_files
		>> /GEIM-VT_Online_files
			T_maximum_L2_relative_error.txt
			T_average_L2_relative_error.txt 
			p_rgh_maximum_L2_relative_error.txt
			p_rgh_average_L2_relative_error.txt 
			U_maximum_L2_relative_error.txt
			U_average_L2_relative_error.txt
```

The absolute and relative error are computed as
```math
\epsilon_M = \frac{|| \phi-\mathcal{I}_M[\phi]||_{L^\infty}}{||\phi||_{L^2}}\qquad \qquad {\phi\in\{T, \mathbf{u}, p_{rgh}\}}
```