# Offline 

Offline phase of the Generalised Empirical Interpolation Method Vectorial Treatment applied to the vector field $(T,\mathbf{u},p_{rgh})$ using temperature sensors only.

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
			GEIM-VTsolverDict  <--- Dictionary needed for the input parameters				
	>> /Folder_2
```

The *GEIM-VTsolverDict* must be put inside *./Study_case/Folder_1/system/*

An example of *GEIM-VTsolverDict* can be found in *application/GEIM-VT/GEIM-VT_Offline*, which requires the following entries:
```
Offline_parameters
{
	error      0.001;            <---- relative L_2 error desidered
	MaxSensorsNumber 20;         <---- Max number of Sensors
	SensorsVariance 0.0001;      <---- variance of the Gaussian function describing the sensors
	foldersList ("Folder_1" 
		     "Folder_2") ;   <---- List of folder names containig the snapshots
	SensorsPositions ( (0.1 0.5 0.3)
			   (0.6 0.8 0.9)
			    ...
			   (0.25 0.41 0.9)  
			  )	     <---- List of admissible sensor locations
																	
} 
```

If 'SensorsPositions' is missing, as many sensors as cell centers will be created. 
**However, for big mashes this option may be very RAM consuming !!!**
It has been implemented a way to tacke these problem, that samples every n cell (the maximum size can be changed according to the available capabilities).
	
## Usage

Inside *./Study_case/Folder_1* launch 
```bash
GEIM-VT_Offline
```
To include folder "0" use 
```bash
GEIM-VT_Offline -withZero
```

# Results

The Magic Functions, Magic Sensors and all the other .txt files are saved in a separate folder called *GEIM-VT_s_(SensorsVariance)*, which has the classical OpenFOAM structure.

```
	>> ./Study_case
		>> /Folder_1  				 		
		>> /Folder_2		
		>> /GEIM-VT_s_0.0001		
			>> /0
				TGEIM-VTMagicFunction0
				TGEIM-VTMagicFunction1
				TGEIM-VTMagicFunction2
				...
				p_rghGEIM-VTMagicFunction0
				p_rghGEIM-VTMagicFunction1
				p_rghGEIM-VTMagicFunction2
				...
				UGEIM-VTMagicFunction0
				UGEIM-VTMagicFunction1
				UGEIM-VTMagicFunction2
				...
				TMagicSensor0
				TMagicSensor1
				TMagicSensor2
				...
				  						
			>> /system
			>> /constant
			>> /GEIM-VT_Offline_files
				Chosen_fields.txt
				Lebesgue_constant.txt
				max_relative_error.txt
				GEIM-VT_T_Offline_Coeffs_avarage_values.txt <---- file needed for 'TR-GEIM-VT'
				GEIM-VT_T_Offline_Coeffs_std.txt            <---- file needed for 'TR-GEIM-VT'


```

The absolute and relative error are computed as
```{math}
\epsilon_M = \text{arg max}_{\phi\in\{T, \mathbf{u}, p_{rgh}\}} \frac{|| \phi-\mathcal{I}_M[\phi]||_{L^\infty}}{||\phi||_{L^2}}
```