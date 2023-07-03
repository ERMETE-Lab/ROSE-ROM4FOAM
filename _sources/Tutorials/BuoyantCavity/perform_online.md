# How to perform the online phase?

The BuoyantCavity problem is parameter-dependent and before performing the online phase some changes in the directory should be made.

The structure of the folder would look like
```
BuoyantCavity
    |---> allrun.py
    |---> BaseCase.py
    |---> TrainSet
    |    |
    |    |---> Case_000_*
    |    |---> Case_001_*
    |    |---> ...
    |    |---> POD_T
    |    |---> POD_U
    |    |---> EIM_T
    |    |---> ...
    |---> TestSet
    |    |
    |    |---> Case_000_*
    |    |---> Case_001_*
    |    |---> ...
```

Change the name of `TrainSet` to `ROM` (or any other name that you like), create then in `ROM` a new folder named `TrainSet` and move all the train snapshots into this: these operations can be done in the terminal with the path set in `BuoyantCavity`
```bash
mv TrainSet ROM
cd ROM
mkdir TrainSet
mv Case_00* TrainSet/.
```
Then, we can also move the test snapshots into `ROM` 
```bash
mv ../TestSet .
```
Therefore, the structure of the repository will look like
```
BuoyantCavity
    |---> allrun.py
    |---> BaseCase.py
    |---> ROM
    |    |
    |    |---> TrainSet
    |    |    |
    |    |    |---> Case_000_*
    |    |    |---> Case_001_*
    |    |    |---> ...
    |    |---> TestSet
    |    |    |
    |    |    |---> Case_000_*
    |    |    |---> Case_001_*
    |    |    |---> ...
    |    |---> POD_T
    |    |---> POD_U
    |    |---> EIM_T
    |    |---> ...
```
## POD
Change directory to `POD_T` and make sure that the file `test_folders.txt` generated with the test snapshots is in `POD_T/system`. The **POD-Online** solver takes as input the test snapshots and projects them into the reduced space, spanned by the POD, and computes the reconstruction error (average and maximum). To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows
```
Online_parameters
{
	field      T;			
	BasisNumber   40;		
	foldersList  (#include "test_folders.txt") ;	
}
```
Now in the terminal (path: `BuoyantCavity/ROM/POD_T`) execute the following command
```
ScalarPOD_Online
```
For the vector field $\mathbf{u}$, change the field name and use the *VectorialPOD_Online* solver.

The output of these solvers is a series of text files into *<field_name>_reconstruction_POD_files* with the average and maximum absolute and relative error.

## EIM
Change directory to `EIM_T` and make sure that the file `test_folders.txt` generated with the test snapshots is in `EIM_T/system`. The **EIM-Online** solver takes as input the test snapshots and use their information at the magic points location to find the interpolant, then computes the reconstruction error (average and maximum). To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows
```
Online_parameters
{
	field      T;		
	mfNumber   40;	
	foldersList  (#include "test_folders.txt");
}
```
Now in the terminal (path: `BuoyantCavity/ROM/EIM_T`) execute the following command
```
ScalarEIM_Online
```
For the vector field $\mathbf{u}$, change the field name and use the *VectorialEIM_Online* solver.

## GEIM - clean data
Change directory to `GEIM_s2_0.0001` and make sure that the file `test_folders.txt` generated with the test snapshots is in `GEIM_s2_0.0001/system`. The **GEIM-Online** solver takes as input the test snapshots and evaluates them at the sensors location through the magic sensors and this information is used to find the interpolant. Then, the solver computes the reconstruction error (average and maximum). To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows
```
Online_parameters
{
	field        T;			
	msNumber     40;		
	foldersList  (#include "test_folders.txt") ;	
}
```
Now in the terminal (path: `BuoyantCavity/ROM/GEIM_s2_0.0001`) execute the following command
```
ScalarGEIM_Online
```
In the offline phase, we have generated two additional folders considering $s^2= 0.0004$ and $s^2=0.0016$, repeat the procedure to get the online reconstruction.

## PBDW - clean data
Change directory to `PBDW_T_WeakGreedy_s_0.0004` and make sure that the file `test_folders.txt` generated with the test snapshots is in `PBDW_T_WeakGreedy_s_0.0004/system`. The **PBDW-Online** solver takes as input the test snapshots and evaluates them at the sensors location through the basis sensors and this information is used to find the reconstruction. Then, the solver computes the reconstruction error (average and maximum). To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows
```
Online_parameters
{
	field T;				 
	MaxSensors 40;			
	BasisNumber 20;				
	sensorsFolder  "PBDW_T_WeakGreedy_s_0.0004";	
	foldersList  (#include "test_folders.txt");		
}
```
The input `sensorsFolder` is used to tell the solver from where the sensors have to be loaded.
Now in the terminal (path: `BuoyantCavity/ROM/PBDW_T_WeakGreedy_s_0.0001`) execute the following command
```
ScalarPBDW_Online
```

## GEIM - noisy data
Change directory to `GEIM_s2_0.0004` and make sure that the file `test_folders.txt` generated with the test snapshots is in `GEIM_s2_0.0004/system`. The main difference with respect to the standard solver is the addiction of noise to the measurement vector $\mathbf{y}\in\mathbb{R}^M$ represented by
\begin{equation*}
    y_m = v_m(u(\cdot;\,\boldsymbol{\mu})) + \epsilon_m
\end{equation*}
where $\epsilon_m$ models random noise as a random variable, i.i. with a zero-mean Gaussian distribution $\sim \mathcal{N}(0,\sigma^2)$. To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows (same as clean data)
```
Online_parameters
{
	field        T;			
	msNumber     40;		
	foldersList  (#include "test_folders.txt") ;	
}
```
Now in the terminal (path: `BuoyantCavity/ROM/GEIM_s2_0.0004`) execute the following command
```
ScalarGEIM_Online -noise 0.01
```
to have $\sigma = 0.01$.

## PBDW - noisy data
Change directory to `PBDW_T_WeakGreedy_s_0.0004` and make sure that the file `test_folders.txt` generated with the test snapshots is in `PBDW_T_WeakGreedy_s_0.0004/system`. The main difference with respect to the standard solver is the addiction of noise to the measurement vector $\mathbf{y}\in\mathbb{R}^M$ represented by
\begin{equation*}
    y_m = v_m(u(\cdot;\,\boldsymbol{\mu})) + \epsilon_m
\end{equation*}
where $\epsilon_m$ models random noise as a random variable, i.i. with a zero-mean Gaussian distribution $\sim \mathcal{N}(0,\sigma^2)$. To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows (same as clean data)
```
Online_parameters
{
	field T;				 
	MaxSensors 40;			
	BasisNumber 20;				
	sensorsFolder  "PBDW_T_WeakGreedy_s_0.0004";	
	foldersList  (#include "test_folders.txt");		
}
```
The input `sensorsFolder` is used to tell the solver from where the sensors have to be loaded.
Now in the terminal (path: `BuoyantCavity/ROM/PBDW_T_WeakGreedy_s_0.0001`) execute the following command
```
ScalarPBDW_Online -noise 0.01
```
to have $\sigma = 0.01$.

## TR-GEIM - noisy data
Change directory to `GEIM_s2_0.0004` and make sure that the file `test_folders.txt` generated with the test snapshots is in `GEIM_s2_0.0004/system`. This solver is needed in presence of noisy data. To execute the solver for the scalar field $T$, the correspondent `dictionary` must be set up as follows
```
Online_parameters
{
	field       T;                 
	msNumber    40;               
	foldersList  (#include "test_folders.txt") ; 
	noise_std   0.01;	      
	N_Repeated_Experiments 5;    
}
in which `noise_std` is the standard deviation $\sigma$ of the random noise and `N_Repeated_Experiments` tells the solver to compute the interpolant more times to have statically robust results.
```
Now in the terminal (path: `BuoyantCavity/ROM/GEIM_s2_0.0004`) execute the following command
```
ScalarTRGEIM
```
A new folder named `T_TR-GEIM_files` is created to store the results.

## GEIM-VT
Change directory to `GEIM-VT_s_0.0004` and make sure that the file `test_folders.txt` generated with the test snapshots is in `GEIM-VT_s_0.0004/system`. To execute the solver, the correspondent `dictionary` must be set up as follows
```
Online_parameters
{
	msNumber 40;         
	foldersList (#include "test_folders.txt") ;                 
} 
```
Now in the terminal (path: `BuoyantCavity/ROM/GEIM-VT_s2_0.0004`) execute the following command
```
GEIM-VT_Online
```