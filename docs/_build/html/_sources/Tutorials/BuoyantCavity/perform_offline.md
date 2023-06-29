# How to perform the offline phase?

Once the snapshots have been generated, the ROM can be executed to find the basis functions and/or the basis sensors.

At first, open the terminal and set yourself in `Tutorials/BuoyantCavity/TrainSet` then (if you have used the same intervals for $Re$ and $Ri$ of the previous page) copy the folders' list into the first folder, e.g.
```bash
cp train_folders.txt Case_000_Re15.00_Ri0.20/system/.
```
To execute the offline algorithms, change directory to `Case_000_Re15.00_Ri0.20`. In this tutorial we are going to learn how to perform POD, EIM, GEIM and PBDW onto the scalar -$T, p$- and vector (when possible) -$\mathbf{u}$- fields.

## POD
This algorithm is based on the Singular Value Decomposition (SVD) and it can be performed both for scalar and vector fields. The first thing we need to set up is the **dictionary** of the solver for POD as follows for $T$
```
Offline_parameters
{
    field    T;
    accuracy      0.99999999999;
    maxBasis 50;
    foldersList  (#include "train_folders.txt") ;
} 
```
We are now ready to launch the POD solver for scalar field, by simply typing in the terminal
```bash
ScalarPOD_Offline
```
The solver creates a folder, named `POD_T` (in this case) containing the POD modes and some text files with eigenvalues, the train error and the modal coefficients. For further details, see the specific `README.md` file for this solver. For the pressure $p$ it is sufficient to write *p* instead of *T* in the dict.

For vector fields, the dictionary should be changed as
```
Offline_parameters
{
    field    U;
    accuracy      0.99999999999;
    maxBasis 50;
    foldersList  (#include "train_folders.txt") ;
} 
```
then, the vectorial solver must be launched
```bash
VectorPOD_Offline
```

## EIM
This algorithm can be performed both for scalar and vector fields. The first thing we need to set up is the **dictionary** of the solver for EIM as follows for $T$
```
Offline_parameters
{
    field    T;
    error      1e-5;
    maxBasis 50;
    foldersList (#include "train_folders.txt");
} 
```
We are now ready to launch the EIM solver for scalar field, by simply typing in the terminal
```bash
ScalarEIM_Offline
```
The solver creates a folder, named `EIM_T` (in this case) containing the EIM magic function and some text files with the train error, the modal coefficients, the Lebesgue constant $\Lambda_M$ and the magic points (in `constant`). For further details, see the specific `README.md` file for this solver. For the pressure $p$ it is sufficient to write *p* instead of *T* in the dict.

For vector fields, the dictionary should be changed as
```
Offline_parameters
{
    field    U;
    error      1e-5;
    maxBasis 50;
    foldersList (#include "train_folders.txt");
} 
```
then, the vectorial solver must be launched
```bash
VectorEIM_Offline
```

## GEIM
This algorithm can be performed for scalar fields only. 

```{note}
The extension to vector fields is a matter of future studies and development, see {cite}`Carolina_Tesi`.
```

The first thing we need to set up is the **dictionary** of the solver for GEIM as follows for $T$
```
Offline_parameters
{
    field      T;
    error      1e-8;
    MaxSensorsNumber 50;
    SensorsVariance 0.0001;
    foldersList (#include "train_folders.txt") ;
    //SensorsPositions
} 
```
```{note}
The input `SensorsVariance` is the parameter $s^2$ of the Gaussian kernel centred in $\mathbf{x}_k$ of the sensors, i.e.
\begin{equation*}
g_k = g(\mathbf{x}-\mathbf{x}_k; s) =\displaystyle \frac{e^{-\frac{\|{\mathbf{x}-\mathbf{x}_k}\|_2^2}{2s^2}}}{\displaystyle \int_\Omega e^{-\frac{\|{\mathbf{x}-\mathbf{x}_k}\|_2^2}{2s^2}}\, d\Omega}
\end{equation*}
```
The input `SensorsPositions` is optional, if nothing is entered all the points in the mesh are used for the centers $\mathbf{x}_k$, otherwise only some locations are considered: the input would be as an example
```
SensorsPositions (
    (1., 2., 0.5)
    (1., 3., 0.5)
    ...
)
```
We are now ready to launch the GEIM solver for scalar field, by simply typing in the terminal
```bash
ScalarGEIM_Offline
```
The solver creates a folder, named `GEIM_T` (in this case) containing the GEIM magic function and magic sensors, in addiction to some text files with the train error, the modal coefficients, their average and standard deviation (for TR-GEIM), the Lebesgue constant $\Lambda_M$ and the GEIM matrix $B$. For further details, see the specific `README.md` file for this solver. For the pressure $p$ it is sufficient to write *p* instead of *T* in the dict.

## PBDW
This algorithm can be performed for scalar fields only. 

The first thing we need to set up is the **dictionary** of the solver for PBDW as follows for $T$
```
Offline_parameters
{
    field      T;
    BasisNumber   50;       // maximum size of the reduced basis
    foldersList (#include "train_folders.txt") ;
    
    // The following are specific for the sensors;
    MaxSensorsNumber 60;    
    SensorsVariance 0.0004;
    /*
    SensorsPositions 
    (
		(0 3.09 0.29) 
	);
    */
} 
```
The PBDW is a general framework able to incorporate different methods to generate the basis functions and basis sensors. The default option is the Weak-Greedy coupled with SGREEDY.

For the sensors, see the remarks made for GEIM.

We are now ready to launch the PBDW solver for scalar field, by simply typing in the terminal
```bash
ScalarPBDW_Offline
```
The solver creates a folder, named `PBDW_T_WeakGreedy_s_0.0004` (in this case) containing the PBDW basis function and basis sensors, in addiction to some text files with the train error, the modal coefficients, output of the sensor placement algorithms. For further details, see the specific `README.md` file for this solver. For the pressure $p$ it is sufficient to write *p* instead of *T* in the dict.

## GEIM-VT
This algorithm can be used for indirect reconstruction purposes, estimating the fields (velocity $\mathbf{u}$ and pressure $p$ for this implementation) from the measuremnts of only one scalar field (temperature $T$ in this case). 

The first thing we need to set up is the **dictionary** of the solver for GEIM as follows for $T$
```
Offline_parameters
{
	error      1e-8;            
	MaxSensorsNumber 50;         
	SensorsVariance 0.0004;      

	foldersList (#include "train_folders.txt") ;         
	// SensorsPositions ( (0.1 0.5 0.3) (0.6 0.8 0.9) (0.25 0.41 0.9) );           
} 
```
```{note}
The input `SensorsVariance` is the parameter $s^2$ of the Gaussian kernel centred in $\mathbf{x}_k$ of the sensors, i.e.
\begin{equation*}
g_k = g(\mathbf{x}-\mathbf{x}_k; s) =\displaystyle \frac{e^{-\frac{\|{\mathbf{x}-\mathbf{x}_k}\|_2^2}{2s^2}}}{\displaystyle \int_\Omega e^{-\frac{\|{\mathbf{x}-\mathbf{x}_k}\|_2^2}{2s^2}}\, d\Omega}
\end{equation*}
```
The input `SensorsPositions` is optional, if nothing is entered all the points in the mesh are used for the centers $\mathbf{x}_k$, otherwise only some locations are considered.

We are now ready to launch the GEIM solver for scalar field, by simply typing in the terminal
```bash
GEIM-VT_Offline
```
The solver creates a folder, named `GEIM-VT_s_0.0004` (in this case) containing the GEIM magic functions for all the fields and magic sensors for all the fields temperature, in addiction to some text files with the train error, the average and standard deviation of the reduced coefficients(for TR-GEIM) and the Lebesgue constant $\Lambda_M$. For further details, see the specific `README.md` file for this solver.