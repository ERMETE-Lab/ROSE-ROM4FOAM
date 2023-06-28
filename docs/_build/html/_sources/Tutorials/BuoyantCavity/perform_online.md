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

## PBDW - clean data

## (TR-)GEIM - noisy data