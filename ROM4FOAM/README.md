# ROM4Foam

This repository collects the codes for Hybrid Data Assimilation and Reduced Order Modelling Techniques for OpenFOAM.

## Installation notes

To compile all the applications, it's necessary to have [OpenFOAM-v6](https://openfoam.org/version/6/) installed:

- Clone or download the repository in your machine
- Open the terminal in the repo directory and execute

```bash
./Allwmake.sh
```

To clean all ROM4FOAM solvers installed, execute
```bash
./Allwclean.sh
````

If only a single solver is necessary (the *src/MOR* dependencies must be installed), open the terminal in the directory of the desired solved and execute
```bash
wmake
````
whereas to clean
```bash
wclean
````
---

## Structure of the repository

### *applications* folder
This folder contains several algorithms, most of them divided in the typical Offline-Online decomposition:
- **Empirical Interpolation Method** (EIM), for synthetic data only.
- **Generalized Empirical Interpolation Method** (GEIM), for synthetic data only (polluted or not by random noise). A regularization method is implemented based on the Tikhonov regularization (TR-GEIM)

- **Parametrized-Background Data-Weak** (PBDW) formulation, for synthetic data only (polluted or not by random noise). The reduced space can be built either with WeakGreedy or GEIM, whereas the sensors selection is performed through GEIM itself or SGREEDY.
- **Proper Orthogonal Decomposition** (POD), for synthetic data only. Moreover, the POD with Interpolation (POD-I) is implemented.
- **Vectorial Treatment for GEIM** (GEIM-VT), for synthetic data only (polluted or not by random noise), regularization is implemented through Tikhonov only.

### *src* folder
In the *src* folder, there are some libraries containing useful functions able to:
- Compute norms in $L^2$, $L^\infty$ or $H^1$ or the scalar product in $L^2$ (in **MOR** folder);
- Some functions necessary for the Proper Orthogonal Decomposition (POD), e.g., the calculation of the eigenvalues/eigenvectors (in **POD** folder).

## Issues
In case of problems, send an email to stefano.riva@polimi.it or carolina.introini@polimi.it