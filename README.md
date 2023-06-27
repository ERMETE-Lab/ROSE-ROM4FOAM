# ROM4Foam

This repository collects the codes for Hybrid Data Assimilation and Reduced Order Modelling Techniques for [OpenFOAM-v6](https://openfoam.org/version/6/).

In each folder, some instructions on how to use the code are reported. Further details in the [docs](https://rose-polimi.github.io/ROM4FOAM/intro.html).

## How to cite
If you use the codes in this repository, please cite the papers:

1. Carolina Introini, Simone Cavalleri, Stefano Lorenzi, Stefano Riva, and An- tonio Cammi. Stabilization of Generalized Empirical Interpolation Method (GEIM) in presence of noise: A novel approach based on Tikhonov regularization. Computer Methods in Applied Mechanics and Engineering, 404:115773, 2023. doi: [https://doi.org/10.1016/j.cma.2022.115773](https://doi.org/10.1016/j.cma.2022.115773).
2. Carolina Introini, Stefano Riva, Stefano Lorenzi, Simone Cavalleri, and Antonio Cammi. Non-intrusive system state reconstruction from indirect measurements: A novel approach based on hybrid data assimilation methods. Annals of Nuclear Energy, 182:109538, 2023. doi: [https://doi.org/10.1016/j.anucene.2022.109538](https://doi.org/10.1016/j.anucene.2022.109538).
3. Stefano Riva, C. Introini, S. Lorenzi, and A. Cammi, “Hybrid data assimilation methods, Part I: Numerical comparison between GEIM and PBDW”, Annals of Nuclear Energy, vol. 190, p. 109864, 2023. doi: [https://doi.org/10.1016/j.anucene.2023.109864](https://doi.org/10.1016/j.anucene.2023.109864).
4. Stefano Riva, C. Introini, S. Lorenzi, and A. Cammi, “Hybrid data assimilation methods, Part II: Application to the DYNASTY experimental facility”, Annals of Nuclear Energy, vol. 190, p. 109863, 2023. doi: [https://doi.org/10.1016/j.anucene.2023.109864](https://doi.org/10.1016/j.anucene.2023.109863).

They are listed in bibtex format in `refs.bib`.

## Installation notes

Clone the repository on your machine (suggested to be in `~/OpenFOAM/<username>-6/applications/.`) and follow the following instructions.

To compile all the applications, execute in the terminal
```{bash}
bash ./Allwmake.sh
````
To clean every solver, execute in the terminal
```{bash}
bash ./Allwclean.sh
````

Otherwise, change the directory to the desired application to compile and execute `wmake`.

## Structure of the repository

### *applications* folder
This folder contains several algorithms, most of them divided in the typical Offline-Online decomposition:
- **Empirical Interpolation Method** (EIM), for synthetic data only.
- **Generalized Empirical Interpolation Method** (GEIM), for synthetic data only (polluted or not by random noise). The Tikhonov regularisation is used (TR-GEIM) to treat noisy measures.

- **Parametrized-Background Data-Weak** (PBDW) formulation, for synthetic data only (polluted or not by random noise). The reduced space can be built either with WeakGreedy or GEIM, whereas the sensors selection is performed through GEIM itself (or SGREEDY?).
- **Proper Orthogonal Decomposition** (POD), for synthetic data only. Moreover, the POD with Interpolation (POD-I) is implemented.
- The **Vectorial Treatment for GEIM** (GEIM-VT), for synthetic data only (polluted or not by random noise), regularization is implemented through Tikhonov only.

### *src* folder
In the *src* folder, there are some libraries containing useful functions able to:
- Compute norms in $L^2$, $L^\infty$ or $H^1$ or the scalar product in $L^2$ (in **MOR** folder);
- Some functions necessary for the Proper Orthogonal Decomposition (POD), e.g., the calculation of the eigenvalues/eigenvectors (in **POD** folder).

## Issues
In case of problems, send an email to stefano.riva@polimi.it or carolina.introini@polimi.it
