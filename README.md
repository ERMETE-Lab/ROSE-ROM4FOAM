# ROM4Foam

This repository collects the codes for Hybrid Data Assimilation and Reduced Order Modelling Techniques for OpenFOAM.

## Installation notes

To compile all the applications, execute `./Allwmake.sh` in the terminal.
To clean every solver `./Allwclean.sh` in the terminal.

Otherwise, change the directory to the desired application to compile and execute `wmake`.

## How to cite
If you use the codes in this repository, please cite the papers:

1. Carolina Introini, Simone Cavalleri, Stefano Lorenzi, Stefano Riva, and An- tonio Cammi. Stabilization of Generalized Empirical Interpolation Method (GEIM) in presence of noise: A novel approach based on Tikhonov regularization. Computer Methods in Applied Mechanics and Engineering, 404:115773, 2023. doi: [https://doi.org/10.1016/j.cma.2022.115773](https://doi.org/10.1016/j.cma.2022.115773).
2. Carolina Introini, Stefano Riva, Stefano Lorenzi, Simone Cavalleri, and Antonio Cammi. Non-intrusive system state reconstruction from indirect measurements: A novel approach based on hybrid data assimilation methods. Annals of Nuclear Energy, 182:109538, 2023. doi: [https://doi.org/10.1016/j.anucene.2022.109538](https://doi.org/10.1016/j.anucene.2022.109538).
3. Stefano Riva, C. Introini, S. Lorenzi, and A. Cammi, “Hybrid data assimilation methods, Part I: Numerical comparison between GEIM and PBDW”, Annals of Nuclear Energy, vol. 190, p. 109864, 2023. doi: [https://doi.org/10.1016/j.anucene.2023.109864](https://doi.org/10.1016/j.anucene.2023.109864).
4. Stefano Riva, C. Introini, S. Lorenzi, and A. Cammi, “Hybrid data assimilation methods, Part II: Application to the DYNASTY experimental facility”, Annals of Nuclear Energy, vol. 190, p. 109863, 2023. doi: [https://doi.org/10.1016/j.anucene.2023.109864](https://doi.org/10.1016/j.anucene.2023.109863).

Here listed in bibtex format:

```{=latex}
@article{INTROINI2023115773,
title = {{Stabilization of Generalized Empirical Interpolation Method (GEIM) in presence of noise: A novel approach based on Tikhonov regularization}},
journal = {Computer Methods in Applied Mechanics and Engineering},
volume = {404},
pages = {115773},
year = {2023},
issn = {0045-7825},
doi = {https://doi.org/10.1016/j.cma.2022.115773},
url = {https://www.sciencedirect.com/science/article/pii/S0045782522007290},
author = {Carolina Introini and Simone Cavalleri and Stefano Lorenzi and Stefano Riva and Antonio Cammi},
keywords = {GEIM, Tikhonov Regularization, Noise stabilization, Model order reduction, Data assimilation},
}

@article{INTROINI2023109538,
title = {{Non-intrusive system state reconstruction from indirect measurements: A novel approach based on Hybrid Data Assimilation methods}},
journal = {Annals of Nuclear Energy},
volume = {182},
pages = {109538},
year = {2023},
issn = {0306-4549},
doi = {https://doi.org/10.1016/j.anucene.2022.109538},
url = {https://www.sciencedirect.com/science/article/pii/S0306454922005680},
author = {Carolina Introini and Stefano Riva and Stefano Lorenzi and Simone Cavalleri and Antonio Cammi},
keywords = {Indirect reconstruction, GEIM, POD-I, Reduced order modelling, Data assimilation},
}

@article{RIVA2023109864,
title = {Hybrid data assimilation methods, Part I: Numerical comparison between GEIM and PBDW},
journal = {Annals of Nuclear Energy},
volume = {190},
pages = {109864},
year = {2023},
issn = {0306-4549},
doi = {https://doi.org/10.1016/j.anucene.2023.109864},
url = {https://www.sciencedirect.com/science/article/pii/S0306454923001834},
author = {Stefano Riva and Carolina Introini and Stefano Lorenzi and Antonio Cammi},
keywords = {GEIM, PBDW, Model order reduction, Data assimilation, Thermal hydraulics},
}

@article{RIVA2023109863,
title = {Hybrid Data Assimilation methods, Part II: Application to the DYNASTY experimental facility},
journal = {Annals of Nuclear Energy},
volume = {190},
pages = {109863},
year = {2023},
issn = {0306-4549},
doi = {https://doi.org/10.1016/j.anucene.2023.109863},
url = {https://www.sciencedirect.com/science/article/pii/S0306454923001822},
author = {Stefano Riva and Carolina Introini and Stefano Lorenzi and Antonio Cammi},
keywords = {GEIM, PBDW, Indirect Reconstruction, Model Order Reduction, Data Assimilation, DYNASTY},
}

```

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
