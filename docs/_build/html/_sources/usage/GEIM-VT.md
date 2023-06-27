# GEIM-VT - Generalised Empirical Interpolation Method Vectorial Treatment

The [Generalised Empirical Interpolation Method - Vectorial Treatment](https://link.springer.com/chapter/10.1007/978-88-470-2592-9_13) was presented in {cite}`Maday2015_GEIM`, a summary of the method can be found in {cite}`Introini2023_IR`. In this repository, the algorithm is implemented in OpenFOAM.

There are 3 folders containing the version of the solver for scalar field:

- GEIM-VT_Offline (generation of the magic function and points)
- GEIM-VT_Online (online reconstruction of the field using **synthetic data**)
- TR-GEIM-VT (online reconstruction of the field using **synthetic data** polluted by noise, regularised version)

The last is a novel regularised approach proposed in {cite}`Introini2023_IR` which adopts the Tikhonov regularisation to retrieve stability of the method.