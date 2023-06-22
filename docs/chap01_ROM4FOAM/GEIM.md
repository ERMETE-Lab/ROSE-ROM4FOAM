# GEIM - Generalised Empirical Interpolation Method

The [Generalised Empirical Interpolation Method](https://link.springer.com/chapter/10.1007/978-88-470-2592-9_13) was firstly presented in {cite}`Maday2013` and in literature several extension and additional studies are present. In this repository, the algorithm is implemented in OpenFOAM, to scalar field only.

There are 3 folders containing the version of the solver for scalar field:

- ScalarGEIM_Offline (generation of the magic function and points)
- ScalarGEIM_Online (online reconstruction of the field using **synthetic data**)
- ScalarTR-GEIM (online reconstruction of the field using **synthetic data** polluted by noise, regularised version)

The last is a novel regularised approach proposed in {cite}`Introini2023_TRGEIM` which adopts the Tikhonov regularisation to retrieve stability of the method.

Here we report the algorithm for scalar fields (without regularisation).

```{image} ../images/chap1/GEIM-algo.png
:alt: GEIM-algo
:class: bg-primary mb-1
:width: 1000px
:align: center
```