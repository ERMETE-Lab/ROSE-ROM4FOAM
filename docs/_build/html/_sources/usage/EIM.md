# EIM - Empirical Interpolation Method

The [Empirical Interpolation Method](https://www.sciencedirect.com/science/article/pii/S1631073X04004248) was firstly presented in {cite}`MadayEIM_2006`, this repository implemented in OpenFOAM, to both scalar and vector fields only, extended in {cite}`Silva2021`.

There are 4 folders containing the version of the solver for scalar and vector field, divided into offline (generation of the magic function and points) and online (creation of the synthetic data and field estimation).

- ScalarEIM_Offline
- ScalarEIM_Online
- VectorialEIM_Offline
- VectorialEIM_Offline

Here we report the algorithm for scalar fields.

```{image} ../images/chap1/EIM-algo.png
:alt: NRGlogo
:class: bg-primary mb-1
:width: 1000px
:align: center
```