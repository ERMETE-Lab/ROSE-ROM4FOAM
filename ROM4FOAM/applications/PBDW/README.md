# Parametrized-Background Data-Weak (PBDW)

This folder collects the codes for the Parametrized-Background Data-Weak (PBDW) formulation for OpenFOAM, applied to scalar fields only, divided into the typical Offline and Online phase.

During the Online Phase, random noise can be introduced as an option of the solver, moreover the general formulation (with the regularization parameter) is directly implemented.

In each folder there is a file with the instructions on how to use each solver.

The *Alia* folder contains a solver able to compute the matrices A and K of the PBDW formulation, necessary for the analysis of the *inf-sup* error constant.

## Essential Bibliography
- Maday, Y., Patera, A., Penn, J., and Yano, M. (2014). A parameterized- background data-weak approach to variational data assimilation: formulation, analysis, and application to acoustics. International Journal for Numerical Methods in Engineering.
- Maday, Y., Patera, A. T., Penn, J. D., and Yano, M. (2015). PBDW state estimation: Noisy observations; configuration-adaptive background spaces; physical interpretations. ESAIM: Proc.
- Maday, Y. and Taddei, T. (2019). Adaptive PBDW approach to state estimation: Noisy observations; user-defined update spaces. SIAM Journal on Scientific Computing.
- Taddei, T. (2016). Model order reduction methods for data assimilation; state estimation and structural health monitoring.
- S. Riva, C. Introini, S. Lorenzi, and A. Cammi, “Hybrid Data Assimilation Methods (Part I): Numerical Comparison between GEIM and PBDW,” SSRN Electronic Journal, November 2022, doi: 10.2139/ssrn.4313614.