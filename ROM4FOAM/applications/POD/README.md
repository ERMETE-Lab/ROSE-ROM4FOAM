# Proper Orthogonal Decomposition (POD)

This folder collects the codes for the Proper Orthogonal Decomposition (POD) for OpenFOAM, applied to scalar and vector fields only, divided into the typical Offline and Online phase.

During the "standard" Online Phase, the true solution is assumed to be known and the reduced basis coefficients are computed through a $L^2$ projection onto the reduced basis $\left(\phi_n\right)_{n=1}^N$.

In each folder there is a file with the instructions on how to use each solver.

The *Alia* folder has some additional solvers for specific purposes, as the calculation of the correlation matrix (to be analyzed in MATLAB or other environments to deal with matrices) and the application of Local-POD on some specific locations.

## Essential Bibliography
- Berkooz, G., Holmes, P., and Lumley, J. L. (1993). The proper orthogonal decomposition in the analysis of turbulent flows. Annual Review of Fluid Mechanics.
- Chatterjee, A. (2000). An introduction to the proper orthogonal decom- position. Current Science.
- Cordier, L. and Bergmann, M. (2008). Proper orthogonal decomposition: an overview.
- Karhunen–loève procedure for gappy data. JOSA A.
- Tezzele, M., Demo, N., Mola, A., and Rozza, G. (2022). An integrated data-driven computational pipeline with model order reduction for industrial and applied mathematics.
