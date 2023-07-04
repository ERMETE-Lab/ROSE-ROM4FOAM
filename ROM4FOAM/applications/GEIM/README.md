# Generalized Empirical Interpolation Method (GEIM)

This folder collects the codes for the Generalized Empirical Interpolation Method (GEIM) for OpenFOAM, applied to scalar fields only, divided into the typical Offline and Online phase.

During the Online Phase, random noise can be introduced as an option of the solver, moreover regularization is implemented using TR-GEIM.

In each folder there is a file with the instructions on how to use each solver.

The *Alia* folder is used for extras.

## Essential Bibliography
- Argaud, J.-P., Bouriquet, B., de Caso, F., Gong, H., Maday, Y., and Mula, O. (2018). Sensor placement in nuclear reactors based on the generalized empirical interpolation method. Journal of Computational Physics.
- Argaud, J.-P., Bouriquet, B., Gong, H., Maday, Y., and Mula, O. (2017). Stabilization of (G)EIM in presence of measurement noise: Application to nuclear
reactor physics.
- Maday, Y. and Mula, O. (2013). A generalized empirical interpolation method: Application of reduced basis techniques to data assimilation. Springer INdAM Series.
- Maday, Y., Mula, O., Patera, A. T., and Yano, M. (2015). The generalized empirical interpolation method: Stability theory on hilbert spaces with an application to the stokes equation. Computer Methods in Applied Mechanics and Engineering
- Maday, Y., Mula, O., and Turinici, G. (2016). Convergence analysis of the generalized empirical interpolation method. SIAM Journal on Numerical Analysis. 

### Regularization
- C. Introini, S. Cavalleri, S. Lorenzi, S. Riva, and A. Cammi, “Stabilization of Generalized Empirical Interpolation Method (GEIM) in presence of noise: A novel approach based on Tikhonov regularization,” Computer Methods in Applied Mechanics and Engineering, vol. 404, p. 115773, 2023, doi: 10.1016/j.cma.2022.115773.