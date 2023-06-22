# PBDW - Parameterised-Background Data-Weak formulation

The [Parameterised-Background Data-Weak](https://onlinelibrary.wiley.com/doi/10.1002/nme.4747) (PBDW) was introduced in {cite}`MadayPBDW` as a practical algorithm to general variational data assimilation 
\begin{equation*}
u^\star_\xi = \text{arg min}_{u\in\mathcal{U}}\;\xi ||u-u^{bk}||_{L^2(\Omega)}^2+\frac{1}{M}\sum_{m=1}^M\left(v_m(u)-y_m^{obs}\right)^2
\end{equation*}

The repository implements the algorithm in OpenFOAM, applied to scalar fields only, the details of the implemented version of the formulation can be found in {cite}`Riva2023_Part1`.

There are 2 folders containing the offline and online phase of the algorithm.

- ScalarPBDW_Offline
- ScalarPBDW_Online