# PBDW - Parameterised-Background Data-Weak formulation

The [Parameterised-Background Data-Weak](https://onlinelibrary.wiley.com/doi/10.1002/nme.4747) (PBDW) was introduced in {cite}`MadayPBDW` as a practical algorithm to general variational data assimilation 
\begin{equation*}
u^\star_\xi = \text{arg }\min\limits_{u\in\mathcal{U}}\;\xi ||u-u^{bk}||_{L^2(\Omega)}^2+\frac{1}{M}\sum_{m=1}^M\left(v_m(u)-y_m^{obs}\right)^2
\end{equation*}

... aggiungere come si arriva a formulazione lineare.


The repository implements the algorithm in OpenFOAM, applied to scalar fields only, the details of the implemented version of the formulation can be found in {cite}`Riva2023_Part1`.

There are 2 folders containing the offline and online phase of the algorithm.

- ScalarPBDW_Offline
- ScalarPBDW_Online

The PBDW is a general framework to combine data and mathematical models approximated through reduced basis techniques, hence it can accomodate different techniques to generate the basis functions and the basis sensors.

In this work, the default option is given by the couple WeakGreedy+SGREEDY, alternatevely the greedy procedure of GEIM is used.

## WeakGreedy algorithm
The rationale behind this algorithm is quite similar to the GEIM one, and the main difference between the two stands in the generation of the basis functions. 

The first iteration starts by selecting the generating function
\begin{equation*}
    u_1(\mathbf{x}) = \text{arg } \max\limits_{u\in\mathcal{U}} \|{u(\mathbf{x};\mu)}\|_{L^2(\Omega)}
\end{equation*}
and the correspondent basis function as
\begin{equation*}
	\zeta_1(\mathbf{x}) = \frac{u_1(\mathbf{x})}{\|{u_1(\mathbf{x})}\|_{L^2(\Omega)}}\qquad\qquad Z_1 = \text{span}\{\zeta_1\}
\end{equation*}
Then the main loop, where $2\leq M \leq M_{max}$, begins: the generating function is selected as the one maximizing the error
\begin{equation*}
	u_M(\mathbf{x}) = \text{arg }\max\limits_{u\in\mathcal{U}} E_{M-1}[u] = \text{arg }\max\limits_{u\in\mathcal{U}} \left\|{u(\mathbf{x};\mu)-\sum_{i = 1}^{M-1}z_i(\mu)\cdot \zeta_i(\vec{x})}\right\|_{L^2(\Omega)}
\end{equation*}
where $z_i = \left(u, \zeta_i\right)_{L^2(\Omega)}$ and $E_{M-1}$ is the reconstruction error at the $M-$th iteration. The generating function is later orthonormalized with respect to the basis, using the Gram-Schmidt procedure and accordingly the reduced space will be defined as
\begin{equation*}
	Z_M = \text{span}\{\zeta_1, \dots, \zeta_M\}
\end{equation*}
$E_{M-1}$ may be replaced by a proper error estimator $\Delta_{M-1}\geq E_{M-1}$, which allows a speed up of the calculations. 
## SGREEDY algorithm