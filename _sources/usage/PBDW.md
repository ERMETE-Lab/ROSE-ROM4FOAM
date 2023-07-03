# PBDW - Parameterised-Background Data-Weak formulation

The [Parameterised-Background Data-Weak](https://onlinelibrary.wiley.com/doi/10.1002/nme.4747) (PBDW) was introduced in {cite}`MadayPBDW` as a practical algorithm to general variational data assimilation 
\begin{equation*}
u^\star_\xi = \text{arg }\min\limits_{u\in\mathcal{U}}\;\xi ||u-u^{bk}||_{L^2(\Omega)}^2+\frac{1}{M}\sum_{m=1}^M\left(v_m(u)-y_m^{obs}\right)^2
\end{equation*}

It has been shown in {cite}`Taddei_phdThesis` that this can be written as a weak formulation to be later converted into a linear system of small dimension. The state estimation can be written as a linear combination in the following way
\begin{equation*}
    u_\xi^\star (\mathbf{x}) = \sum_{m=1}^M \gamma_m\cdot g_m(\mathbf{x})+\sum_{n=1}^N z_n\cdot \zeta_n(\mathbf{x}),
\end{equation*}
in which the first summation represents the correction term related to the measurements, whereas the latter is the part arising from the reduced basis approximation of the snapshots space. The coefficients are the solution of the following linear system
\begin{equation*}
	\left[ 
	\begin{array}{ccc}
		\xi M I+A& & K  \\ & & \\
		K^T & & 0
	\end{array}
	\right] \cdot
	\left[ 
	\begin{array}{c}
		\boldsymbol{\gamma}_\xi^\star \\  \\\mathbf{z}_\xi^\star
	\end{array}
	\right]   =
	\left[ 
	\begin{array}{c}
		\mathbf{y}^{obs} \\  \\\mathbf{0}
	\end{array}
	\right]
\end{equation*}
provided the following definitions: let $A\in\mathbb{R}^{M\times M}$ and $K\in\mathbb{R}^{M\times N}$ matrices, defined as
\begin{equation*}
	\begin{split}
		A_{mm'}&=\left(g_m,\,g_{m'}\right)_{L^2(\Omega)}\qquad \qquad\quad\;\; m,m' = 1, \dots M,\\
		K_{mn}&=\left(g_m,\,\zeta_{n}\right)_{L^2(\Omega)}= v_m(\zeta_n)\qquad m= 1, \dots M\quad n= 1, \dots, N.
	\end{split}
\end{equation*}

The algorithm is implemented in OpenFOAM, applied to scalar fields only, the details of the implemented version of the formulation can be found in {cite}`Riva2023_Part1`.

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
	u_M(\mathbf{x}) = \text{arg }\max\limits_{u\in\mathcal{U}} E_{M-1}[u] = \text{arg }\max\limits_{u\in\mathcal{U}} \left\|{u(\mathbf{x};\mu)-\sum_{i = 1}^{M-1}z_i(\mu)\cdot \zeta_i(\mathbf{x})}\right\|_{L^2(\Omega)}
\end{equation*}
where $z_i = \left(u, \zeta_i\right)_{L^2(\Omega)}$ and $E_{M-1}$ is the reconstruction error at the $M-$th iteration. The generating function is later orthonormalized with respect to the basis, using the Gram-Schmidt procedure and accordingly the reduced space will be defined as
\begin{equation*}
	Z_M = \text{span}\{\zeta_1, \dots, \zeta_M\}
\end{equation*}
$E_{M-1}$ may be replaced by a proper error estimator $\Delta_{M-1}\geq E_{M-1}$, which allows a speed up of the calculations. 

## SGREEDY algorithm
This algorithm maximizes the *inf-sup* constant $\beta_{N,M}$ in a greedy manner {cite}`HAIK2023115868, Maday2015_GEIM`, the main difference with respect to GEIM is that this procedure works also for $M>N$, hence we can place more sensors $M$ than the number of basis function $N$ used; furthermore SGREEDY is equivalent to GEIM are equivalent if $M=N$. The details are reported in following algorithm.

```{image} ../images/chap1/SGREEDY-algo.png
:alt: NRGlogo
:class: bg-primary mb-1
:width: 1000px
:align: center
```

````{prf:theorem} Inf-Sup theorem
:label: inf-sup-theorem
The inf-sup constant $\beta_{N,M}$ is the square root of the minimum eigenvalue of the following problem
\begin{equation*}
    K^TA^{-1}K\mathbf{z}_n = \lambda_n Z\mathbf{z}_n \qquad\qquad n = 1, \dots N
\end{equation*}
where the matrices are defined as
\begin{equation*}
\begin{split}
    A_{mm'}&=\left(g_m,\,g_{m'}\right)_{L^2(\Omega)}\qquad m,m' = 1, \dots M\\
    K_{mn}&=\left(g_m,\,\zeta_{n}\right)_{L^2(\Omega)}\qquad m= 1, \dots M\quad n= 1, \dots, N\\
	Z_{nn'} &= \left(\zeta_n,\zeta_{n'}\right)_{L^2(\Omega)} \qquad n,n' = 1, \dots, N
    \end{split}
\end{equation*}
````