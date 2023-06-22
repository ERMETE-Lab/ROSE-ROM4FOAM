# POD - Proper Orthogonal Decomposition

The Proper Orthogonal Decomposition (POD) is considered as the state-of-the-art in Reduced Order Modelling {cite}`MOR_2020Book`, especially in fluid-dynamics and in nuclear reactors applications. The algorithm is based on the Singular Value Decomposition of the snapshots matrix, which can be linked to the correlation matrix $C\in\mathbb{R}^{N_s\times N_s}$:

\begin{equation*}
C_{nm} = \int_\Omega u_n\cdot u_m \, d\Omega\qquad n,m = 1,\dots, N_s
\end{equation*}

and its eigenvalue problem $C\lambda_n = \lambda_n \boldsymbol{\eta}_n $. The POD modes are then defined with the following

\begin{equation*}
\phi_n(\mathbf{x})= \frac{1}{\sqrt{\lambda_n}}\sum_{i=1}^{N_s} \eta_{n,i} u_i(\mathbf{x}) \qquad n = 1, \dots, N
\end{equation*}

which provides also the orthonormality of the modes with respect to the inner product in $L^2$.

The online phase consists in two different version of the reconstruction
\begin{equation*}
u\simeq \sum_{n=1}^N \alpha_n\,\phi_n\qquad \qquad \alpha_n = \int_\Omega u\,\phi_n\,d\Omega
\end{equation*}
the reduced coefficients $\alpha_n$ are computed by projection given some test snapshots or by interpolation of the coefficients through suitable maps $\alpha_n = \mathcal{F}(\alpha_{n,train})$, this version is known as POD-I (POD with Interpolation).

There are 6 folders containing the version of the solver for scalar and vector field, divided into offline (generation of the modes) and online (reconstruction of the field).

- ScalarPOD_Offline
- ScalarPOD_Online
- VectorialPOD_Offline
- VectorialPOD_Online
- ScalarPODInterp_Online
- VectorialPODInterp_Online