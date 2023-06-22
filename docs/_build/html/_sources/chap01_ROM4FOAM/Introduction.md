# Introduction

This library is conceived to be a collection of solvers to be used in OpenFOAM-v6.

Three algorithms have been implemented (with their extensions):

1. Proper Orthogonal Decomposition (POD): with projection for the online phase
2. Empirical Interpolation Method (EIM)
3. Generalised Empirical Interpolation Method (GEIM)
4. Parameterised-Background Data-Weak (PBDW) formulation

Moreover, some useful routines have been implemented in the src/MOR/MOR folder to compute the following:

- Scalar product in $L^2$ for scalar and vector fields:
\begin{equation*}
\langle \phi, \psi\rangle = \int_\Omega \phi \cdot \psi\,d\Omega \qquad \qquad  
\langle \mathbf{u}, \mathbf{v}\rangle = \int_\Omega  \mathbf{u} \cdot \mathbf{v}\,d\Omega
\end{equation*}

- Norm in $L^2$ for scalar and vector fields:
\begin{equation*}
\|\phi\|_{L^2}^2=\int_\Omega \phi^2\,d\Omega \qquad \qquad 
\|\mathbf{u}\|_{L^2}^2=\int_\Omega \mathbf{u}\cdot \mathbf{u}\,d\Omega
\end{equation*}

- Norm in $H^1$ for scalar and vector fields:
\begin{equation*}
\|\phi\|_{H^1}^2=\int_\Omega \phi^2\,d\Omega  + \int_\Omega \nabla\phi\cdot \nabla \phi\,d\Omega \qquad \qquad 
\|\mathbf{u}\|_{H^1}^2=\int_\Omega \mathbf{u}\cdot \mathbf{u}\,d\Omega  + \int_\Omega \nabla\mathbf{u}: \nabla \mathbf{u}\,d\Omega
\end{equation*}

- Norm in $L^\infty$ for scalar and vector fields:
\begin{equation*}
\|\phi\|_{L^\infty} =\max\limits_\Omega |\phi|\qquad \qquad 
\|\mathbf{u}\|_{L^\infty} =\max\limits_\Omega \|\mathbf{u}\|_2
\end{equation*}

## Installation notes

To compile all the application execute
```console
./Allwmake.sh
``` 
in the terminal and to clean every solver
```console
./Allclean.sh
```
in the terminal. Otherwise, execute *wmake* and/or *wclean* in the folder of the solver.