# Differentially Heated Cavity

The side-driven differentially heated encloures are characterised by one (or two) vertical wall moving at costant velocity with different temperature imposed. This kind of problem shows simple domain with complex fluid dynamics strongly coupled with the buouyancy effects {cite}`Saha2018`.

Some modelling assumptions are made:
- the flow is considered laminar, two-dimensional and steady-state;
- the properties of the fluid are assumed constant except in the body force term, where the density varies linearly with temperature according to the Boussinesq approximation;
- the fluid inside the cavity is assumed Newtonian, while viscous dissipation effects are negligible;
- the Prandtl number is $Pr=0.71$ (air-filled cavity).

The viscous incompressible flow and the temperature distribution inside the cavity are governed by the Navier–Stokes and the energy equations, respectively. In light of the assumptions mentioned above, the parametrized PDEs expressed in properly non-dimensional form read, denoting with $\Omega=[0,1]^2\subset\mathbb{R}^2$ the domain:

\begin{equation*}
\left\{
\begin{aligned}
        \nabla \cdot \mathbf{u}&=0  \quad &\text{ in } \; \Omega\\
        (\mathbf{u} \cdot \nabla)\mathbf{u} - \nu \Delta \mathbf{u}+ \nabla p - \mathbf{g}\,\beta(T-T_\infty) &=0 & \text{ in }\; \Omega \\
        \mathbf{u} \cdot \nabla T - \alpha \Delta T&= 0\quad& \text{ in } \; \Omega \\
        \mathbf{u}(0,y)=\mathbf{u}(1,y)&= (0,1)  \\
        \mathbf{u}(x,0)=\mathbf{u}(x,1)&= (0,0) \\
        T(0,y)=T_H, \; T(1,y)&=T_C  \\
        \left.\frac{\partial T}{\partial y} \right|_{y=0}= \left.\frac{\partial T}{\partial y} \right|_{y=1}&=0 
\end{aligned}
\right.
\end{equation*}

This kind of problems are governed by some dimensionless groups, in addiction to the Prandtl number (specific for a fluid):
- **Reynolds number** $Re$ describing the ratio between inertia and viscous forces
\begin{equation*}
Re = \frac{U\,L}{\nu}
\end{equation*}
- **Grashof number** $Gr$ describing the ratio between buoyancy and viscous forces
\begin{equation*}
Gr = \frac{g\beta(T_H-T_C)L^3}{\nu^2}
\end{equation*}
- **Richardson number** $Ri$ measures the importance of gravitational effects on the fluid dynamics
\begin{equation*}
Ri = \frac{g\beta(T_H-T_C)L}{U^2} = \frac{Gr}{Re^2}
\end{equation*}

In this work, the parameter space has been sampled in the following domain
\begin{equation*}
Re\in[15,\,150]\qquad \qquad Ri\in[0.2, \,5]
\end{equation*}

The geometrical domain has been discretized using a uniform $128\times 128$ grid and the steady-state solver for buoyant, laminar flow of incompressible fluids *buoyantBoussinesqSimpleFoam* has been used.

## Generation of the snapshots
In the [`ROM4FOAM/Tutorials/BuoyantCavity`](https://github.com/ROSE-Polimi/ROM4FOAM) folder, there are a subfolder `BaseCase` and a script `allrun.py` which can be used to solve this problem for various values of $Re$ and $Ri$; in particular, it's possible to specify which solutions are desired. 

The **Train Set** of snapshots is given by
\begin{equation*}
Re\in[15:5:150]\qquad \text{ and } \qquad Ri\in[0.2:0.4:5]
\end{equation*}
correspondent in Python to
```python
dRe = 5.
dRi = 0.4

Re = np.arange(15,  150+dRe/2, dRe)
Ri = np.arange(0.2,   5+dRi/2, dRi)
```
By executing this python code (requires `numpy, pandas, matplotlib` and `os`) with
```bash
python allrun.py
```
The snapshots are generated into `TrainSet` folder (created by the script itself, if required). Then, execute another python script
```bash
python create_folder_list.py
```
to generate a `train_folders.txt` file containing a list of the names of the folders with the snapshots.