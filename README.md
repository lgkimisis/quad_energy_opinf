The repository includes the data and codes for the numerical example of the paper "On the representation of energy-preserving quadratic operators with application to Operator Inference" by L. Gkimisis, I. Pontes, P. Goyal and P. Benner (https://arxiv.org/abs/2503.10824).

%----------------------------------------------------------------------%

(xpos_nonl2D.txt, ypos_nonl2D.txt, usol2_nonl2D.txt) correspond to the x,y coordinate data and the solution of a 2D Bugers' equation on a 1x1 square with periodic boundary conditions

$$\frac{\partial u}{\partial t}=0.2 u  \nabla \cdot u + 0.002 \nabla^2{u}.$$

and initial condition $$u(t = 0) = cos(2\pi x)cos(2\pi y)$$. The timestep is 0.01 s and the simulation time is 4s. The numerical solution was obtained by integrating the discretized using ode45 in MATLAB2024a.

%----------------------------------------------------------------------%

OpInf_2D_Burg.m performs sequential Operator Inference without constraints on the quadratic term, for several reduced dimensions r=[5,10,15,20,25,30]

%----------------------------------------------------------------------%

OpInf_Energy_2D_Burg.m performs sequential Operator Inference for several reduced dimensions r=[5,10,15,20,25,30], with constraints on the quadratic term, such that 

$$\mathbf{x}^\top\mathbf{H} \left( \mathbf{x} \otimes \mathbf{x} \right) =0, \forall \mathbf{x} \in \mathbb{R}^n.$$

by enforcing

$$H=[H_1, \dots H_n],$$ with $$H_i = - H_i^T, \quad \forall i \in [1, \dots n]$$

%----------------------------------------------------------------------%

Comparison_Plots.m creates the plots in the paper.

%----------------------------------------------------------------------%
