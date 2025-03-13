The repository includes the data and codes for the paper "On the representation of energy-preserving quadratic operators with application to Operator Inference" by L. Gkimisis, I. Pontes, P. Goyal and P. Benner.
%----------------------------------------------------------------------%
-(xpos_nonl2D.txt, ypos_nonl2D.txt, usol2_nonl2D.txt) correspond to the x,y coordinate data and the solution of a 2D Bugers' equation on a 1x1 square with periodic boundary conditions. 
The timestep is 0.01 s and the simulation time is 4s.

\begin{equation}
    \label{burg}
\frac{\partial u}{\partial t}=0.2 u  \nabla \cdot u + 0.002 \nabla^2{u},
\end{equation}
%----------------------------------------------------------------------%
-OpInf_2D_Burg.m performs sequential Operator Inference without constraints on the quadratic term, for several reduced dimensions r=[5,10,15,20,25,30]
%----------------------------------------------------------------------%
-OpInf_Energy_2D_Burg.m performs sequential Operator Inference for several reduced dimensions r=[5,10,15,20,25,30], with constraints on the quadratic term, such that 
 \begin{equation}
 	\label{energy_pres}
 	\bx(t)^\top\mathbf{H} \left( \bx(t) \otimes \bx(t) \right) =0, \; \forall \bx(t) \in \mathbb{R}^n.
 \end{equation}
%----------------------------------------------------------------------%
-Comparison_Plots.m creates the plots in the paper.
%----------------------------------------------------------------------%
