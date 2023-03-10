# Stability_Fish_Hydro
This readme file will help replicate the results in the paper titled: "Stability of schooling patterns of a fish pair swimming against a flow".

Directory “/Code/” contains all the codes required to replicate the results of this work.

The programs and functions are detailed below:
## “SampleSimulations.m”: 
This code calculates and plots sample trajectories for specified initial conditions.
## “Find_AllRoots_alph.m”:
This code finds the equilibria of the system of equations at different values of alpha for a given value of lambda.
## “Plot_equil_vs_alph.m”:
This code plots the cross-stream coordinate (xi) of stable/unstable equilibria of the system as a function of the flow parameter alpha. 
## “Find_AllRoots_lam_alph.m”:
This code finds the equilibria of the system of equations at different values of alpha and lambda.
## “Plot_heatmap.m”:
This code plots the heatmap of the stable configurations of the system in alpha-lambda plane.
## “Write_nat_freq.m”:
This code computes the natural frequencies of stable equilibria at different values of alpha-lambda.
## “Plot_nat_freq.m”:
This code plots the contour plots of natural frequencies of stable equilibria in alpha-lambda plane.
## “uniqueroots.m”:
This is a function that finds the unique solutions of the system of equations among all the solutions numerically obtained
## “equationmap_Tdip.m”:
This function checks whether a given state of the system is an equilibrium solution of the system and returns 1 if it is and 0 if not.
## “ODEfive_Tdip.m”:
This function returns the system governing equations.
## “func_jacob_Tdip.m”:
This function returns the Jacobian matrix at a given state of the system in a vectorized form.
## “jac_Tdip.m”:
This function calculates and classifies the eigenvalues of the Jacobian matrix of the system at a given state.
## “plot_lines.m”:
This function plots the stable (green) and unstable (red) equilibria.
## “func_SOLVE5eq_var4_Tdip.m”:
This function returns the system five equations for finding the solutions/equilibria of the system.
