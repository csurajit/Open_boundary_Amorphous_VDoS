# Open boundary Amorphous VDoS

## Generating open boundary solid droplets in two and three dimensions that uses LJ interaction between bi-disperse particles 

"main_2d.c" generates zero pressure small temperature (T=0.2) liquid configurations under periodic boundary conditions in two dimensions. It uses subroutines: "intialize_2d.c" that generates a square lattice configuration of bi-disperse particles and associates random velocities to particles corresponding to a temperature where the square lattice will melt. Then we generate equilibrated configurations at T=0.55. The system is cooled to T=0.2 at a small cooling rate of 10^-4. Then NPT simulation is performed to get zero-pressure liquid configurations. It uses other subroutines, "verlet_pbc_2d.c/cell_pbc_2d.c" to generate neighbors list and and "force_pbc_nn_2d.c" to calculate the forces using the neighbor list. 
Solid droplets in two dimensions are generated using three minimization protocols in open boundary conditions: 1. Damping dynamics: "damp_dynamics_2d.c", 2.FIRE minimization: "cut_FIRE_2_OB_2d.c". These two minimizations use two subroutines, verlet_2d.c' and 'force_nn_2d.c' to generate neighbors list and to calculate the forces using neighbor list respectively.  The Damped Dynamics protocol uses a dissipative viscous drag in the absence of a temperature bath, FIRE procedure involves numerically integrating a dynamical equation with variable viscous damping and, additionally, a gradient director. Another minimization protocol that we have used is the Conjugate gradient minimization protocol in "cut_out_CG_2d.c". These minimization protocols start by cutting out liquid droplets from zero pressure liquid configurations. Perform the energy minimization, then generate the Hessian matrix which is then diagonalized using "LAPACKE" to get the eigenvalues and eigenvectors. 

"in.3d_KABLJM_PBC_zero_pressure" is the LAMMPS script we have used to generate large system-size three dimensional zero-pressure liquid configurations.

"main_3d.c" generates zero pressure small temperature (T=0.2) liquid configurations under periodic boundary condition in three dimensions. It uses subroutines: "intialize_3d.c" that generates a qubic lattice configuration of bi-disperse particles and associate random velocities to particles corresponding to a temperature where the square lattice will melt. Then we generate equilibrated configurations at T=0.55. The system was then cool to T=0.2 at a small cooling rate of 10^-4. Then NPT simulation is performed to get zero pressure liquid configurations. It uses other subroutines  "verlet_pbc_3d.c/cell_pbc_3d.c" to generate neighbors list and and "force_pbc_nn_3d.c" to calculate the forces using the neighbor list.
Solid droplets are generated using three minimization protocol in open boundary conditions, "damp_dynamics_3d.c" these use two subroutines, 'verlet_3d.c' and 'force_nn_3d.c' to generate neighbours list and to calculate the forces using neighbour list respectively. This minimization protocols start by cutting out liquid droplets by first reading the zero pressure liquid configurations. Perform the energy minimzation by solving equation of motion in the presence of a damping parameter then generate the hessian matrix which is then diagonalized using "LAPACKE" to get the eigenvalues and eigenvectors.

## Distribution of low frequency modes

"dos_cdos_pr_2d.c" reads the first 100 non-zero eigenvalues and measures the probability density doing logarithm binning. An identical  code for the analysis of VDoS in three dimensional systems "dos_cdos_pr_3d.c".

## Analysis of the two dimensional solids structures generated using different minimization protocols 

"mode_analysis_2d.c" perform the anharmonic stability analysis of the minimum energy structure along the direction specified by the first non-zero mode. It reads the inherent structure and the particle displacement corresponding to the first non-zero mode and performs the measure of "B2, B3 and B4" by taking derivatives of the energy funtional.
"dist_beta_2d.c" reads the values of "B2, B3 and B4" and measure the distribution of "\beta = B3/\sqrt{3 B2 B4}" which quantify the stability of the modes.
"displacement_radial_theta_DD.c" measures the radial and tangential component of displacements incurred by the particles in the minimization protocols at different radial distances from the center of mass of the droplets. 

## Effect of surface on the vibrations of open boundary solids

"stress_dist_2d.c" measure the stress distribution within the solid using coarse-grainning boxes. "stress_dist_3d.c" measure the stress distribution with in three dimensional solid.
"parti_ratio_sur.c" measure the contribution of surface particles in the vibrations. Surface particles are identified from the radial distribution of pressure. The contribution of the particles on the surface of the solid to a particular mode is estimated from the displacement vector of the surface particles in the corresponding vibration. "parti_ratio_sur_3d.c" is for three dimensional systems.
"sur_cont_analysis.c" measure the probablity distribution of the quantity that characterize the surface contribution. "sur_cont_analysis_3d.c" measure the histogram for three dimensional systems.
"dos_diff_sc_region.c" measure the probability distribution of the frequencies for surface and bulk localized vibrations. "dos_diff_sc_region_3d.c" is for three dimensional systems.
