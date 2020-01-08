# Computing Transport Coefficients from Molecular Dynamics

This code uses output files from LAMMPS molecular dynamics simulations to compute transport coefficients $L^{ij}$ in electrolyte solutions, where ![equation](http://www.sciweavers.org/upload/Tex2Img_1578524990/render.png).

Here, $i$ and $j$ are ionic species in solution, $c_i = \frac{N_i}{V}$ is the concentration of species $i$, $\mathbf{v}_i$ is the velocity of species $i$, and $\mathbf{v}$ is the mass-averaged velocity of all species in the system.

For the purpose of this computation, we assume $\mathbf{v} = 0$, i.e. there is negligible center of mass velocity of the simulation box. Note that while $\mathbf{v}$ will impact $L^{ij}$, it will not influence bulk, experimentally measurable transport quantities which can be computed from $L^{ij}$, such as the ionic conductivity. We also note that $\mathbf{v}_i$ is technically an average over all instances of species $i$ in the system (e.g. over every cation of a given type in the electrolyte): $\mathbf{v}_i = \frac{\sum_\alpha^{N_i} v_{i,\alpha}(t)}{N_i}$. Thus, we rewrite our original expression for $L^{ij}$ as follows:

$L^{ij} = \frac{V}{3 k_B T} \int_0^\infty \big<\frac{N_i N_j}{V^2}\frac{\sum_\alpha^{N_i} v_{i,\alpha}(t)}{N_i} \cdot \frac{\sum_\beta^{N_j} v_{j,\beta}(0)}{N_j}\big>$

$L^{ij} = \frac{1}{3 k_B T V} \int_0^\infty \big<\sum_\alpha v_{i,\alpha}(t) \cdot \sum_\beta v_{j,\beta}(0)\big>$

For computational ease, this final form of the equation is the one evaluated in this code.

### Preparation of input files

To run this code, you will need dump files from a LAMMPS simulation outputting the velocities of the ions in solution (one file for each type of ion). These can be prepared using the following lines in the LAMMPS input script:

            # Calculate velocities of anions
            group anion type 2  # creat group containing all anions
            dump anion_dump anion custom 1 v_anion.out id vx vy vz  # dump velocities every time step

            # Calculate velocities of cations
            group cation type 3
            dump cation_dump cation custom 1 v_cation.out id vx vy vz

### Examples

See the Example notebook for a walkthrough of the computation of $L^{ij}$ for a simulation of LiCl in DMSO. This example also shows how to use the computed $L^{ij}$ to compute ionic conductivity. The output files which should be generated from the example are in the folder example_output.
