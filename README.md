# Computing Transport Coefficients from Molecular Dynamics

This code uses output files from LAMMPS molecular dynamics simulations to compute transport coefficients ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525075/render.png) in electrolyte solutions, where ![equation](http://www.sciweavers.org/upload/Tex2Img_1578524990/render.png).

Here, ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525099/render.png) and ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525118/render.png) are ionic species in solution, ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525141/render.png) is the concentration of species ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525099/render.png), ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525179/render.png) is the velocity of species ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525099/render.png), and ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525212/render.png) is the mass-averaged velocity of all species in the system.

For the purpose of this computation, we assume ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525242/render.png), i.e. there is negligible center of mass velocity of the simulation box. Note that while ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525212/render.png) will impact ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525075/render.png), it will not influence the final ionic conductivity. We also note that ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525179/render.png) is technically an average over all instances of species ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525099/render.png) in the system (e.g. over every cation of a given type in the electrolyte): ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525344/render.png). Thus, we rewrite our original expression for ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525075/render.png) as follows:

![equation](http://www.sciweavers.org/upload/Tex2Img_1578525394/render.png)

![equation](http://www.sciweavers.org/upload/Tex2Img_1578525413/render.png)

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

See the Example notebook for a walkthrough of the computation of ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525075/render.png) for a simulation of LiCl in DMSO. This example also shows how to use the computed ![equation](http://www.sciweavers.org/upload/Tex2Img_1578525075/render.png) to compute ionic conductivity. The output files which should be generated from the example are in the folder example_output.
