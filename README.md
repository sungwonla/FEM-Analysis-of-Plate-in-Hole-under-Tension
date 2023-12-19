# FEM-Analysis-of-Plate-in-Hole-under-Tension
FEM code written in MATLAB for analysis and visualization of a hole in a semi-infinite plate, subject to tension
(edited from Prof. Allan Bower's code) (ENGN 2340 Final Project)

Main template: fem_template_hole.m 

Running fem_template_hole.m will call the mesh file "mesh_hole.m", and will produce outputs for visualization of the FEM solution of the plate with a hole, as well as the stress concentration factor and the L2 energy norm. The mesh is generated with quarter symmetry to save computational costs. To change the mesh size, edit the variable "meshn" in line 7 of "fem_template_hole.m". 
