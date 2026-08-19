.. dielcav documentation master file, created by
   sphinx-quickstart on Fri May 29 09:27:12 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

dielcav documentation
=====================

dielcav is a Python module for the analysis of molecular dynamics trajectories of polar liquids, with the aim of extracting their static and dynamic dielectric response. 

It implements the virtual cavity method to obtain a radius and time dependent Kirkwood correlation factor gK(r,t). Provided that the simulation box is large enough this allows faster convergence of the dielectric permittivity than the standard method relying on the total dipole moment of the simulation box. Moreover it makes it possible to disangle self and cross dipolar correlations.

If you have questions or remarks, do not hesitate to send an email to marceau.henot[at]cea.fr

If you use dielcav in your research, please cite the following paper: Hénot J. Chem. Phys. 163, 124501 (2025); doi: 10.1063/5.0289314.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   example
   dielcav
