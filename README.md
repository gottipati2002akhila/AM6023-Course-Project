# AM6023-Course_Project
This repository contains the files and codes developed as part of the course project for AM6023: Geometry and Mechanics of Materials. The project focuses on the buckling and post-buckling behavior of a floating elastic beam. The analysis involves solving nonlinear equations, generating bifurcation plots, and visualizing buckling modes using Python scripts and Fortran-based files. Continuation algorithms implemented through AUTO software are used for numerical studies.

The repository includes two main Python scripts. The first, buckling.py, solves the equilibrium equations and plots the buckling thresholds and deformation modes for the floating elastica. The second, modes.py, generates detailed visualizations of the buckling modes, capturing both symmetric and antisymmetric configurations.

Additionally, two folders—Floating Elastica Varying P and Floating Elastica Varying Delta—contain the input files and main Fortran code used for continuation studies. In the first folder, the compressive load P is the input parameter, and the output is the corresponding deformation values. In the second folder, the deformation parameter delta is used as the input, and the output is the corresponding compressive load values. Both folders include files named c.floating and floating.f90 that are designed to run within AUTO, providing a seamless workflow for bifurcation analysis.

The outputs of these computations include bifurcation plots that map the relationship between input and output parameters. These results are generated using AUTO, a widely-used software package for continuation and bifurcation analysis of nonlinear systems.

