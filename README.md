# Reconstruction-Using-Principal-Component-Analysis
# Getting the functional form of the dependent variable in terms of the independent variables. 
<br>
Author: Ranbir

* This program gives the function of dependent vaiable in terms of independent variable when the array of 
* dependent and the indendent variable is available. 

# Theory and the principle are explained in the following papers:
1. https://arxiv.org/abs/2004.01393
2. https://arxiv.org/abs/2211.13608

# input files:
1. The fortran array of dependent and independent variable, with the error in dependent variable (example: beta_data.dat). 
2. The table of patch points in the parameter space (example: input.dat)

# output files:
1. Eigen-vector, eigen-values of the PCA analysis. 
2. The coefficients of the final functional form of the dependent variable in terms of the independent variable.

# Following are the pre-requisite for the source_code_pcaReconstruction.f90 to run.
1. OpenMP
2. Lapack and Blas
3. PyMC
