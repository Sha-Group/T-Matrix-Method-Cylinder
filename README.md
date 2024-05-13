# T-Matrix-Method-Cylinder
This Matlab code calculates the scattering field of multiple multilayer 2-D plasma cylinders using T-matrix method.

![License](https://img.shields.io/badge/license-GPL3.0-orange)
[![Visits Badge](https://badges.strrl.dev/visits/Shuai-Yuan-1997/T-Matrix-Method-Cylinder)](https://github.com/Shuai-Yuan-1997/T-Matrix-Method-Cylinder)
## Publication
S. S. A. Yuan, Z. H. Lin, L. -B. Lv, S. -J. Hao and W. E. I. Sha, “[Investigating the scattering characteristics of artificial field-aligned irregularities based on T-Matrix algorithm](https://ieeexplore.ieee.org/document/10058168),” IEEE J. Multiscale Multiphys. Comput. Tech., vol. 8, pp. 147-157, 2023.
## MATLAB files
T_matrix_for_Multiple_Cylinders.m: The T-matrix of a single multilayer cylinder is calculated by generalized reflection coefficient method, then a T-matrix equation is formed and solved for the multiple scattering problem. This method can also be used for approximately estimating the RCS of long and thin 3-D cylinders.  

translation_matrix.m: Translation matrix from Hankel to Bessel function

translation_matrix2.m: Translation matrix from Hankel to Hankel function

Notice: Large number of Mie series is needed for larger-scale problem. T-matrix method is suitable for many compact (~10 wavelength) scatterers.
## Other information
Repository Author: Shuai S. A. Yuan

Email: shuaiyuan1997@zju.edu.cn
