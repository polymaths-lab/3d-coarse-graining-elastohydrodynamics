# 3d-coarse-graining-elastohydrodynamics

Code related to "The three-dimensional coarse-graining formulation of interacting elastohydrodynamic filaments and multi-body microhydrodynamics"

Publication:
https://royalsocietypublishing.org/doi/10.1098/rsif.2023.0021

Running "main.m" will solve for a relaxating filament that is also driven by a time-dependent intrinsic curvature. 

Can edit "get_params_relaxing_filament.m" to change the filament parameters, for example the number of segments or sphere per segments.

Can edit "calc_RHS.m" to change the chosen travelling wave of curvature, or change the reference state of the filament.
