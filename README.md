# nnparafac2
Decomposition of multi-dimensional arrays with drift in one mode. As described by Cohen and Bro:

Cohen, Jeremy E., and Rasmus Bro. "Nonnegative PARAFAC2: A flexible coupling approach." International Conference on Latent Variable Analysis and Signal Separation. Springer, Cham, 2018.

Features an estimator for the signal-to-noise ratio for each slice, which is used to determine an appropriate value for the coupling constant, mk. This is an original contribution to the algorithm.

Currently a dependency on fcnnls.m as a non-negative least-squares solver: https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e5925d25-4a80-11e4-9553-005056977bd0/0b2678a0-fdea-4929-820c-ef5a41b26b46/previews/source/fcnnls.m/index.html 

Syntax:
[Bk,A,Dk,Bs,ssr,timeOut] = nnparafac2(Xk,R,Bki*,Ai*,Dki*,Bsi*)

Where Xk is a cell array or tensor of K samples, where each sample is an I times J matrix of I observations and J variables, and I is not necessarily the same across K samples iff the input Xk is a cell array. Commonly used to model chromatographic drift in hyphenated chromatographic experiments with multiple samples. R is the number of factors, to be decided by the user. All inputs marked by a * are optional values for initialisation of the model. If none are provided, the initialisation proceeds following several random initialisations. Otherwise, if at least one is provided, all other initialisations are randomised once.
