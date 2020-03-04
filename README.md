# HankelSingularValues    
Code implements an efficient algorithm for computing precisely the Hankel singular values of linear singularly perturbed systems. The algorithm is obtained in terms of reduced-order equations and avoids numerical ill-conditioning typically associated with singularly perturbed systems when the singular perturbation parameter is very small.

Implementation is based on the following.
* Kodra, Skataric, Gajic "Finding Hankel singular values for singularly perturbed linear continuous-time systems," IET Control Theory, 2017.

If you use this code, please cite as follows.

	@article{7898625,
      		title   = {{Finding Hankel singular values for singularly perturbed linear continuous-time systems}},
      		author  = {K. {Kodra} and M. {Skataric} and Z. {Gajic}},
      		journal = {IET Control Theory Applications},
      		year    = {2017},
	                volume  = {11},
        	        number  = {7}
}

## 1. Instructions

Running _main.m_ generates the results presented in the paper. See reference in _main.m_ for corrsponding results in paper.

## 2. Dependencies

MATLAB's Control System Toolbox is required to run the code.

