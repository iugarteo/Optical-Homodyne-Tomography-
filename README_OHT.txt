Requirements: 
	- numpy module
	- scipy module
	- matplotlib module

Variables: 
	- n_cutoff : the cutoff number for the MaxLik algorithm
	- n_iterations: number of iterations to maximize likelihood

Input: 
	 name=[] #rescaled data for reconstruction

Attention!
--> there are two fock_wave() functions.
	use the simple one for reconstructivos were n_cutoff>=30.
	use the iterative versión for large photon numbers. 