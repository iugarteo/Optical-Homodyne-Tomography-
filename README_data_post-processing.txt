Requirements: 
	- module numpy
	- module matplotlib
	- module scipy

Input:
	names=[vac, state, cohe] 

The code does not reconstruct the pase of the 'state', but instead retrieves the pase of the coherent state and selects a length that corresponds to a 2pi pase shift. 


Variables: 
	window_upper_bound
	window_lower_bound

both functions are used in the function clean() to apply a window function in the FFT domain. 

	n_reduction

The number of simples is reduced for the curve fitting of the coherent state. It can be changed accordingly to the size of the sample.

	 
_______________________________________________________
FUNCTIONS:

	read() : reads the raw data and choses the 'linear' section
	save(): saves the data into a .txt file

	
	phase_estimation(): retrieves phase value from coherent light measurement
	--> subfunctions: 
		my_sin(): sin function for curve fit
		norm(): normalization function for curve fit
		reduction(): reduces the size of the sampled data to easen the curve fit
		phase(): reconstructs the pase
		two_pi_adjust(): calculates where do we get a 2pi pase shift. The accuracy may need to be adjusted
        
	
	clean(): cleans the signal by applying a window function in f domain
	--> subfunctions: 
		fourier_transform(): performs a FFT
		f_window(): window function with chosen boundaries
		inverse_fourier_transform(): performs a iFFT


	norm_variance(): normalizes the variance to the chosen convention

	n_cohe(): calculates the average photon number of a coherent state
