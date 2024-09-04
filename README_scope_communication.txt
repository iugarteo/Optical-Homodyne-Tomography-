Requirements: 

- USB connection to the oscilloscope
- Keysight Connection expert
- module pyvisa
- module numpy
- module matplotlib

The initial code was created using he Command Expert program of Keysight. 

Working principle: 
- the function turn_on() displays the channels that we want to capture and adjust the scale, offset, time scale accordingly. 

- the data is adquired by means of data_adq(). 
 The channels that want to be adquired need to be digitalized first. Otherwise the máximum points mode cannot be used and information Will be lost.  The function pream() is used to get the scaling factors of the data. 

- Once the data is adquired it needs to be reescaled by means of data_adjust()

- The retrieved information is saved into a .txt file with the function write()

--> The program adquires three channels containing:
    1) output signal of the 775nm cavity to evaluate the stabilization
    2) trigger signal
    3) BHD output
It would be helpful to implement some conditions here so that the data is not adquired if the stabilization is very bad or the BHD is unbalanced. 


PARAMETERS THAT NEED TO BE ADJUSTED: 
	time_range
	channel3_mV_range
	chan3_offset
	chan2_offset
	chan2_V_range
	chan1_mV_range
	chan1_offset
	time_offset


