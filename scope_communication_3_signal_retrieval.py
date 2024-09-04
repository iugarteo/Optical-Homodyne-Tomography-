# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:53:57 2024

Imported functions from the Command Expert program of Python


@author: IDOIA
"""
# NOTE: the default pyvisa import works well for Python 3.6+
# if you are working with python version lower than 3.6, use 'import visa' instead of import pyvisa as visa


import pyvisa as visa
import time
import numpy as np
import matplotlib.pyplot as  plt

rm = visa.ResourceManager()
DSO_X_3034G = rm.open_resource('USB0::0x2A8D::0x1764::MY62510103::0::INSTR')
DSO_X_3034G.write('*CLS')

#Values that have to be changed and optimized for data adquisition

time_range=1.0e-3
channel3_mV_range=60#
chan3_offset=0.0
chan2_offset=2.0
chan2_V_range=20.0
chan1_mV_range=20.0
chan1_offset=1
time_offset=-600e-6

#Choose the appropiate labelling of the measurements
what='AA'
#what='vacuum_trigger_test_ramp30'
what='3.02-squeezing-400hz'
what='test'
#what='vac_200hz'
#Number of a1dquisitions to be made:
n_meas=3
#_______________________________________________________________________________________________________________

def turn_on():
    DSO_X_3034G.write(':CHANnel2:DISPlay %d' % (1))
    DSO_X_3034G.write(':CHANnel3:DISPlay %d' % (1))
    DSO_X_3034G.write(':CHANnel1:OFFSet %G' % (chan1_offset))
    DSO_X_3034G.write(':CHANnel2:OFFSet %G' % (chan2_offset))
    DSO_X_3034G.write(':CHANnel3:OFFSet %G' % (chan3_offset))
    DSO_X_3034G.write(':CHANnel2:RANGe %G' % (chan2_V_range))
    DSO_X_3034G.write(':CHANnel1:RANGe %G mV' % (chan1_mV_range))
    DSO_X_3034G.write(':CHANnel3:RANGe %G mV' % (channel3_mV_range))
    DSO_X_3034G.write('*WAI')
    DSO_X_3034G.write(':CHANnel2:IMPedance %s' % ('FIFTy'))
    DSO_X_3034G.write(':CHANnel3:IMPedance %s' % ('FIFTy'))
    DSO_X_3034G.write(':CHANnel1:IMPedance %s' % ('FIFTy'))
    DSO_X_3034G.write(':CHANnel2:PROBe %G' % (1.0))
    DSO_X_3034G.write(':CHANnel3:PROBe %G' % (1.0))
    DSO_X_3034G.write(':CHANnel1:PROBe %G' % (1.0))
    #DSO_X_3034G.write(':TIMebase:SCALe %G' % (1.4e-08))
    DSO_X_3034G.write(':TIMebase:RANGe %G' % (time_range))
    DSO_X_3034G.write(':TIMebase:POSition %G' % (time_offset)) #a time delay of -3ms is added so that we do not getb the edge in every data adquisition
    print("Turn on rutine executed")
 


def pream():
    temp_values = DSO_X_3034G.query_ascii_values(':WAVeform:PREamble?')
    format = int(temp_values[0])
    type = int(temp_values[1])
    points = int(temp_values[2])
    count = int(temp_values[3])
    xincrement = temp_values[4]
    xorigin = temp_values[5]
    xreference = int(temp_values[6])
    yincrement = temp_values[7]
    yorigin = temp_values[8]
    yreference = int(temp_values[9])
    values=[xincrement,xorigin,xreference, yincrement,yorigin,yreference]
    return values

def data_adq():
    #I am going to analyze the values we gha 
    #it is important to digituze to have the mXIMUM AMOUNT OF DATA POINTS
    DSO_X_3034G.write(':DIGitize %s,%s' % ('CHAN2', 'CHAN3', 'CHAN1'))
    DSO_X_3034G.write(':WAVeform:POINts:MODE %s' % ('MAXimum'))
    pointsMode = DSO_X_3034G.query(':WAVeform:POINts:MODE?')
    DSO_X_3034G.write(':WAVeform:SOURce %s' % ('CHANnel2'))
    DSO_X_3034G.write(':WAVeform:FORMat %s' % ('BYTE'))
    value2=pream()
    data2 = DSO_X_3034G.query_binary_values(':WAVeform:DATA?','B',False)
    
    DSO_X_3034G.write(':DIGitize %s,%s' % ('CHANnel3', 'CHAN2', 'CHAN1'))
    DSO_X_3034G.write(':WAVeform:SOURce %s' % ('CHANnel3'))
    DSO_X_3034G.write(':WAVeform:FORMat %s' % ('BYTE'))
    value1=pream()
    data1 = DSO_X_3034G.query_binary_values(':WAVeform:DATA?', 'B', False)

    DSO_X_3034G.write(':DIGitize %s,%s' % ('CHANnel3', 'CHAN2', 'CHAN1'))
    DSO_X_3034G.write(':WAVeform:POINts:MODE %s' % ('MAXimum'))
    DSO_X_3034G.write(':WAVeform:SOURce %s' % ('CHANnel1'))
    DSO_X_3034G.write(':WAVeform:FORMat %s' % ('BYTE'))
    value3 = pream()
    data3 = DSO_X_3034G.query_binary_values(':WAVeform:DATA?', 'B', False)


    

    
    return [data1,value1,data2,value2, data3, value3]



#I should implement a condition here for the saved values: taking only the relevant ones (not messed with the piezo)

def write(data1,data2,data3,data4,n):
    a = open('%i' %n+'_'+'%s' %what+'.txt', 'w')
    for i in range(len(data1)):
        a.write("%e,%e,%e, %e "% (data1[i], data2[i], data3[i], data4[i])+"\n")
    a.close()
    
def data_adjust(data,value):
    tme=[]
    wfm=[]
    for t in range (len(data)):
        tme.append((t*value[0])+value[1])
    for d in data:
       # wfm.append((d*value[1])+value[0])
        wfm.append(((d-value[5])*value[3])+value[4])
      #  voltage = ((values[i] - y_reference) * y_increment) + y_origin
    return [tme,wfm]
        

#______________________________________________________________________________________________________________   
    

for n in range(n_meas):
    print(n)
    turn_on()
    signal, s_values, trigger,t_values, error, t_error=data_adq()
    
    
    
    time,signal=data_adjust(signal, s_values)
    time_stab, stab_error=data_adjust(error, t_error)
        
    write(time,signal,trigger,stab_error, n) #To save the data in a txt file
    
    temp_values = DSO_X_3034G.query_ascii_values('*OPC?')
    complete = int(temp_values[0])
    
    if complete==0:
        DSO_X_3034G.write('*WAI')
        


#plot waveform data
    #plt.plot(time,signal)
   # plt.plot(time,trigger)
    #plt.title('Captured data')
    #plt.xlabel('Time (sec)')
    #plt.ylabel('Voltage (V)')
    #plt.show()
    #plt.savefig('%s.png' %n_meas)


print('-----------------------------------------------------------------------------------')
print('DONE')
DSO_X_3034G.write(':RUN')
DSO_X_3034G.close()
rm.close()

