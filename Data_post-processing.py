"""
Here all the functions from the jupyter file "Step_by_step_analysis.ipynb"
are gathered to run the code independently.
"""


import time
import numpy as np
import matplotlib.pyplot as  plt
import math
from scipy.optimize import curve_fit




#_____the global variables needed:________________
window_upper_bound=19e6
window_lower_bound=0e3 #these were taken form the vacuum & EN FFT comparison
n_reduction=1000 #this is necesary for a good curve fitting

names=['0_vac800.txt','4_sqva800badverzmuch.txt', '0_cohe800hz.txt']



#_________All the required functions here:______________
def read( what):
    a = open( '%s' %what, 'r')
    time=[]
    meas=[]
    trigger=[]
    error=[]
    for line in a:
        lines = [i for i in line.split(",")]
        time.append(float(lines[0]))
        meas.append(float(lines[1]))
        trigger.append(float(lines[2]))
        error.append(float(lines[3]))
    a.close()
    minimoa=np.argmin(trigger)
    return time[:minimoa], meas[:minimoa], error[:minimoa]

# This function normalizes the function variance to 0.5 as set by our convention

def norm_variance(signal, std):
    normalized_signal = (signal) /(np.sqrt(2)*std)
    return normalized_signal

#Orain funtxion berria definitukot coherent statetik hasita fase balioak ateatzeko :)
def my_sin(x, freq, amplitude, phase, offset):
    emaitza = []
    for t in x:
        emaitza.append(np.sin(t * float(freq) + float(phase)) * float(amplitude) + float(offset))
    return emaitza

def norm(x):
    X_norm=[]
    X_norm = (x) / (np.max(x) - np.min(x))
    return X_norm

def save(n, data1, data2, what):
    a = open('READY_'+ '_' + '%s' % what, 'w')
    for i in range(len(data1)):
        a.write("%e,%e " % (data1[i], data2[i]) + "\n")
    a.close()
    print('------------> GORDETA ')

def reduction(signal):
    n=0
    copy=[]
    while n<(len(signal)):
        copy.append(signal[n])
        n=n+n_reduction
    return copy

def phase(f, offset, time):
    phase=[]
    phase = [float(f) * x + float(offset) for x in time]
    return phase

def twopi_adjust(phase):
    print('2pi bilatzen')
    save=0
    for x in range(len(phase)-1):
        if phase[x]==phase[-1]:
            pass
            #print('nada', x)
        elif abs(abs(phase[-1]-phase[x])-2*math.pi)<=1e-5 :
            save=x
            print(x)
            #print('puntu honetan dugu ziklo osoa:', x)
            #print(abs(phase[-1] - phase[x])-2*math.pi)
    return phase[save:]

def phase_estimation(signal, time):
    signal_check = norm(signal)
    tope = math.ceil(1.5 * len(signal) / 10)
    t = time[tope:]
    c = norm(signal[tope:])
    t = reduction(t)
    c = reduction(c)
    fit = curve_fit(my_sin, t, c)

    # recreate the fitted curve using the optimized parameters
    data_fit = my_sin(time[tope:], fit[0][0], fit[0][1], fit[0][2], fit[0][3])
    plt.plot(signal_check[tope:], '.')
    plt.plot(data_fit, label='after fitting')
    plt.legend()
    plt.show()
    # calculate the phase and reduce the arrays to a phase change of 2pi only
    theta = phase(fit[0][0], fit[0][2], time[tope:])
    theta = twopi_adjust(theta)
    #  print('2pi check: ', theta[-1]-theta[0])
    new_length = len(time) - len(theta)
    signal = signal[new_length:]

    return signal, theta


def fourier_transform(signal, time):
    freqs = []
    power_spectrum = []
    sn_value = np.array(signal)
    fourier = np.fft.fft(sn_value)
    dt = time[1] - time[0]
    print('Time resolution of the measurement: ', dt)
    # fs=1/dt  #this is the frequency sampling used in the data retrieval
    freqs = np.fft.fftfreq(sn_value.size, d=dt)
    print('The BW of the measured data is: ', np.max(freqs))
    # power_spectrum = np.abs(fourier)
    return fourier, freqs


def inverse_fourier_transform(power_spectrum):
    signal = []
    signal = np.fft.ifft(power_spectrum)
    return signal


def f_window(signal, f):
    windowed_signal = []
    for i in range(len(signal)):
        if np.abs(f[i]) < window_upper_bound and np.abs(f[i]) > window_lower_bound:
            windowed_signal.append(signal[i])
        else:
            windowed_signal.append(0)
    return windowed_signal


def clean(signal, t):
    pow, f = fourier_transform(signal, t)
    pow_win = f_window(pow, f)
    signal_t = inverse_fourier_transform(pow_win)
    return signal_t

def n_cohe(cohe):
    maximum=np.argmax(cohe)
    average=np.average(cohe[maximum-500:maximum+500])
    print('The photon number of the measured coherent state is: ', (average**2)/2)
#______________________________________________
# We adquire the four states and reduce them to the same size preserving the last second data


states=[]
errors=[]
for i in names:
    time=[]
    signal=[]
    error=[]
    time,signal, error=read(i)
    print(len(time))
    time=time[len(time)-250000:]
    signal = signal[len(signal) - 250000:]
    error=error[len(error) - 250000:]
    states.append(signal)
    errors.append(error)

print(len(time))
print(len(states[0]), len(states[1]), len(states[2]))
error=errors[1]
plt.plot(time, states[1], '.', color='tomato')
plt.plot(time, error, '.', color='blue')
plt.xlabel('Time (s)',  fontsize=11)
plt.ylabel('Output voltage (V)',  fontsize=11)
plt.show()



fig = plt.figure()
gs = fig.add_gridspec(3, 1, hspace=0.1, wspace=0)
(ax1), (ax2), (ax3) = gs.subplots(sharex='col', sharey='row')
fig.suptitle('Retrieved  data')
ax1.plot(time, states[0], label='Vacuum')
ax2.plot(time, states[1], label='Squeezed vacuum')
ax3.plot(time, states[2], label='Displaced squeezed light')
ax3.set_xlabel('Time')
plt.show()
"""
# We do a frequency filtering to only consider a flat section
# The filtering boundaries can be adjusted with the variables set at the beginning of the code.
#They were chosen considering the section where vacuum acts as a white noise.
#The lower bound might need to be adjusted not to filter out the oscillations of the coherent state.

"""
states=[clean(x,time) for x in states]



fig = plt.figure()
gs = fig.add_gridspec(3, 1, hspace=0.1, wspace=0)
(ax1), (ax2), (ax3) = gs.subplots(sharex='col', sharey='row')
fig.suptitle('Cleaned  data')
ax1.plot(time, states[0], label='Vacuum')
ax2.plot(time, states[1], label='Squeezed vacuum')
ax3.plot(time, states[2], label='Displaced squeezed light')
ax3.set_xlabel('Time')
plt.show()


#we calculate the phase using the cleaned coherent state
whatever, theta =phase_estimation( states[2], time)

#disp, theta_2 =phase_estimation( states[3][:len(time)-20000], time[:len(time)-20000])
time=time[len(time)-len(theta):]
error=error[len(error)-len(theta):]
error_max=np.max(error)
error_min=np.min(error)
print('the change in the error signal is: ', error_max-error_min)
print('Studied time window: ', time[-1]-time[0])
#time_2=time[len(time)-len(theta_2):]

states[:3]=[x[len(x)-len(theta):] for x in states[:3]]

#now we re-scale it for the vacuum variance to match the convention
vac_std=np.std(states[0])
states=[norm_variance(x, vac_std) for x in states]
#disp=norm_variance(disp, vac_std)


plt.plot(theta, states[1], '.',color='tomato')
plt.xlabel('$\Theta$',  fontsize=11)
plt.ylabel('$\hat{X}_{\\theta}$', fontsize=11)
plt.show()

print('security check: variance of vacuum', np.var(states[0]))
#print('security check: variance of coherent', np.var(states[1][0:1000]))

print(len(states))
for i in range(len(states)):
    print(i)
    save(n,states[i], theta ,names[i])



fig = plt.figure()
gs = fig.add_gridspec(3, 1, hspace=0.1, wspace=0)
(ax1), (ax2), (ax3) = gs.subplots(sharex='col', sharey='row')
fig.suptitle('Output data')
ax1.plot(time, states[0], label='Vacuum')
ax2.plot(time, states[1], label='Squeezed vacuum')
ax3.plot(time, states[2], label='Displaced squeezed light')
ax3.set_xlabel('Time')
plt.show()


""""
plt.subplot(4, 1, 1)
plt.plot(theta, states[0], label='vacuum')
plt.legend(loc='upper right')
plt.subplot(4, 1, 2)
plt.plot(theta, states[1], label='coherent')
plt.legend(loc='upper right')
plt.subplot(4, 1, 3)
plt.plot(theta, states[2], label='squeezed vacuum')
plt.legend(loc='upper right')
plt.subplot(4, 1, 4)
plt.plot(theta_2, disp, label='disp_sq')
plt.xlabel('$\Theta(t) +\\theta_0$')
plt.legend(loc='upper right')
plt.show()
"""
print('The states are ready to be processed by the MaxLik algorithm! :)')



