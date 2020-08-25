import SDFT  # import your C file creation
import numpy as np

'''
Set sampling rate(fs_), frequency base(base_frequency), angle of signal(theta), frequency of signal(freq_) and you can begin.
'''
ndata_ = 750 # how many data you use to calculate to a frquecny , angle , voltage
fs_ = 7500  # sample rate, must > base_frequency * 2
theta = 30 * np.pi / 180  # the angleof signal(voltage) you wnat to create(Degree to Radian)
base_frequency = 60  # frequency base of country , area or you want to it be (Hz)
freq_ = 59  # the frequency of signal(voltage) you want to creat (Hz)
N_ = fs_/base_frequency  # N is the number of signal(voltage) data that is used to calculate a DFT value
t = np.arange(0,1,1.0/fs_)  # generating 0 to 1 second array per 1/fs_ second, [0, 1/fs, 2/fs, ...., 1]
vi = np.cos(2*np.pi*freq_*t + theta)[0:ndata_].astype('complex128')
# signal(voltage) array, you can add some harmonics or noises behinde vi array (type is complex)
shift = 0  # shift points of vi array

freq, Vrms, Vangle = SDFT.Smart_DFT(ndata_, fs_, N_, vi, shift, base_frequency)
# call the function, Smart_DFT of the file, SDFT that cython created.
print 'frequency =', freq, 'Hz'
print 'magnitude of signal(voltage) in RMS value = ', Vrms, 'Volt'
print 'angle of signal(voltage) in radian = ', Vangle, 'rad'
