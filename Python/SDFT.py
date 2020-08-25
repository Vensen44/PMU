import numpy as np
np.set_printoptions(threshold=np.inf)

def Smart_DFT(fs = 7500, vi=[0]*7500, ndata=750, base_frequency=60, n=125, shift=0):
    '''
    Smart DFT is the precise version of DFT in the case of frequency's variation can't be ignored,
    it use numerical solution to calculate the frequency , phase angel , and magnitude of the signal(voltage)
    that you give to.
    This function can use the vi array to calculate to 10 sets of frequency, angle, and magnitude of the signal(voltage)
    , base on the base_frequency.

    fs : sampling rate, fs = base_frequency*n, fs must be the multiples of base_frequency, or the frequency,
        phase angle ,and magnitude of the signal(voltage) you want to calculate will be wrong.
    vi : the samples of voltage.
    ndata : the number of data which every time we take from array vi in "for count" loop,
        ndata = fs/data_output_per_second.
    base_frequency : frequency base of country, area or you want it to be.
    n : the number of voltage data that is used to calculate to a DFT value.
    shift : shift points of vi array.
    '''


    m = ndata - n - 2
    v = np.array(vi[shift:ndata+shift])
    _Xn = 0
    _Xna = [None] * (ndata - N+1)
    for k in range (0,n):
        _Xn = _Xn + v[k] * np.exp(complex(0,-1) * 2 * np.pi / n * k)  # Discrete Fourier Transform(DFT)
    _Xna[0] = _Xn / n * 2  # every array which doing one time DFT will change to a complex number

    for r in range (2,(ndata - n+2)):
        # create many new Xna(DFT value)
        _Xn = (_Xn - v[r - 2]) * np.exp(complex(0, 1) * 2 * np.pi / n) + v[r + n - 2] * np.exp(
            complex(0, -1) * 2 * np.pi / n * (n - 1))
        _Xna[r - 1] = _Xn / n * 2
        if r > m + 2:
            # do only one time in "for r" loop
            # r - M = 3
            ##### if you want to use the Harmonic_analysis
            fr, _Ar = harmonic_analysis(_Xna, r, m, fs)
            ##### if you don't want to use the Harmonic_analysis
            #_X = (_Xna[0]+_Xna[2])/(2*_Xna[1])
            #o = np.arccos(_X.real)
            #a = np.exp(complex(0, 1) * o)
            #_Ar1 = (_Xna[1]*a - _Xna[0])/(a*a-1)
            #fr1 = o*fs/(2*np.pi)
            #pr1_abs = abs(_Ar1)*n*np.sin(np.pi*(fr1 - base_frequency)/fs)/np.sin(
            #    np.pi * (fr1 - base_frequency) / base_frequency) / np.sqrt(2)  # RMS value of voltage
            #pr1_theta = np.angle(_Ar1) - np.pi / (base_frequency * n) *\
            #    ((fr1 - base_frequency) * (n - 1))  # angle of voltage in radians
            #print (fr1, pr1_abs, pr1_theta)

    fr_abs = np.abs(fr)  # absolute value
    pr_abs = np.abs(_Ar) * n * np.sin(np.pi * (fr_abs - base_frequency) / (base_frequency * n)) / np.sin(
        np.pi * (fr_abs - base_frequency) / base_frequency) / np.sqrt(2)  # RMS value of voltage
    pr_theta = np.angle(_Ar) - np.pi / (base_frequency * n) * ((
        fr_abs - base_frequency) * (n - 1))  # angle of voltage in radians
    return fr_abs[0], pr_abs[0], pr_theta[0]


def harmonic_analysis(Xna=[0], r=0, m=0, fs_=7500):
    '''
    harmonic_analysis is being used in the condition that voltage data is received by Analog to Digital Converter,
    so the signal(voltage) data will combine with some noises or harmonics appearing in wire or in power system witch
    are much smaller than the main signal(smaller than 10% of main signal).

    Xna : the DFT values.
    r : count parameter, r = m+3.
    m : count parameter, m = r-3.
    '''

    Xna = np.array(Xna)
    x1 = Xna[r - m - 3]
    x2 = Xna[r - m - 2]
    _A1 = [None] * m  # the array of values of DFT_r+1 (X_hat_r+1)
    _B1 = [None] * m  # the array of values of DFT_r + DFT_r+2 (X_hat_r + X_hat_r+2)
    for k in range(0, m):
        _A1[k] = 2 * Xna[r-m+k-1]  # k=0 => r-M+k-1 = 2
        _B1[k] = Xna[r-m+k-2] + Xna[r-m+k]  # k=0 => r-M+k-2 = 1, r-M+k = 3

    _A1 = np.array([_A1])  # plus [] to let A1 be a column vector
    _B1 = np.array(_B1)
    _X = np.dot(_A1.conj(), _B1) / np.dot(_A1, _A1.conj().transpose())  # minimum variance of error
    # X will be "a" number ,not a matrix or array
    # A1.conj().transpose() is conjugate transpose

    if _X.real >= 1:
        # in python 1 may = 1.000000000000001 that arccos can't deal with it
        o = 0
    elif _X.real <= -1:
        # in python -1 may = -1.000000000000001 that arccos can't deal with it
        o = np.pi
    else:
        o = np.arccos(_X.real)

    a = np.exp(complex(0, 1) * o)  # if o = 0 => a =1

    return (o*fs_)/(2*np.pi), (x2*a-x1)/(a*a-1)  # frequency and Ar


def make_signal(fs_=7500, freq_=60, angle=0):
    '''
    This function is being used to make a wave.
    You can also add some harmonics or noises at return.

    fs_: sampling rate
    freq_: the frequency of the wave you create
    angle: the angle of the wave you create in degree
    '''

    theta = angle * np.pi/180  # radian
    t = np.arange(0, 1, 1.0 / fs_)  # generating 0 to 1 second array per 1/fs_ second, [0, 1/fs_, 2/fs_, .... , 1]
    return np.cos(2 * np.pi * freq_ * t + theta)
    # signal(voltage) array, you can add some harmonics or noises on the return signal


################################ main code #####################################
fs = 7500  # sampling rate, must > base_frequency * 2
theta = 30  # the angle of signal(voltage) you want to create(Degree)
frequency_siganl = 59  # the frequency of the signal(voltage) you create (Hz)
vi = make_signal(fs, frequency_siganl, theta)  # making a wave using the parameter you give
data_per_set = int(fs * 0.1)  # 10% of signal(voltage) data, must > fs/base_frequency
base_frequency = 60  # frequency base of country , area or you want ot it be (Hz)
N = int(fs/base_frequency)  # N is the number of signal(voltage) data that is used to calculate a DFT value
shift = 0  # shift points of vi array
freq, Vrms, theta = Smart_DFT(fs, vi, data_per_set, base_frequency, N, shift)

print("frequency = ", freq[0], 'Hz')
print("magnitude of signal(voltage) in RMS value = ", Vrms[0], 'Volt')
print("angle of signal(voltage) in radian = ", theta[0], 'rad')

'''
You can find out that if you add a harmonic to make_signal() that is make by the frequency of base_frequency*N , for
N = 1,2,3,... , and the sampling rate is multiple of base_frequency(complying with Nyquist frequency) , then the 
harmonic signal will not affect the result of frequency calculated by Smart_DFT().
ex: in make_signal():
return np.cos(2*np.pi*freq_*t+theta) + np.cos(2*np.pi*(freq_*3)*t+theta)
that means vi = np.cos(2*np.pi*freq_*t+theta) + np.cos(2*np.pi*(freq_*3)*t+theta)
the result of frequency calculated by Smart_DFT() and harmonic_analysis() will be freq_ a.k.a frequency_siganl .

But this method can't calculate the frequency_siganl that is the multiple of base_frequency.
(frequency_siganl = base_frequency * N , N = 2,3,4...)
ex:
frequency_siganl = 100 (Hz)
base_frequency = 50 (Hz)
the result will be wrong.

And if fs(the sampling rate) is not the multiple of base_frequency(complying with Nyquist frequency) , the result of 
angle and magnitude will be wrong.
ex:
base_frequency = 60 (Hz)
fs = 2000 (Hz)
The result fo frequency will be correct.
But the result of angle and magnitude will be wrong.

Let me set the data_per_set to be 10% of the vi array(data_per_set = 0.1 * fs = ndata). In Smart_DFT(),
vi[shift:ndata+shift](shift = 0) gets samples from the first variable in the vi array, so if I set the ndata = 0.1 * fs,
the angle calculated by Smart_DFT() function will be the angle that I set in make_signal().
But if I use vi[i:ndata+i](i != 0 & 750), then the angle calculated by Smart_DFT() will not be the angle I set in
make_signal(), because there is a shift(i points) of the first sample point.
If I set frequency_siganl to be a number of multiples of 10 (ex : 10, 20, 30, ...) and
vi[ndata:ndata+ndata](ndata = 750 = fs * 0.1), then I can get the angle result as the original angle we set, because
the first sample point is the starting point of the waves we generate in make_signal(). And if frequency_siganl is not a
number of multiples of 10 (ex : 68, 21, 103, ...), and the samples is same as
vi[ndata:ndata+ndata](ndata = 750 = fs * 0.1), then the angle will not be the angle I set in make_signal().
'''