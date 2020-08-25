import numpy as np
cimport numpy as np
from libc.math cimport acos,fabs,sqrt,M_PI # import math-function from C library

cdef extern from "<complex.h>": 
    # import complex exponentail from c library
    double complex cexp(double complex)

cpdef tuple Smart_DFT(int ndata, int fs, int N, np.ndarray[np.complex128_t, ndim=1] vi, int shift, float base_frequency):
    '''
    cpdef: cython and python file can call it , much faster than def in python.
    
    Smart DFT is the precise version of DFT in the case of frequencyls variation can't be ignored,
    it use numerical solution to calculate the frequency, phase angle, and magnitude of the signale(voltage)
    that you give to.
    This function can use v array to calculate to 1 sets of frequency, angle, and magnitude of the signale(voltage)
    , base on the base_frequency.

    ndata : the number of data you want to calculate to a frequency, angle, voltage.
    fs : sampling rate, fs = base_frequency*N, fs must be the multiples of base_frequency, or the frequency,
        phase angle, and magnitude of the signal(voltage) you wnat to calculate will be wrong.
    N : the number of voltage data that is used to calculate to a DFT value.
    v : the samples of voltage.
    shift : shift points of vi array
    base_freuqncy : frequency base of country, area or you wnat it to be.
    '''
    cdef:  # cdef: only cython can call it, fastest
        np.ndarray[np.complex128_t, ndim=1] v  # shift of vi array
        int g  # loop index
        double complex fr  # frequency
        double complex Ar  # can calculate the phase angle and magnitude of voltage
        int M = ndata - N - 2  # condition of if/else and loop limit value
        double complex Xn = 0  # DFT vaule
        np.ndarray[np.complex128_t, ndim=1] Xna  # DFT array, storaging DFT data
        int k  # loop index
        int r  # loop index
        float fr_abs  # frequency output
        float pr_abs  # Vrms output
        float pr_theta  # theta output (radian)
    v = np.zeros(ndata, dtype = np.complex128)
    Xna = np.zeros((ndata - N + 1) ,dtype = np.complex128)
    
    for h in range(shift,ndata+shift):  # shift of vi array
        v[h-shift] = vi[h]

    for k in range (0,N):
        Xn = Xn + v[k] * cexp(-1j * 2 * M_PI / N * k) # Discrete Fourier Transform(DFT)
    Xna[0] = Xn / N * 2 # every array which doing one time DFT will change to a number

    for r in range (2,(ndata-N+2)):
        # create many new Xna(DFT value)
        Xn = (Xn - v[r-2]) * cexp(1j * 2 * M_PI / N) + v[r+N-2] * cexp(-1j * 2 * M_PI / N * (N-1))
        Xna[r-1] = Xn / N * 2
        if (r > M+2) :
            # do only one time in "for r" loop
            # r - M = 3
            fr, Ar = harmonic_analysis(Xna, r, M, fs)

    fr_abs = np.abs(fr) # absolute value
    pr_abs = np.abs(Ar) * N * np.sin(np.pi*(fr_abs - base_frequency) / (base_frequency * N))/ np.sin(
            np.pi * (fr_abs - base_frequency) / base_frequency) / sqrt(2)  # RMS value of voltage
    pr_theta = np.angle(Ar) - np.pi / (base_frequency * N) * ((
        fr_abs - base_frequency) * (N-1)) # angle of voltage in radians
    return(fr_abs , pr_abs , pr_theta)

cdef tuple harmonic_analysis(np.ndarray[np.complex128_t, ndim=1] Xna, int r, int m, int fs_):
    '''
    cdef : only c file can call it, the fastest.
    
    harmonic_analysis is being used in the condition that voltage data is received by Analog to Digital Converter,
    so the signal(voltage) data will combine with some noise or harmonics appearing in wire or power system witch
    are much smaller than the main signale(smaller than 10% of main siganle).

    Xna : the DFT values.
    r : count parameter, r = m+3.
    m : count parameter, m = r-3.
    '''
    cdef:
        double complex x1 = Xna[r-m-3]
        double complex x2 = Xna[r-m-2]
        np.ndarray[np.complex128_t, ndim=1] _A1  # the array of values of DFT_r+1 (X_hat_r+1)
        np.ndarray[np.complex128_t, ndim=1] _B1  # the array of values of DFT_r + DFT_r+2 (X_hat_r + X_hat_r+2)
        int k  # loop parameter 
        int t  # loop parameter
        double complex X  # for frequency, X = X_up/X_down
        double complex X_up = 0  # inner product parameter
        double complex X_down = 0  # inner product parameter
        double o  # acos of the real part of X
        double complex a  # for frequency
        double X_real  # real part of X 
    _A1 = np.zeros(r, dtype = np.complex128)
    _B1 = np.zeros(r, dtype = np.complex128)

    for k in range(0,m):
        _A1[k] = 2 * Xna[r-m+k-1]
        _B1[k] = Xna[r-m+k-2] + Xna[r-m+k]
    _A1_conj = _A1.conjugate()

    for t in range(0,m):
    # do "product" , much faster than call the python function
        X_up = X_up + _B1[t] * _A1_conj[t]
        X_down = X_down +_A1_conj[t] * _A1[t]
    X = X_up / X_down
    X_real = X.real
    
    if X_real >= 1:
    # in python 1 may = 1.000000000000001 that arccos can't deal with it
        o = 0
    elif X_real <= -1:
    # in python -1 may = -1.000000000000001 that arccos can't deal with it
        o = M_PI
    else:
        o = acos(X_real)
    a = cexp(1j*o)

    return (o*fs_)/(2*np.pi), (x2*a-x1)/(a*a-1) # frequency and Ar
