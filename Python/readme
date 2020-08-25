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
