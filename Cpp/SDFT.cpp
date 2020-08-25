#include "pch.h"
#include <iostream>
#include <complex>

#define M_PI 3.14159265358979323846

using namespace std;

/**
 * conjugate_ can let the array of complex number to be the conjugate of itself.
 * ex. [1+2j, 8-5j] to [1-2j, 8+5j]
 */
void conjugate_(complex<double> *matrix_)
{	// conjugate of matrix_
	int rows = sizeof(matrix_);  // rows

	for (int i = 0; i < rows; i++)
	{
		matrix_[i] = conj(matrix_[i]);
	}
}


/**
 * harmonic_analysis is being used in the condition that voltage data is received by Analog to Digital Converter, so the signal(voltage) data will combine with some noises
 * or harmonics appearing in wire or in power system witch are much smaller than the main signal(smaller than 10% of main signal).
 *
 * Xna : the DFT values.
 * r : count parameter, r = m+3.
 * m : count parameter, m = r-3.
 * fs_ : sampling rate.
 * fr : calculated frequency.
 * Ar : calculated voltage parameter, can be used to derive magnitude and angle.
 */
void harmonic_analysis(complex<double> *Xna, const int r, int m, int fs_, complex<double> &fr, complex<double> &Ar)
{
	complex<double> x1 = Xna[r - m - 3];
	complex<double> x2 = Xna[r - m - 2];
	complex<double> *_A1 = new complex<double>[r];  // the array of values of DFT_r + 1 (X_hat_r + 1)
	complex<double> *_B1 = new complex<double>[r]; // the array of values of DFT_r + DFT_r + 2 (X_hat_r + X_hat_r + 2)
	int k; // loop parameter
	int t;  // loop parameter
	complex<double> *_A1_conj = new complex<double>[r];
	complex<double> X;  // for frequency, X = X_up / X_down
	complex<double> X_up = 0;  // inner product parameter
	complex<double> X_down = 0;  // inner product parameter
	double o;  // acos of the real part of X
	complex<double> a;  // for frequency
	double X_real;  // real part of X

	for (int k = 0; k < m; k++)
	{
		_A1[k] = complex<double>(2, 0) * Xna[r - m + k - 1];
		_B1[k] = Xna[r - m + k - 2] + Xna[r - m + k];
	}

	for (int a = 0; a < r; a++)  // copy
		_A1_conj[a] = _A1[a];
	
	conjugate_(_A1);


	for (int t = 0; t < m; t++)
	{	//do "product", much faster than call the python function
		X_up = X_up + _B1[t] * _A1[t];
		X_down = X_down + _A1[t] * _A1_conj[t];
	}
	X = X_up / X_down;
	X_real = X.real();

	if (X_real >= 1)  // in python 1 may = 1.000000000000001 that arccos can't deal with it
		o = 0;
	else if (X_real <= -1)// in python - 1 may = -1.000000000000001 that arccos can't deal with it
		o = M_PI;
	else
		o = acos(X_real);

	a = exp(complex<double>(0, 1 * o));
	fr = (o*fs_) / (2 * M_PI);  // frequency
	Ar = (x2*a - x1) / (a*a - complex < double>(1, 0));  // Ar
}


/**
 * Smart DFT is the precise version of DFT in the case of frequency's variation can't be ignored, it use numerical solution to calculate the frequency , phase angel ,
 * and magnitude of the signal(voltage) that you give to.
 * This function can use the vi array to calculate to 10 sets of frequency, angle, and magnitude of the signal(voltage), base on the base_frequency.
 *
 * fs : sampling rate, fs = base_frequency*n, fs must be the multiples of base_frequency, or the frequency,
 *	phase angle ,and magnitude of the signal(voltage) you want to calculate will be wrong.
 * vi : the samples of signal(voltage).
 * ndata : the number of data which every time we take from array vi in "for count" loop, ndata = fs/data_output_per_second.
 * N : the number of voltage data that is used to calculate to a DFT value.
 * base_frequency : frequency base of country, area or you want it to be.
 * shift : shift points of vi array. 
 * fr_abs : calculated frequency in Hz.
 * pr_abs : calculated magnitude of signal(voltage) in RMS value.
 * pr_theta : calculated angle of signal(voltage) in radian.
 */
void Smart_DFT(float fs, complex<double> vi[7500], int ndata, int N, float base_frequency, int shift, double &fr_abs, double &pr_abs, double &pr_theta)
{
	
	complex<double> *v = new complex<double>[ndata];
	complex<double> fr;  // frequency
	complex<double> Ar;  // can calculate the phase angle and magnitude of voltage
	int M = ndata - N - 2;  // condition of if / else and loop limit value
	complex<double> Xn = 0;  // DFT vaule
	const int Xna_size = ndata - N + 1;
	complex<double> *Xna = new complex<double>[Xna_size];  // DFT array, storaging DFT data
	complex<double> power_index;

	for (int h = shift; h < ndata+shift; h++)
		v[h - shift] = vi[h];

	for (int k = 0; k < N; k++)
	{
		power_index = complex<double>(0, -1*2*M_PI/N*k);
		Xn = Xn + v[k] * exp(power_index); // Discrete Fourier Transform(DFT)
	}
	Xna[0] = Xn / complex<double>(N, 0) * complex<double>(2, 0); // every array which doing one time DFT will change to a number

	for (int r = 2; r < (ndata - N + 2); r++)
	{	// create many new Xna(DFT value)
		Xn = (Xn - v[r - 2]) * exp(complex<double>(0, 1 * 2 * M_PI / N)) + v[r + N - 2] * exp(complex<double>(0, -1 * 2 * M_PI / N * (N - 1)));
		Xna[r - 1] = Xn / complex<double>(N, 0) * complex<double>(2, 0);

		if (r > M + 2)
		{	// do only one time in "for r" loop
			// r - M = 3
			harmonic_analysis(Xna, r, M, fs, fr, Ar);
		}
		
	}
	fr_abs = abs(fr); // absolute value
	pr_abs = abs(Ar) * N * sin(M_PI*(fr_abs - base_frequency) / (base_frequency * N)) / sin(
		M_PI * (fr_abs - base_frequency) / base_frequency) / sqrt(2);  // RMS value of voltage
	pr_theta = arg(Ar) - M_PI / (base_frequency * N) * ((
		fr_abs - base_frequency) * (N - 1)); // angle of voltage in radians
}



/* Set sampling rate(fs_), frequency base(base_frequency), angle of signal(theta), frequency of signal(freq_) and you can begin*/
int main()
{
	double freq, Vrms, Vangle;
	float fs_ = 7500;  // sampling rate, must > base_frequency * 2
	float base_frequency = 60;  // frequency base of country, area or you want to it be(Hz)
	int ndata_ = fs_ * 0.1; // how many data you use to calculate to a frquecny, angle, and voltage(must > fs_ / base_frequency)
	double theta = 30 * (M_PI / 180);  // the angle of signal(voltage) you want to creat(degree to radian)
	float freq_ = 59;  // the frequency of signal(voltage) you want to creat(Hz)
	int N_ = ceil(fs_ / base_frequency);  // N is the number of signal(voltage) data that is used to calculate a DFT value
	float *t = new float [fs_];  // generating 0 to 1 second array per 1 / fs_ second, [0, 1 / fs, 2 / fs, ...., 1]
	complex<double> *vi = new complex<double> [fs_];
	int shift_ = 0;  // shift points of vi array

	for (int i = 0; i< fs_; i++)
	{	// signal(voltage) array, you can add some harmonics or noises behinde vi array(type is complex)
		vi[i] = complex<double> (cos(2 * M_PI * freq_ * (i / fs_) + theta), 0);
	}

	Smart_DFT(fs_, vi, ndata_, N_, base_frequency, shift_, freq, Vrms, Vangle);
	
	cout << "Frequency(Hz) = " << freq  << endl;
	cout << "Voltage(RMS) = " << Vrms  << endl;
	cout << "Angle(Radian) = " << Vangle  << endl;
}
