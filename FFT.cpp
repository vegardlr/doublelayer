/*
 *  fftri.cpp
 *  
 *  Created by Jan Trulsen on 20.11.08.
 *  Copyright 2008 __UiO__. All rights reserved.
 */

#include "FFT.h"



FFT::FFT()
{
	Nstore = 0;
}


FFT::~FFT()
{
  delete [] s;
	delete [] c;
}

void FFT::forward(double *w, const double *y, const int N)
{
	// forward discrete fast Fourier transform of real valued y
	int NN, Nh, nh, n, ni, nii, l, m, ii, i, j, k, groups, levels;
	double a, cn, sn;

	double *x = new double [N];

	Nh = nh = N/2;
	if (N != Nstore) {
		if (Nstore != 0) { delete [] c, delete [] s; }
		c = new double [N], s = new double [N];
		double t1, t2;
		
		t2 = (t1 = tan(M_PI/N))*t1;
		c[0] = 1., s[0] = 0., c[1] = cn = (1.-t2)/(1.+t2), s[1] = sn = 2.*t1/(1.+t2);
		for (i = 2; i < Nh; ++i) { 
			j = i-1, c[i] = c[j]*cn-s[j]*sn, s[i] = c[j]*sn+s[j]*cn;
		}
		Nstore = N;
		//cout<<"Nstore = "<<Nstore<<endl;
	}
	
	levels = 0, j = N;
	while ((j /= 2) && !(j%2)) ++levels;
	NN = 2<<levels;
//	cout<<"NN levels "<<NN<<' '<<levels<<endl;
	if (NN != N) cout<<"WARNING NN != N /n";

	for (i = 0, j = i+Nh; i < Nh; ++i, ++j) {		// fold to N/2 level 1
		w[i] = y[i]+y[j], w[j] = y[i]-y[j];
		x[i] = 0., x[j] = w[j]*s[i], w[j] *= c[i];
	}

	for (l = 2, groups = 1, ni = 1; l <= levels; ++l) {	// keep on folding
		n = nh, nh /= 2, groups *= 2, ni *= 2;
		for (m = 1, k = 0; m <= groups; ++m, k += n) {
			for (ii = 0; ii < nh; ++ii) {
				i = ii+k, j = i+nh, nii = ni*ii;
				a = w[j], w[j] = w[i]-w[j], w[i] += a;
				a = x[j], x[j] = x[i]-x[j], x[i] += a;
				a = x[j], x[j] = w[j]*s[nii]+a*c[nii], w[j] = w[j]*c[nii]-a*s[nii];
			}
		}
	}
	
	for (i = 0; i < N; i += 2) {   // last folding
		j = i+1, x[j] = x[i]+x[j], x[i] = w[i]+w[j];
	}

	// ordering and normalization
	a = 1./N, w[0] *= a; w[Nh] = 0.;
	for (l = 1, n = N, k = 1; l <=levels; ++l, n /= 2)
		for (i = n/2; i < N; i += n) {
			w[k] = a*x[i], w[N-k] = a*x[i+1], ++k;
		}
	
	delete [] x;
}

void FFT::inverse(double *y, const double *w, const int N)
{
	// inverse discrete fast Fourier transform with "w[N-k] = w[k]^*"
	int NN, Nh, nh, n, ni, nii, l, m, ii, i, j, k, groups, levels;
	double a, cn, sn;

	double *x = new double [N];

	Nh = nh = N/2;
	if (N != Nstore) {
		if (Nstore != 0) { delete [] c, delete [] s; }
		c = new double [N], s = new double [N];
		double t1, t2;
		
		t2 = (t1 = tan(M_PI/N))*t1;
		c[0] = 1., s[0] = 0., c[1] = cn = (1.-t2)/(1.+t2), s[1] = sn = 2.*t1/(1.+t2);
		for (i = 2; i < Nh; ++i) { 
			j = i-1, c[i] = c[j]*cn-s[j]*sn, s[i] = c[j]*sn+s[j]*cn;
		}
		Nstore = N;
		cout<<"Nstore = "<<Nstore<<endl;
	}
	
	levels = 0, j = N;
	while ((j /= 2) && !(j%2)) ++levels;
	NN = 2<<levels;
	if (NN != N) cout<<"WARNING NN != N /n";
	
	x[0] = w[0], x[Nh] = y[0] = y[Nh] = 0.;
	for (i = 1; i < Nh; ++i) {  // first folding using "w[N-k] = w[k]^*"
		j = i+Nh;
		x[i] = w[i]+w[Nh-i], y[i] = w[N-i]-w[Nh+i];
		x[j] = w[i]-w[Nh-i], y[j] = w[N-i]+w[Nh+i];
		a = x[j], x[j] = y[j]*s[i]+a*c[i], y[j] = y[j]*c[i]-a*s[i];
	}	
	
	for (l = 2, groups = 1, ni = 1; l <= levels; ++l) {	// keep on folding
		n = nh, nh /= 2, groups *= 2, ni *= 2;
		for (m = 1, k = 0; m <= groups; ++m, k += n) {
			for (ii = 0; ii < nh; ++ii) {
				i = ii+k, j = i+nh, nii = ni*ii;
				a = y[j], y[j] = y[i]-y[j], y[i] += a;
				a = x[j], x[j] = x[i]-x[j], x[i] += a;
				a = x[j], x[j] = y[j]*s[nii]+a*c[nii], y[j] = y[j]*c[nii]-a*s[nii];
			}
		}
	}
	
	for (i = 0; i < N; i += 2) {   // last folding
		j = i+1, a = x[i], x[i] = a+x[j], x[j] = a-x[j];
	}

	// ordering
	int *ia = new int[N/4];

	y[0] = x[0], y[Nh] = x[1], y[1] = x[Nh], y[Nh+1] = x[Nh+1]; 
	k = Nh/2;
	i = 2, y[i] = x[k], y[i+Nh] = x[k+1];
	++i, y[i] = x[k+Nh], y[i+Nh] = x[k+Nh+1];
	for (n = Nh/2, j = 1, ia[0] = 2; n >= 4; n /= 2, j *= 2) {
		k = 2*j;
		for (i = j-1, l = 2*j-1; i >= 0; --i, l -= 2) {
			ia[l] = ia[i]+2*k, ia[l-1] = ia[i]+k;
		}
		for (k = n/2, l = 0; k < Nh; k += n, ++l) {
			i = ia[l], y[i] = x[k], y[i+Nh] = x[k+1];
		++i, y[i] = x[k+Nh], y[i+Nh] = x[k+1+Nh];
		}
	}
	
	delete [] ia;
	delete [] x;
}

void FFT::spectrum(double *p, const double *w, const int N)
{
	// Power spectrum p of FT w: "p = w*conj(w)"
	
	p[0] = w[0]*w[0];
	for (int i = N/2-1; i > 0; --i) p[i] = 2.*(w[i]*w[i]+w[N-i]*w[N-i]);
}

