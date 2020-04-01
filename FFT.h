/*
 *  fftri.h
 *  
 *  Created by Jan Trulsen on 20.11.08.
 *  Copyright 2008 __UiO__. All rights reserved.
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
	
class FFT {
// Fast Fourier Transform of real valued function
	int Nstore;
	double *s, *c;
	
 public:
	FFT();
	~FFT();
	void forward(double *w, const double* y, const int N);
	void inverse(double *y, const double* w, const int N);
	void spectrum(double *p, const double* w, const int N);
};
	
using namespace std;
