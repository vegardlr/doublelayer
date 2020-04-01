/*
 *  Utils.h
 *  picMC
 *
 *  Created by Jan Trulsen on 22.04.08.
 *  Copyright 2008 __UiO__. All rights reserved.
 *
 */

#ifndef UTILS_H
#define UTILS_H

#define POW2(x) ((x)*(x))

#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>

typedef double real;                            // define datatype real

void ErrMsg(const char* msg, const char* msgdef);


class Histogram {
  // administrate histogram in one, two or three dimensions
	int dim;
	long counts, nbins, nx, ny, nz;
	double lxm, lym, lzm, bx, by, bz;
	int* histogram;
 public:
	Histogram(const double xl, const double xu, const int nix);
	Histogram(const double xl, const double xu, const int nix,
						const double yl, const double yu, const int niy);
	Histogram(const double xl, const double xu, const int nix,
						const double yl, const double yu, const int niy,
						const double zl, const double zu, const int niz);
	~Histogram();
	void reset();
	void record(const double x);
	void record(const double *data, const long start, 
							const long stride, const long N);
	void record(const double *data, const long startx, const long starty,
							const long stride, const long N);
	void record(const double *data, const long startx, const long starty, 
							const long startz, const long stride, const long N);
	void print(const char text[]);
};


class Random
{
  static const int BUCKETSIZE = 1000;
  double PRIMEROOT, last, *bucket;
 public:
	int iBoxMueller;

 public:
  Random(double seed);
  ~Random() { delete [] bucket;};

  double Uniform();
	double BoxMueller();
};


class Inversion
{
  // Draw random numbers (1D or 2D) using numeric direct inversion
	
	private:
		int dim;
		long Lx, Ly, Lz, *lyi, *uyi, *lzi, *uzi;
		double lx, ux, dx, ly, uy, dy, lz, uz, dz, 
					 *Fx, *Fy, *Fz, *lby, *uby, *lbz, *ubz;

  public:
		Inversion(const double xl, const double xu, const int Kx, 
							double (*fx)(double), double& fnorm);
		Inversion(const double xl, const double xu, const int Kx, 
							const double yl, const double yu, const int Ky,
							double (*lfy)(double), double (*ufy)(double), 
							double (*fxy)(double,double), double& fnorm);
		Inversion(const double xl, const double xu, const int Kx, 
							const double yl, const double yu, const int Ky,
							const double zl, const double zu, const int Kz,
							double (*lfy)(double), double (*ufy)(double), 
							double (*lfz)(double,double), double (*ufz)(double,double), 
							double (*fxyz)(double,double,double), double& fnorm);
		~Inversion();
										 
		void draw(double& x);
		void draw(const long N, double* x, const bool quietx = 1);
		void draw(double& x, double& y);
		void draw(const long N, double* x, double* y, const bool quietx = 1);
		void draw(double& x, double& y, double& z);
		void draw(const long N, double* x, double* y, double* z, const bool quietx = 1);

//	private:
		double clint(const int k0, const int k1, const long n0, const long n1,
						const double lb, const double ub, const double lv, const double dv, 
						const double f0, const double f1, const double* f, double* F);
		double interpolate(const double u, const double* F, const long ns,
						const double lb, const double ub, const double lv, const double dv);
		double interpolate(const double u);
		double interpolate(const double u, const int k);
		double interpolate(const double u, const double x);
		double interpolate(const double u, const double x, const double y);
};

double LagrangePol(const double* x, const double* f, const int n, const double X);
double NRLagrangePol(const double* x, const double* f, const int n, const double eps,
									const double u, const double X0);


#endif


