/*
 *  Utils.cpp
 *  picMC
 *
 *  Created by Jan Trulsen on 22.04.08.
 *  Copyright 2008 __UiO__. All rights reserved.
 *
 */

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>
#include "Utils.h"

using namespace std;

//#define TEST

void ErrMsg(const char* msg, const char* msgdef = "")
{
  cerr<<msg<<' '<<msgdef<<endl;
	exit(1);
}




Histogram::Histogram(const double xl, const double xu, const int nix)
{
  // Inside range: [xl,xu>
	// Bins: inside 1,...,nix,  outside 0, nix+1

  if (nix <= 0 || xu <= xl) cout<<"invalid histogram range, bins : "
																<<xl<<" "<<xu<<" "<<nix<<endl;
  dim = 1, counts = 0L;
	bx = (xu-xl)/nix;
	lxm = xl-bx;
	nx = nix,	nbins = nx;
	
	histogram = new int[nbins];
	for (long i = 0; i < nbins; ++i) histogram[i] = 0;
}


Histogram::Histogram(const double xl, const double xu, const int nix,
										 const double yl, const double yu, const int niy)
{
  // Inside range: [xl,xu> x [yl,yu>
	// Bins: inside (0,...,nix-1) x (0,...,niy-1)

  if (nix <= 0 || xu <= xl || niy <= 0 || yu <= yl) 
		cout<<"invalid histogram range, bins : "
				<<xl<<" "<<xu<<" "<<nix<<" "<<yl<<" "<<yu<<" "<<niy<<endl;
  dim = 2, counts = 0L;
	bx = (xu-xl)/nix, by = (yu-yl)/niy;
	lxm = xl-bx, lym = yl-by;
	nx = nix, ny = niy, nbins = nx*ny;
	
	histogram = new int[nbins];
	for (long i = 0; i < nbins; ++i) histogram[i] = 0;
}


Histogram::Histogram(const double xl, const double xu, const int nix,
										 const double yl, const double yu, const int niy,
										 const double zl, const double zu, const int niz)
{
  // Inside range: [xl,xu> x [yl,yu> x [zl,zu> 
	// Bins: (0,...,nix-1) x (0,...,niy-1) x (0,...,niz-1),  

  if (nix <= 0 || xu <= xl || niy <= 0 || yu <= yl || niz <= 0 || zu <= zu) 
		cout<<"invalid histogram range, bins : "<<xl<<" "<<xu<<" "<<nix<<" "
				        <<yl<<" "<<yu<<" "<<niy<<" "<<zl<<" "<<zu<<" "<<niz<<endl;
																						
  dim = 3, counts = 0L;
	bx = (xu-xl)/nix, by = (yu-yl)/niy, bz = (zu-zl)/niz;
	lxm = xl-bx, lym = yl-by, lzm = zl-bz;
	nx = nix, ny = niy, nz = niz, nbins = nx*ny*nz;
	
	histogram = new int[nbins];
	for (long i = 0; i < nbins; ++i) histogram[i] = 0;
}


Histogram::~Histogram()
{
  delete [] histogram;
}


void Histogram::reset()
{
  counts = 0L;
	for (long i = 0L; i < nbins; ++i) histogram[i] = 0;
}



void Histogram::record(const double x)
{
  // add single record to histogram
  long i;
	
	i = (long)((x-lxm)/bx)-1L;        // -1 to avoid floor !
	if (i >= 0L && i < nx) {
		histogram[i]++;
		counts++;
	} else cout<<i<<" "<<x<<endl;
}


void Histogram::record(const double* data, const long start, 
											 const long stride, const long N)
{
  // add multiple records to histogram
  int i;
	
	for (long n = start+(N-1)*stride; n >= start; n -= stride) {
		i = (long)((data[n]-lxm)/bx)-1L;        // -1 to avoid floor !
		if (i >= 0L && i < nx) {
			histogram[i]++;
			counts++;
		} 
	}
}
 

void Histogram::record(const double* data, const long startx, const long starty,
											 const long stride, const long N)
{
  // add multiple records to histogram
  long i, j;
	
	for (long n = (N-1)*stride; n >= 0; n -= stride) {
		i = (long)((data[n+startx]-lxm)/bx)-1L;
		if (i >= 0L && i < nx) {
			j = (long)((data[n+starty]-lym)/by)-1L;
			if (j >= 0L && j < ny) {
				histogram[i+nx*j]++;
				counts++;
			}
		}
	}
}
 

void Histogram::record(const double* data, const long startx, const long starty,
											 const long startz, const long stride, const long N)
{
  // add multiple records to histogram
  long i, j, k, nxy = nx*ny;
	
	for (long n = (N-1)*stride; n >= 0; n -= stride) {
		i = (long)((data[n+startx]-lxm)/bx)-1L;
		if (i >= 0L && i < nx) {
			j = (long)((data[n+starty]-lym)/by)-1L;
			if (j >= 0L && j < ny) {
				k = (long)((data[n+startz]-lzm)/bz)-1L;
				if (k >= 0L && k < nz) {
					histogram[i+nx*j+nxy*k]++;
					counts++;
				}
			}
		}
	}
}
 

void Histogram::print(const char text[])
{
	cout<<text<<" "<<counts<<endl;
	if (dim == 1) {
		long i;
		for (i = 0; i < nx; ++i) cout<<histogram[i]<<" "; cout<<endl;
	}
	if (dim == 2) {
		long i, j, ii;
		for (j = ny-1; j >= 0; --j) {
			for (i = 0, ii = i+nx*j; i < nx; ++i, ii += 1) 
				cout<<histogram[ii]<<" "; cout<<endl;
		}
	}
}




Random::Random(double seed)
{
  iBoxMueller = 0;
	
  PRIMEROOT = sqrt(2.0);
	bucket = new double[BUCKETSIZE];
  last = seed;
  for (int n = 0; n < BUCKETSIZE; ++n) {
    last += PRIMEROOT; last -= floor(last);
		bucket[n] = last;
  }
}



double Random::Uniform()
{
  // Quiet uniform random numbers in (0.,1.)

  int draw = int(drand48()*BUCKETSIZE);
  double r = bucket[draw];
  last += PRIMEROOT, last -= floor(last); 
  bucket[draw] = last;

  return r;
}


double Random::BoxMueller()
{
  static double u;
	
  if (iBoxMueller != 0) {
    iBoxMueller = 0; 
    return u;
  }
  else {
    double x, y, r2;

    for (;;) {
			x = 2.*Uniform()-1.,  y = 2.*Uniform()-1.;
			r2 = x*x+y*y;
			if (r2 < 1.) {
				r2 = sqrt(-2.*log(r2)/r2);
				iBoxMueller = 1, u = y*r2;
 				return x*r2;
			}
    }
  }
}
      


Inversion::Inversion(const double xl, const double xu, const int Kx, 
										 double (*fx)(double), double& fnorm)
{
/*
 *	Prepare draw of r.v. x in fx(x) by direct numerical inversion
 *	Maximum (compact) x domain: (xl,xu) divided in KK cells
 *	Returns Fnorm = \int_xl^xu fx(x) dx via trapes method
*/

  int i, i0, i1;
	long l0, l1;
	double x, f0, f1, *f;
	
	dim = 1;
	lx = xl, ux = xu, Lx = Kx, dx = (ux-lx)/(Kx-1);

	Fx = new double[Lx];

  f = new double[Lx];
	
	// Sample pdf fx(x)
	i0 = 0, i1 = Kx-1, l0 = i0, l1 = i1;
	f0 = f1 = 0.;
	for (i = i0, x = lx; i <= i1; ++i, x += dx) f[i] = fx(x);
	
	// Evaluate cumulative Fx(x) = \int^x fx(x) dx
	fnorm = clint(i0,i1,l0,l1, lx,ux,lx,dx, f0,f1,f, Fx)*dx;

	delete [] f;
}


Inversion::Inversion(const double xl, const double xu, const int Kx, 
										 const double yl, const double yu, const int Ky,
										 double (*lfy)(double), double (*ufy)(double), 
										 double (*fxy)(double, double), double& fnorm)
{

/*
 *	Prepare draw of r.v. (x,y) in fxy(x,y) by direct numerical inversion
 *	Maximum (compact) xy rectangle: (xl,xu) x (yl,yu) divided in KK x LL cells
 *	Lower and upper boundary for y: lfy(x) and ufy(x)
 *  Method: draw x from fx(x) = \int fxy(x,y) dy in (lx,ux)
 *          and y given x from fxy(y|x) in (lfy(x),ufy(x))
 *	Returns fnorm = \int \int fxy(x,y) dydx
 */

  int i, i0, i1, j, j0, j1;
	long l0, l1, m0, m1;
	double x, y, f0, f1, *f;

	dim = 2;
	lx = xl, ux = xu, Lx = Kx, dx = (ux-lx)/(Kx-1);
	ly = yl, uy = yu, Ly = Ky, dy = (uy-ly)/(Ky-1);
	
	Fx = new double[Lx];
	Fy = new double[Lx*Ly];
	lby = new double[Lx];
	uby = new double[Lx];
	lyi = new long[Lx];
	uyi = new long[Lx];
	
	f = new double[max(Lx,Ly)];
	for (long l = Lx-1; l >= 0; --l) Fx[l] = 0.;
	for (long l = Lx*Ly-1; l >= 0; --l) Fy[l] = 0.;
	
  // Sweep over x-direction
	i0 = 0, i1 = Kx-1, l0 = 0, l1 = i1;
	for (i = i0, x = lx; i < Kx; ++i, x += dx) {
	  lby[i] = lfy(x);
		uby[i] = ufy(x);

		// Sweep over y-direction
		if (uby[i] >= lby[i]) {               // =  NB !

			m0 = Ky*i + (j0 = (int)ceil((lby[i]-ly)/dy));
			m1 = Ky*i + (j1 = (int)floor((uby[i]-ly)/dy));
			lyi[i] = j0, uyi[i] = j1;
			
			// Sample pdf
			y = ly+dy*j0;
			f0 = (y > lby[i]) ? fxy(x,lby[i]) : 0.;
			for (j = j0; j <= j1; ++j, y += dy) f[j] = fxy(x,y);
			y = ly+dy*j1;
			f1 = (y < uby[i]) ? fxy(x,uby[i]) : 0.;
			
			// Evaluate cumulative Fy(y|x) = \int^y fxy(x,y) dy
			Fx[i] = clint(j0,j1,m0,m1, lby[i],uby[i],ly,dy, f0,f1,f, Fy);

		}	 //  endif lby < uby		
/*
		if (i%10 == 0) {		
			cout<<i<<" "<<m0<<"  "<<f0<<" "; for (j = 0; j < 101; j += 10) cout<<f[j0+j]<<" "; 
			cout<<f1<<"   "<<lby[i]<<" "<<uby[i]<<endl;
			for (j = 0; j < 101; j += 10) cout<<Fy[m0+j]<<" "; 
			cout<<" "<<Fx[i]*dy<<endl;
		}
*/
	}  //  endfor i

  // Evalutate cumulative Fx(x) = \int^x(\int\int fxy(x,y) dy)dx
	f0 = f1 = 0.;
	for (i = i0; i <= i1; ++i) f[i] = Fx[i];
	fnorm = clint(i0,i1,l0,l1, lx,ux,lx,dx, f0,f1,f, Fx)*dx*dy;
/*	
	cout<<"Fx "<<endl;
	for (i = i0; i <= i1; i += 10) cout<<Fx[i]<<" "; cout<<endl;
	cout<<"Fy "<<endl;
	for (j = Ly-1; j >= 0; j -= 10) {
		cout<<j<<"  ";
		for (i = i0; i <= i1; i += 10) cout<<Fy[Ly*i+j]<<" "; cout<<endl;
	}
*/
	delete [] f;
}



Inversion::Inversion(const double xl, const double xu, const int Kx, 
										 const double yl, const double yu, const int Ky,
										 const double zl, const double zu, const int Kz,
										 double (*lfy)(double), double (*ufy)(double), 
										 double (*lfz)(double,double), double (*ufz)(double,double), 
										 double (*fxyz)(double, double, double), double& fnorm)
{

/*
 *	Prepare draw of r.v. (x,y,z) in fxyz(x,y,z) by direct numerical inversion
 *	Maximum (compact) xyz volume: (xl,xu) x (yl,yu) x (zl,zu)
 *					divided in Kx x Ky x Kz cells
 *	Lower and upper boundary for y: lfy(x) and ufy(x)
 *	Lower and upper boundary for z: lfz(x,y) and ufz(x,y)
 *  Method: draw x from fx(x) = \int fxy(x,y) dy in (lx,ux),
 *          y given x from fxy(y|x) = \int fxyz(x,y,z) dz in (lfy(x),ufy(x))
 *          and z given x,y from fxyz(x,y,z) in (lfy(x),ufy(x))x(lfz(x,y),ufz(x,y)
 *	Returns fnorm = \int \int \int fxyz(x,y,z) dzdydx
 */

  int i, i0, i1, j, j0, j1, k, k0, k1;
	long l0, l1, m, m0, m1, n0, n1;
	double x, y, z, f0, f1, *f;

	dim = 3;
	lx = xl, ux = xu, Lx = Kx, dx = (ux-lx)/(Kx-1);
	ly = yl, uy = yu, Ly = Ky, dy = (uy-ly)/(Ky-1);
	lz = zl, uz = zu, Lz = Kz, dz = (uz-lz)/(Kz-1);

	Fx = new double[Lx];
	Fy = new double[Lx*Ly];
	Fz = new double[Lx*Ly*Lz];
	lby = new double[Lx];
	uby = new double[Lx];
	lbz = new double[Lx*Ly];
	ubz = new double[Lx*Ly];
	lyi = new long[Lx];
	uyi = new long[Lx];
	lzi = new long[Lx*Ly];
	uzi = new long[Lx*Ly];
	
	f = new double[max(Lz,max(Lx,Ly))];
	for (long l = Lx-1; l >= 0; --l) Fx[l] = 0.;
	for (long l = Lx*Ly-1; l >= 0; --l) Fy[l] = 0.;
	for (long l = Lx*Ly*Lz-1; l >= 0; --l) Fz[l] = 0.;
	
	// Sweep over x-direction
	i0 = 0, i1 = Kx-1, l0 = i0, l1 = i1;
	for (i = i0, x = lx; i <= i1; ++i, x += dx) {
	  lby[i] = lfy(x);
		uby[i] = ufy(x);

    // Sweep over y-direction
		if (uby[i] > lby[i]) {

			m0 = Ky*i + (j0 = (int)ceil((lby[i]-ly)/dy));
			m1 = Ky*i + (j1 = (int)floor((uby[i]-ly)/dy));
			lyi[i] = j0, uyi[i] = j1;
			
			for (j = j0, m = m0, y = ly+dy*j0; j <= j1; ++j, ++m, y += dy) {
				lbz[m] = lfz(x,y);
				ubz[m] = ufz(x,y);
				
				// Sweep over z-direction
				if (ubz[m] > lbz[m]) {

					n0 = Kz*m + (k0 = (int)ceil((lbz[m]-lz)/dz));
					n1 = Kz*m + (k1 = (int)floor((ubz[m]-lz)/dz));
					lzi[m] = k0, uzi[m] = k1;
					
					// Sample pdf
					z = lz+dz*k0;
					f0 = (z > lbz[m]) ? fxyz(x,y,lbz[m]) : 0.;
					for (k = k0; k <= k1; ++k, z += dz) f[k] = fxyz(x,y,z);
					z = lz+dz*k1;
					f1 = (z < ubz[m]) ? fxyz(x,y,ubz[m]) : 0.;
			
					// Evaluate cumulative Fz(z|x,y) = \int^z fxyz(x,y,z) dz
					Fy[m] = clint(k0,k1,n0,n1, lbz[m],ubz[m],lz,dz, f0,f1,f, Fz);

				}	 //  endif lbz < ubz				
			}  //  endfor j

			// Evalutate cumulative pf Fy(y|x) = \int^y(\int fxyz(x,y,z)dz)dy
			f0 = f1 = 0.;
			for (j = j0, m = m0; j <= j1; ++j, ++m) f[j] = Fy[m];
			Fx[i] = clint(j0,j1,m0,m1, lby[i],uby[i],ly,dy, f0,f1,f, Fy);

		}  // endif lby < uby
	}  // endfor i

  // Evalutate cumulative Fx(x) = \int^x(\int\int fxyz(x,y,z)dzdy)dx
	f0 = f1 = 0.;
	for (i = i0; i <= i1; ++i) f[i] = Fx[i];
	fnorm = clint(i0,i1,l0,l1, lx,ux,lx,dx, f0,f1,f, Fx)*dx*dy*dz;
/*	
	Fx[i0] = 0.;
	if (uby[i0] == lby[i0] && (3.*(f[i0]-f[i0+1])+f[i0+2]) > f[i0]/3.) // inf deriv
		Fx[i0+1] = 2./3.*f[i0];
	else 
		Fx[i0+1] = c9*f[i0]+c19*f[i0+1]-c5*f[i0+2]+c1*f[i0+3];
	for (i = i0+2; i < i1; ++i)
		Fx[i] = Fx[i-1]-c1*(f[i-2]+f[i+1])+c13*(f[i]+f[i-1]);
	if (uby[i1] == lby[i1] && (3.*(f[i1]-f[i1-1])+f[i1-2]) > f[i1]/3.) // inf deriv
		Fx[i1] = Fx[i1-1]+2./3.*f[i1];
	else 
		Fx[i1] = Fx[i1-1]+c1*f[i1-3]-c5*f[i1-2]+c19*f[i1-1]+c9*f[i1];
	
	// Normalize the cumulative F_x(x) distribution
	f[0] = Fx[i1];
	fnorm = f[0]*dx*dy*dz;
	for (i = i0; i <= i1; ++i) Fx[i] /= f[0];
*/

  delete [] f;
}


Inversion::~Inversion()
{
	delete [] Fx;
	if (dim >= 2) {
		delete [] Fy;
		delete [] lby;
		delete [] uby;
		delete [] lyi;
		delete [] uyi;
	}
	if (dim == 3) {
		delete [] Fz;
		delete [] lbz;
		delete [] ubz;
		delete [] lzi;
		delete [] uzi;
	}
}

double Inversion::clint(const int k0, const int k1, const long n0, const long n1,
								const double lb, const double ub, const double lv, const double dv, 
								const double f0, const double f1, const double* f, double* F)
{
	// Evaluate cumulative line integral F(v) = \int^v f(v) dv
	
	int k;
	long n;
	double FF, tdv = 2.*dv;
	const double c1 = 1./24., c2 = 2./24., c4 = 4./24., c5 = 5./24., c9 = 9./24., 
							 c10 = 10./24, c13 = 13./24., c16 = 16./24., c19 = 19./24.;
	
	// NB sjekk k0 = k1
	
	if (k0 > k1) FF = (f0+f1)*(ub-lb)/tdv;
	else {
		F[n0] = (f0+f[k0])*(lv+dv*k0-lb)/tdv;
		if (k1-k0 > 2) {
			F[n0+1] = F[n0]+c9*f[k0]+c19*f[k0+1]-c5*f[k0+2]+c1*f[k0+3];
			for (k = k0+2, n = n0+2; k < k1; ++k, ++n)
				F[n] = F[n-1]-c1*(f[k-2]+f[k+1])+c13*(f[k]+f[k-1]);
			F[n1] = F[n1-1]+c9*f[k1]+c19*f[k1-1]-c5*f[k1-2]+c1*f[k1-3];
		} else if (k1-k0 == 2) {
			F[n0+1] = Fz[n0]+c10*f[k0]+c16*f[k0+1]-c2*f[k1];
			F[n1] = Fz[n1-1]+c10*f[k1]+c16*f[k0+1]-c2*f[k0];
		} else if (k1-k0 == 1) F[n1] = F[n0]+(f[k0]+f[k1])/2.;
	}
	FF = F[n1]+(f1+f[k1])*(ub-lv-dv*k1)/tdv;

	for (n = n0; n <= n1; ++n) F[n] /= FF;

//	if (F[n0] > 0.) F[n0-1] = 0.;
//	if (F[n1] < 1.) F[n1+1] = 1.;
					
	return FF;
}
				



double Inversion::interpolate(const double u)
{
	long i, il, iu, i0;
  double Fl, Fu, x, x0, xx[4], FF[4];
	
	il = 0, Fl = 0.;
	iu = Lx-1, Fu = 1.;
	
	while (iu-il > 1) {
		i = (il+iu)/2;
		if (Fx[i] < u) {
			il = i, Fl = Fx[i];
		} else {
			iu = i, Fu = Fx[i];
		}
	}
	x0 = lx+(il+(u-Fl)/(Fu-Fl))*dx;

	i0 = (il > 0) ? il-1 : il;
	if (iu == Lx-1) i0--;
	for (i = 0; i < 4; ++i) {
		xx[i] = lx+(i0+i)*dx, FF[i] = Fx[i0+i];
	}
	x = NRLagrangePol(xx,FF, 4, 1.E-4, u, x0);
//	cout<<u<<" "<<x<<" "<<jl<<" "<<ju<<" "<<Fl<<" "<<Fu<<" "<<-.5+sqrt(.25+2.*u)<<endl;
	
  return x0;
}



double Inversion::interpolate(const double u, const double x)
{
/*
  int i0, i1;
  long kl, ku, k, m0, m1;
	double Fl, Fu, F, y, dx0, dx1, b;
	
	dx0 = (y = (x-lx)/dx)-(i0 = int(y));
	dx1 = 1.-dx0, i1 = i0+1, m0 = Ly*i0, m1 = Ly*i1;
	
	kl = max(lyi[i0],lyi[i1]);
	ku = min(uyi[i0],uyi[i1]);
	
	if (ku >= kl) {
		Fl = dx1*Fy[m0+kl]+dx0*Fy[m1+kl];
		Fu = dx1*Fy[m0+ku]+dx0*Fy[m1+ku];
		if (u < Fl) {
			b = dx1*lby[i0]+dx0*lby[i1];
			y = b+u/Fl*(ly+kl*dy-b);
		} else if (u >= Fu) {
			b = dx1*uby[i0]+dx0*uby[i1];
			y = b-(1.-u)/(1.-Fu)*(b-ly-ku*dy);
		} else {
			while (ku-kl > 1) {
				k = (kl+ku)/2;
				F = dx1*Fy[m0+k]+dx0*Fy[m1+k];
				if (F < u) {
					kl = k, Fl = F;
				} else {
					ku = k, Fu = F;
				}
			}
			y = ly+(kl+(u-Fl)/(Fu-Fl))*dy;
			// NB
		}
	} else {
		Fu = dx1*uby[i0]+dx0*uby[i1];
		Fl = dx1*lby[i0]+dx0*lby[i1];
		y = Fu+u*(Fu-Fl);
	}
*/

	// new version
  int i0, i1, i, jl, ju, j, l;
  long kl, ku, m, m0, m1;
	double Fl, Fu, yy[2], dx0, dx1, y;
	
	dx0 = (y = (x-lx)/dx)-(i = i0 = int(y));
	dx1 = 1.-dx0, i1 = i0+1, m = m0 = Ly*i0, m1 = Ly*i1;
	for (l = 0; l < 2; ++l) {
		jl = lyi[i], Fl = Fy[m+jl];
		ju = uyi[i], Fu = Fy[m+ju];

		if (ju < jl) yy[l] = (1.-u)*lby[i]+u*uby[i];
		else if (Fl >= u) yy[l] = lby[i]+u/Fl*(ly+jl*dy-lby[i]);
		else if (Fu <= u) yy[l] = uby[i]-(1.-u)/(1.-Fu)*(uby[i]-ly-ju*dy);
		else {
			while (ju-jl > 1) {
				j = (jl+ju)/2;
				if (Fy[m+j] < u) {
					jl = j, Fl = Fy[m+j];
				} else {
					ju = j, Fu = Fy[m+j];
				}
			}
			yy[l] = ly+(jl+(u-Fl)/(Fu-Fl))*dy;
		}
		if (yy[l] > uby[i] || yy[l] < lby[i])
			cout<<"out "<<lby[i]<<" "<<uby[i]<<" "<<yy[l]<<" "<<Fl<<" "<<Fu<<" "<<u<<" "<<ly+jl*dy<<" "<<ly+ju*dy
			<<" "<<jl<<" "<<ju<<endl;
//		cout<<"l "<<l<<" "<<yy[l]<<" "<<jl<<" "<<ju<<endl;
//		for (int k = 0; k < Ly; ++k) cout<<Fy[m+k]<<" "; cout<<endl;
		i = i1, m = m1;
	}
	y = yy[0]*dx1+yy[1]*dx0;
//	cout<<y<<" "<<dx0<<" "<<dx1<<endl;
	
	return y;
}
	


double Inversion::interpolate(const double u, const double x, const double y)
{
/*
  int i0, i1, j0, j1;
  long kl, ku, k, m0, m1, n0, n1, n2, n3;
	double Fl, Fu, F, z, dx0, dx1, dy0, dy1, d0, d1, d2, d3, b;
	
	dx0 = (z = (x-lx)/dx)-(i0 = int(z));
	dx1 = 1.-dx0, i1 = i0+1, m0 = Ly*i0, m1 = Ly*i1;
	dy0 = (z = (y-ly)/dy)-(j0 = int(z));
	dy1 = 1.-dy0, j1 = j0+1;
	n0 = Lz*(m0+j0), n1 = Lz*(m0+j1), n2 = Lz*(m1+j0), n3 = Lz*(m1+j1);
	d0 = dx1*dy1, d1 = dx1*dy0, d2 = dx0*dy1, d3 = dx0*dy0;
	
	kl = max(max(lzi[m0+j0],lzi[m0+j1]),max(lzi[m1+j0],lzi[m1+j1]));
	ku = min(min(uzi[m0+j0],uzi[m0+j1]),min(uzi[m1+j0],uzi[m1+j1]));
	
	if (ku >= kl) {
		Fl = d0*Fz[n0+kl]+d1*Fz[n1+kl]+d2*Fz[n2+kl]+d3*Fz[n3+kl];
		Fu = d0*Fz[n0+ku]+d1*Fz[n1+ku]+d2*Fz[n2+ku]+d3*Fz[n3+ku];
		if (u < Fl) {
			b = d0*lbz[m0+j0]+d1*lbz[m0+j1]+d2*lbz[m1+j0]+d3*lbz[m1+j1];
			z = b+u/Fl*(lz+kl*dz-b);
		} else if (u >= Fu) {
			b = d0*ubz[m0+j0]+d1*ubz[m0+j1]+d2*ubz[m1+j0]+d3*ubz[m1+j1];
			z = b-(1.-u)/(1.-Fu)*(b-lz-ku*dz);
		} else {
			while (ku-kl > 1) {
				k = (kl+ku)/2;
				F = d0*Fz[n0+k]+d1*Fz[n1+k]+d2*Fz[n2+k]+d3*Fz[n3+k];
				if (F < u) {
					kl = k, Fl = F;
				} else {
					ku = k, Fu = F;
				}
			}
			z = lz+(kl+(u-Fl)/(Fu-Fl))*dz;
			// NB
//	cout<<"z "<<x<<" "<<y<<" "<<kl<<" "<<ku<<" "<<u<<" "<<Fl<<" "<<Fu<<" "<<z<<" "<<-.5+sqrt(.25+2.*u)<<endl;
		}
	} else {
		Fu = d0*ubz[m0+j0]+d1*ubz[m0+j1]+d2*ubz[m1+j0]+d3*ubz[m1+j1];
		Fl = d0*lbz[m0+j0]+d1*lbz[m0+j1]+d2*lbz[m1+j0]+d3*lbz[m1+j1];
		z = Fu+u*(Fu-Fl);
	}
*/

  // new version
  int i0, i1, j0, j1;
  long kl, ku, k, m0, m1, n[4];
	double Fl, Fu, F, z, dx0, dx1, dy0, dy1, d0, d1, d2, d3;
	
	dx0 = (z = (x-lx)/dx)-(i0 = int(z));
	dx1 = 1.-dx0, i1 = i0+1, m0 = Ly*i0, m1 = Ly*i1;
	dy0 = (z = (y-ly)/dy)-(j0 = int(z));
	dy1 = 1.-dy0, j1 = j0+1;
	n[0] = Lz*(m0+j0), n[1] = Lz*(m0+j1), n[2] = Lz*(m1+j0), n[3] = Lz*(m1+j1);
	d0 = dx1*dy1, d1 = dx1*dy0, d2 = dx0*dy1, d3 = dx0*dy0;
	
	
	return z;
}
	



double Inversion::interpolate(const double u, const double* F, const long ns,
									const double lb, const double ub, const double lv, const double dv)
// u: r.v., F: cumul. dist., ns: starting index, lb: lower bound, ub: upper bound, 
// lv: absolute lower boundary (corresponding to idex ns), dv: step length
{
  long n, nl, nu;
	double Fl, Fu, v;
	
	nl = ns+(long)ceil((lb-lv)/dv);
	nu = ns+(long)floor((ub-lv)/dv);
	
	if (nu >= nl) {
		Fl = F[nl], Fu = F[nu];
		if (Fl >= u) v = lb+u/Fl*((nl-ns)*dv-lb);
		else if (Fu <= u) v = ub-(1.-u)/(1.-Fu)*(ub-(nu-ns)*dv);
		else {
			while (nu-nl > 1) {
				n = (nl+nu)/2;
				if (F[n] < u) { 
					nl = n, Fl = F[n]; 
				}	else {
					nu = n, Fu = F[n];
				}
			}
			v = lv+dv*(nl-ns+(u-Fl)/(Fu-Fl));
		}
	} else v = (1.-v)*lb+u*ub;
	
	return v;
}
	


void Inversion::draw(double& x)
{
  // Draw single r.v. x in fx(x) in (lx,ux)
	
	x = interpolate(drand48());
}

void Inversion::draw(const long N, double* x, const bool quietx)
{ 
  // Draw N rv x in fx n (lx, ux)
	// if (quietx) quiet x else noisy x
	
	int i, il, iu, i0;
	long n;
	double du, u, xu, x0, xx[4], FF[4];
	
	if (quietx) {

		du = 1./N, iu = 1, xu = lx+dx;
		for (n = 0, u = .5*du; n < N; ++n, u += du) {
			while (Fx[iu] <= u) {
				++iu, xu += dx;
			}
			il = iu-1;
			x0 = xu-(Fx[iu]-u)/(Fx[iu]-Fx[il])*dx;
			i0 = (il > 0) ? il-1 : il;
			if (iu == Lx-1) i0--;
			for (i = 0; i < 4; ++i) {
				xx[i] = lx+(i0+i)*dx, FF[i] = Fx[i0+i];
			}
			x[n] = NRLagrangePol(xx, FF, 4, 1.E-4, u, x0);
		}

	} else for (n = 0; n < N; ++n) x[n] = interpolate(drand48());
}		


void Inversion::draw(double& x, double& y)
{ 
  // Draw rv (x,y) in fxy in (lx,ux)x(lby(x),uby(x))
	
	x = interpolate(drand48());
  y = interpolate(drand48(), x);
}


void Inversion::draw(const long N, double* x, double* y, const bool quietx)
{
  // Draw N rv (x,y) in fxy in (lx,ux)x(lby(x),uby(x))
	// if (quietx) then quiet x else noisy x
	
	int i, il, iu, i0;
	long n;
	double du, u, xu, x0, xx[4], FF[4];

	if (quietx) {

		du = 1./N, iu = 1, xu = lx+dx;
		for (n = 0, u = .5*du; n < N; ++n, u += du) {
			while (Fx[iu] <= u) {
				++iu, xu += dx;
			}
			il = iu-1;
			x0 = xu-(Fx[iu]-u)/(Fx[iu]-Fx[il])*dx;
			i0 = (il > 0) ? il-1 : il;
			if (iu == Lx-1) i0--;
			for (i = 0; i < 4; ++i) {
				xx[i] = lx+(i0+i)*dx, FF[i] = Fx[i0+i];
			}
			x[n] = NRLagrangePol(xx, FF, 4, 1.E-4, u, x0);
			y[n] = interpolate(drand48(), x[n]);
		}

	} else {
	
		for (n = 0; n < N; ++n) {
			x[n] = interpolate(drand48());
			y[n] = interpolate(drand48(), x[n]);
		}
	}
}



void Inversion::draw(double& x, double& y, double& z)
{ 
  // Draw rv (x,y,z) in fxyz in (lx,ux)x(lby(x),uby(x))x(lbz(x,y),ubz(x,y))
	
	x = interpolate(drand48());
	y = interpolate(drand48(), x);
	z = interpolate(drand48(), x, y);
}


void Inversion::draw(const long N, double* x, double* y, double* z, 
										 const bool quietx)
{
  // Draw N rv (x,y,z) in fxyz in (lx,ux)x(lby(x),uby(x))x(lbz(x,y),ubz(x,y))
	// if (quietx) then quiet x else noisy x
	
	int i, il, iu, i0;
	long n;
	double du, u, xu, x0, xx[4], FF[4];

	if (quietx) {

		du = 1./N, iu = 1, xu = lx+dx;
		for (n = 0, u = .5*du; n < N; ++n, u += du) {
			while (Fx[iu] <= u) {
				++iu, xu += dx;
			}
			il = iu-1;
			x0 = xu-(Fx[iu]-u)/(Fx[iu]-Fx[il])*dx;
			i0 = (il > 0) ? il-1 : il;
			if (iu == Lx-1) i0--;
			for (i = 0; i < 4; ++i) {
				xx[i] = lx+(i0+i)*dx, FF[i] = Fx[i0+i];
			}
			x[n] = NRLagrangePol(xx, FF, 4, 1.E-4, u, x0);
			y[n] = interpolate(drand48(), x[n]);
			z[n] = interpolate(drand48(), x[n], y[n]);
		}

	} else {
	
		for (n = 0; n < N; ++n) {
			x[n] = interpolate(drand48());
			y[n] = interpolate(drand48(), x[n]);
			z[n] = interpolate(drand48(), x[n], y[n]);
		}
	}
}




double LagrangePol(const double* x, const double* f, const int n,
									 const double X)
{
  // Lagrange interpolation
	// returns P(X) with P(x[i]) = f[i], i = 0,...,n-1.
	
	int i, j, jn, ni;
	double xi[n], PP, P = 0.;
	
	for (i = 0; i < n; ++i) xi[i] = X-x[i];
	for (i = 0; i < n; ++i) {
		ni = n+i, PP = f[i];
		for (j = i+1; j < ni; ++j) {
			jn = j%n, PP *= xi[jn]/(x[i]-x[jn]);
		}
		P += PP;
	}
	
	return P;
}



double NRLagrangePol(const double* x, const double* f, const int n, 
										 const double eps, const double u, const double X0)
{
	// Newton-Ralphson solver of Lagrange polynom equation, first guess X0
	// returns X; P(X) = u with P(x[i]) = f[i], i = 0,...,n-1, accuracy eps

  int i, j, ni;
	double AA[n], A[n], B[n], dX, P, Pd, F, X = X0, z[2*n-1];

  if (n > 5) {
		cout<<"NRLagrangePol: n>5\n"; exit(1);
	}
	
  for (i = 0; i < n; ++i) {
		z[i] = x[i], A[i] = 0.;
	}
	for (i = n; i < 2*n-1; ++i) z[i] = x[i%n];
	
  for (i = 0; i < n; ++i) {
		ni = n+i, F = f[i], AA[0] = AA[n-1] = 1.;
		for (j = 1; j < n-1; ++j) AA[j] = 0.;
		for (j = i+1; j < ni; ++j) {
			F /= (z[i]-z[j]);
			AA[1] += z[j], AA[n-1] *= z[j];
			if (n > 3) for (int k = j+1; k < ni; ++k) {
				AA[2] += z[j]*z[k];
				if (n == 5) for (int l = k+1; l < ni; ++l) AA[3] += z[j]*z[k]*z[l];
			}
		}

		for (j = 0; j < n; ++j) {
			A[j] += F*AA[j], B[j] = (n-1-j)*A[j], F = -F;
		}
	}
	
	do {
			P = A[0]*X+A[1], Pd = B[0];
			for (i = 2; i < n; ++i) {
				P = P*X+A[i],	Pd = Pd*X+B[i-1];
			}
			dX = (P-u)/Pd, X -= dX;
	} while (abs(dX) > eps);
	return X;
}



#ifdef TEST
// test programs

double fx(double x);
double fx(double x)
{
  return 1.+2.*x;
}

double fxy(double x, double y);
double fxy(double x, double y)
{
  return (1.+x)*(1.+y);
}

double fxyz(double x, double y, double z);
double fxyz(double x, double y, double z)
{
  return (1.+2.*x)*(1.+2.*y)*(1.+2.*z);
}

double lbz(double x, double y);
double lbz(double x, double y)
{
  return 0.;
}


double ubz(double x, double y);
double ubz(double x, double y)
{
  return 1.;
}

double lby(double x);
double lby(double x)
{
  return -1;
}

double uby(double x);
double uby(double x)
{
  return -(1.-x)/2.;
}


double lb(double x);
double lb(double x)
{
  return -1.;
}

double ub(double x);
double ub(double x)
{
  return 1.;
}

double phi(double x);
double phi(double x)
{
  return x-sin(x);
}

double fea(double x, double y);
double fea(double x, double y)
{
  return exp(-.5*y*y+phi(x));
}

double lea(double x);
double lea(double x)
{
  return sqrt(2.*phi(x));
}

double uea(double x);
double uea(double x)
{
  return 5.;
}

double fed(double x, double y);
double fed(double x, double y)
{  
	return exp(-.5*y*y+phi(x)-4.*atan(1.));
}

double led(double x);
double led(double x)
{
  return -5.;
}

double ued(double x);
double ued(double x)
{
  return -sqrt(2.*phi(x));
}

double fer(double x, double y);
double fer(double x, double y)
{  
	return exp(-.5*y*y+phi(x)-4.*atan(1.));
}

double ler(double x);
double ler(double x)
{
  return -sqrt(2.*phi(x));
}

double uer(double x);
double uer(double x)
{
  return 0.;
}

double circle(double x, double y);
double circle(double x, double y)
{
  return 1.;
}

double lc(double x);
double lc(double x)
{
  double y = 1.-x*x;
	if (y <= 0.) return 0.; else return -sqrt(y);
}

double uc(double x);
double uc(double x)
{
  double y = 1.-x*x;
  if (y <= 0.) return 0.; else return sqrt(y);
}



int main()
{
	double fnorm, *x, *y, *z, xx, yy, zz;
	
	long N = 50000;
	x = new double[N];
	y = new double[N];
	z = new double[N];

	double tpi = 4.*atan(1.);

	Histogram enhist(0., 1., 100);
/*	
	Inversion endim(0., 1., 5, fx, fnorm);
	cout<<"endim "<<fnorm<<endl;

  endim.draw(xx);
	cout<<xx<<endl;
	endim.draw(N, x);
	enhist.record(x, 0, 1, N);
	enhist.print();
*/

	Inversion todim(0.,1.,101, -1.,0.,101,lby, uby, fxy, fnorm);

	Histogram xhist(0., 1., 10);
	xhist.reset();

	todim.draw(N, x, y);
	xhist.record(x, 0, 1, N);
	xhist.print("x histogram");
	
	Histogram yhist(-1., 0., 10);
	yhist.reset(); 
	yhist.record(y, 0, 1, N);
	yhist.print("y histogram");
	
	Histogram xyhist(0.,1.,10, -1.,0.,10);
	double *xy = new double(2*N);
	for (long n = 0; n < N; ++n) {
		xy[2*n] = x[n], xy[2*n+1] = y[n];
	}
	xyhist.reset();
	xyhist.record(xy, 0, 1, 2, N);
	xyhist.print("xy-hist");

/*
	Inversion tredim(0.,1.,101, 0.,1.,101, 0.,1.,101, lby, uby,
									lbz, ubz, fxyz, fnorm);
	cout<<"tredim "<<fnorm<<endl;
	
	tredim.draw(xx, yy, zz);

	tredim.draw(N, x, y, z);
	enhist.reset();
	enhist.record(x, 0, 1, N);
	cout<<"x histogram"<<endl;
	enhist.print();
	
	enhist.reset();
	enhist.record(y, 0, 1,N);
	cout<<"y histogram"<<endl;
	enhist.print();

	enhist.reset();
	enhist.record(z, 0, 1,N);
	cout<<"z histogram"<<endl;
	enhist.print();

									

	Inversion circ(-1., 1., 101, -1., 1., 101, lc, uc, circle, fnorm);
	cout<<fnorm<<' '<<tpi<<endl;
*/

/*
	Inversion ea(0., tpi, 100, 0., 5., 100, lea, uea, fea, fnorm);
	ea.draw(N, x, y);
	cout<<"fea "<<fnorm<<endl;
	
	Histogram hist(0., tpi, 100);
	hist.record(x, 0, 1, N);
	hist.print();

	Inversion ed(0., tpi, 100, -5., 0., 100, led, ued, fed, fnorm);
	ed.draw(N, x, y);
	cout<<"fed "<<fnorm<<endl;
	
	hist.reset();
	hist.record(x, 0, 1, N);
	hist.print();

	Inversion er(0., tpi, 100, -5., 0., 100, ler, uer, fer, fnorm);
	er.draw(N, x, y);
	cout<<"fer "<<fnorm<<endl;
	
//	x[0] += 10.;
	hist.reset();
	hist.record(x, 0, 1, N);
	hist.print();
*/
	
	return 0;
}

#endif


