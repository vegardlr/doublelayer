/*
 *  MultiGrid.cpp
 *  pic3D
 *
 *  Created by Jan Trulsen on 05/01/06.
 *  Copyright 2006 __UiO__. All rights reserved.
 *
 */

#include <cstdio>
//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include "Utils.h"
#include "MultiGrid.h"

using namespace std;

// Multigrid solver for the linear or non-linear Poisson equations
//  del^2 u = rho     or    del^2 u + fNZP(u, dfNZP) = rho
// where fNZP(const double* u, double &dfNZP) is a user supplied function
// with the auxilary return value dfNZP = d/du fNZP().
// For the non-linear case define flags MGFAS and NZP,
// for the linear case define flag MGLIN or chose non-linear obtion with fNZP(..) = 0.
// The equations are solved in one, two or three dimensional space according to the 
// value of the global parameter (flag) XDIM = 1, 2 or 3 over domene with regular grid,
// grid spacing h[i], size L[i] = n[i]*h[i] with grid positions X[1] = 0, X[n[0]] = L[0]
// and similar for other directions. Total grid dimensions are n[i]+2 in each direction.
// Multigrid depth ngrid determined by 2^ngrid = LCF(n[i]-1).

// Boundary type are set by global parameters (flags) BCX0, BCX1, ... with any
// possible combination allowed through values 0: Dirichlet, 1: von Neumann, 2: periodic.
// Actual boundary values are supplied through rho(X[0]) and rho(X[n[0]+1]), that is,
// BCX0 = 0: rho(X[0]) = u(X[1]), BCX0 = 1: rho(X[0]) = d/dx u(X[1]) and similar.
// For the periodic case, note that u(X[1]) = u(X[n[0]]) and similar for other directions.
// Note that rho are modified by boundary conditions on return from mglin/mgfas, but may be reset
// by calling delbc2rho(rho, n, h) where n and h are the full grid size and spacings.

// relax(u, rho, n, h) may be run as an independent red-black Gauss-Seidel iteration but will
// require an initial call to addbc2rho(rho, n, h).


MultiGrid::MultiGrid(const int* n, const double* h)
{
  // Initialize multigrid solver for (non-linear) Poisson equation
#ifdef MGLIN
#ifdef MGFAS
  ErrorMessage("specify either MGLIN or MGFAS");
#endif
#ifdef NZP
  ErrorMessage("do not specify both MGLIN and NZP");
#endif
#endif
	
	int i, j, jj, r, nmin;
	
	// find largest common factor of the different (n[i]-1)
	for (i = 1, nmin = n[0]-1; i < XDIM; ++i) if (n[i]-1 < nmin) nmin	= n[i]-1;
	for (i = 0, r = nmin; i < XDIM; ++i) { 
		while ((r = (n[i]-1)%nmin) > 0) {
			if (r == 1) ErrorMessage("largest common factor of (n-1)'s must be > 1");
			while ((r = nmin%r) > 0) {
				if (r == 1) ErrorMessage("largest common factor of (n-1)'s must be > 1");
				nmin = r;
			}
		}
		r = nmin;
	}
	ngrid = 0; while (r >>= 1) ngrid++;
	if (nmin != 1<<ngrid) ErrorMessage("largest common factor of (n-1)'s must be power of 2");
	if (ngrid >= NGMAX) ErrorMessage("increase NGMAX");

	for (i = 0, nnp[ngrid] = 1; i < XDIM; ++i) { 
	  L[i] = n[i]*h[i],	nn[ngrid][i] = n[i], hh[ngrid][i] = h[i], nnp[ngrid] *= (nn[ngrid][i]+2); 
	}
	for (j = ngrid-1, nnp[j] = 1; j > 0; --j) {
		jj = j+1;
		for (i = 0, nnp[j] = 1; i < XDIM; ++i) { 
			nn[j][i] = nn[jj][i]/2+1, hh[j][i] = 2.*hh[jj][i], nnp[j] *= (nn[j][i]+2); 
		}
	}
				
	for (j = 1; j <= ngrid; ++j) {
	  if (j != ngrid) {
			irho[j] = new double[nnp[j]];
			iu[j] = new double[nnp[j]];
		}
		itmp[j] = new double[nnp[j]];
#ifdef MGFAS
		itau[j] = new double[nnp[j]];
#endif
	}
	
	ncycle = 1;
	npre = 3, npost = 3;
	alpha = .33;
}



MultiGrid::~MultiGrid()
{
	for (int j = ngrid; j >= 1; --j) {
#ifdef MGFAS
		delete [] itau[j];
#endif
		delete [] itmp[j];
		if ( j != ngrid) {
			delete [] irho[j];
			delete [] iu[j];
		}
	}
}


#ifdef MGFAS
// Multigrid non-linear Poisson equation solver
void MultiGrid::mgfas(double *u, double* rho, const int *n, const double *h)
{
	int j, jj, jm, jcycle, jpre, jpost;

	// restrict rho and adjust boundary conditions
	for (j = ngrid-1, nnp[j] = 1; j > 0; --j) {
		jj = j+1;
		rstrct(irho[j], (jj != ngrid ? irho[jj] : rho), nn[j], nn[jj]);   
		addbc2rho((jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);
#ifdef MGPRNT		
		prnt("rho ", nn[jj], hh[jj], (jj != ngrid ? irho[jj] : rho));
#endif
	}
	addbc2rho(irho[1], nn[1], hh[1]);
#ifdef MGPRNT
	prnt("rho ", nn[1], hh[1], irho[1]);
#endif	
				
	slvsml(iu[1], irho[1], nn[1], hh[1]);  // initial solution on coarsest grid
#ifdef MGPRNT
	prnt("slvsml ", nn[1], hh[1], iu[1]);
#endif	

	
	for (j = 2; j <= ngrid; ++j) {
		jm = j-1;
	
		// interpolate from coarse to finer grid
		interp((j != ngrid ? iu[j] : u), iu[jm], (j != ngrid ? irho[j] : rho), nn[j], nn[jm]);

		for (jcycle = 0; jcycle < ncycle; ++jcycle) {  // V-cycle loop
			for (jj = j, jm = jj-1; jj > 1; --jm, --jj) {  // downward stroke of V

				// pre-smoothing
				for (jpre = 0; jpre < npre; ++jpre) 
					relax((jj != ngrid ? iu[jj] : u), (jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);
#ifdef MGPRNT
				if (jj == j) prnt("presmoothing ", nn[jj], hh[jj], (jj != ngrid ? iu[jj] : u));
#endif

				lop(itmp[jj], (jj != ngrid ? iu[jj] : u), nn[jj], hh[jj]);
				rstrct(itmp[jm], itmp[jj], nn[jm], nn[jj]);
				rstrct(iu[jm], (jj != ngrid ? iu[jj] : u), nn[jm], nn[jj]);
				lop(itau[jm], iu[jm], nn[jm], hh[jm]);
				matsub(itau[jm], itau[jm], itmp[jm], nn[jm]);
				if (j == jj) trerr = alpha*anorm(itau[jm], nn[jm]);
				rstrct(irho[jm], (jj != ngrid ? irho[jj] : rho), nn[jm], nn[jj]);
				matadd(irho[jm], irho[jm], itau[jm], nn[jm]);
			}
			
			slvsml(iu[1], irho[1], nn[1], hh[1]);  // bottom of V, solve on coarsest grid
#ifdef MGPRNT
			prnt("bottom V u ", nn[1], hh[1], iu[1]);
#endif
			
			for (jj = 2, jm = jj-1; jj <= j; jm = jj, ++jj) {  // upward stroke of V

				rstrct(itmp[jm], (jj != ngrid ? iu[jj] : u), nn[jm], nn[jj]);
				matsub(itmp[jm], iu[jm], itmp[jm], nn[jm]);
				interp(itau[jj], itmp[jm], (jj != ngrid ? iu[jj] : u), nn[jj], nn[jm]); 
				matadd((jj != ngrid ? iu[jj] : u), (jj != ngrid ? iu[jj] : u), itau[jj], nn[jj]);

				// post-smoothing
 				for (jpost = 0; jpost < npost; ++jpost) 
					relax((jj != ngrid ? iu[jj] : u), (jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);
#ifdef MGPRNT
				if (jj == j) prnt("postsmoothing ", nn[jj], hh[jj], (jj != ngrid ? iu[jj] : u));
#endif
			}

			lop(itmp[j], (j != ngrid ? iu[j] : u), nn[j], hh[j]);  
			matsub(itmp[j], itmp[j], (j != ngrid ? irho[j] : rho), nn[j]);
			res = anorm(itmp[j], nn[j]);
					
//			cout<<"jcycle j nfres trerr "<<jcycle<<" "<<j<<" "<<res<<" "<<trerr<<endl;
			if (res < trerr) break;

		}  // end jcycle
	}  // end j
}



	

void MultiGrid::lop(double *res, double *u, const int *n, const double *h)
{
#ifdef NZP
	double dfNZP;
//	double fNZP(const double u, double &dfNZP);
#endif

#if XDIM == 1
	int i, i0, il;
	double A = 1./(h[0]*h[0]), B = -2.*A;

#if BCX0+BCX1 > 3  // periodic
  i0 = 1, il = n[0]-1;
	u[0] = u[il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 2;
	res[1] = 0.;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 1;
	u[0] = u[2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
	res[n[0]] = 0.;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	u[n[0]+1] = u[n[0]-1];
#endif

  for (i = i0; i <= il; ++i)  // interior points
#ifdef NZP
		res[i] = A*(u[i-1]+u[i+1])+B*u[i]+fNZP(u[i], dfNZP);
#else
		res[i] = A*(u[i-1]+u[i+1])+B*u[i];
#endif
#endif  // XDIM == 1


#if XDIM == 2
	long i, i0, il, ii, j, j0, jl, dj = n[0]+2;
#if BCY0 == 1 || BCY1 == 1
	long dj2 = 2*dj;
#endif
	double A = 1./(h[0]*h[0]), B = 1./(h[1]*h[1]), C = -2.*(A+B);

#if BCX0+BCX1 > 3  // periodic 
  i0 = 1, il = n[0]-1;
	for (j = 1, ii = dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 2;
	for (j = 1, ii = dj+1; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 1;
	for (j = 1, ii = dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
	for (j = 1, ii = dj+n[0]; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	for (j = 1, ii = dj+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii-2];
#endif
#if BCY0+BCY1 > 3  // periodic
  j0 = 1, jl = n[1]-1;
	for (i = 1, j = jl*dj; i <= n[0]; ++i) u[i] = u[i+j];
#endif
#if BCY0 == 0  // Diriclet
	j0 = 2;
	for (i = 1, ii = dj+i; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY0 == 1  // von Neumann
	j0 = 1;
	for (i = 1, ii = i; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dj2];
#endif
#if BCY1 == 0  // Diriclet
	jl = n[1]-1;
	for (i = 1, ii = dj*n[1]+i; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY1 == 1  // von Neumann
	jl = n[1];
	for (i = 1, ii = dj*n[1]+dj+i; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dj2];
#endif

  for (j = j0; j <= jl; ++j)  // interior points
		for (i = i0, ii = dj*j+i; i <= il; ++i, ++ii) 
#ifdef NZP
			res[ii] = A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*u[ii]+fNZP(u[ii], dfNZP);
#else
			res[ii] = A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*u[ii];
#endif
#endif  // XDIM == 2


#if XDIM == 3
	long i, i0, il, ii, j, j0, jl, jj, k, k0, kl, kk, dj = n[0]+2, dk = dj*(n[1]+2);
#if BCY0 == 1 || BCY1 == 1
  long dj2 = 2*dj;
#endif
#if BCZ0 == 1 || BCZ1 == 1
	long dk2 = 2*dk;
#endif
	double A = 1./(h[0]*h[0]), B = 1./(h[1]*h[1]), C = 1./(h[2]*h[2]), D = -2.*(A+B+C);

#if BCX0+BCX1 > 3  // periodic
  i0 = 1, il = n[0]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 2;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+1; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii-2];
#endif
#if BCY0+BCY1 > 3  // periodic
  j0 = 1, jl = n[1]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+1, j = jl*dj; i <= n[0]; ++i, ++ii) u[ii] = u[ii+j];
#endif
#if BCY0 == 0  // Diriclet
	j0 = 2;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY0 == 1  // von Neumann
	j0 = 1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dj2];
#endif
#if BCY1 == 0  // Diriclet
	jl = n[1]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj*n[1]+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY1 == 1  // von Neumann
	jl = n[1];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj*n[1]+dj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dj2];
#endif
#if BCZ0+BCZ1 > 3  // periodic
	k0 = 1, kl = n[2]-1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+1, k = kl*dk; i <= n[0]; ++i, ++ii) u[ii] = u[ii+k];
#endif
#if BCZ0 == 0  // Diriclet
	k0 = 2;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk+jj+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCZ0 == 1  // von Neumann
	k0 = 1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dk2];
#endif
#if BCZ1 == 0  // Diriclet
	kl = n[2]-1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk*n[2]+jj+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCZ1 == 1  // von Neumann
	kl = n[2];
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk*n[2]+dk+jj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dk2];
#endif

  for (k = k0, kk = dk*k; k <= kl; ++k, kk += dk)  // interior points
		for (j = j0, jj = kk+dj*j0; j <= jl; ++j, jj += dj) 
			for (i = i0, ii = jj+i; i <= il; ++i, ++ii) 
#ifdef NZP
				res[ii] = A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*(u[ii-dk]+u[ii+dk])
					+D*u[ii]+fNZP(u[ii],dfNZP);
#else
				res[ii] = A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*(u[ii-dk]+u[ii+dk])+D*u[ii];
#endif
#endif  // XDIM == 3
}





void MultiGrid::matsub(double *c, const double *a, const double *b, const int *n)
{
#if XDIM == 1
  long i;
	
	for (i = 1; i <= n[0]; ++i) c[i] = a[i]-b[i];
#endif


#if XDIM == 2
  long i, j, ii, jj, dj = n[0]+2;
	
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) c[ii] = a[ii]-b[ii];
#endif  // XDIM == 2


#if XDIM == 3
  long i, j, k, ii, jj, kk, dj = n[0]+2, dk = dj*(n[1]+2);
	
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk) 
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) c[ii] = a[ii]-b[ii];
#endif
}
	


void MultiGrid::matadd(double *c, const double *a, const double *b, const int *n)
{
#if XDIM == 1
  int i;
	
	for (i = 1; i <= n[0]; ++i) c[i] = a[i]+b[i];
#endif


#if XDIM == 2
  long i, j, ii, jj, dj = n[0]+2;
	
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) c[ii] = a[ii]+b[ii];
#endif  // XDIM == 2


#if XDIM == 3
  long i, j, k, ii, jj, kk, dj = n[0]+2, dk = dj*(n[1]+2);
	
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk) 
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) c[ii] = a[ii]+b[ii];
#endif
}
	


double MultiGrid::anorm(const double *a, const int *n)
{
	double sum = 0.;
#if XDIM == 1
  int i;
	
	for (i = 1; i <= n[0]; ++i) sum += a[i]*a[i];

	return sqrt(sum)/n[0];
#endif


#if XDIM == 2
  long i, j, ii, jj, dj = n[0]+2;
	
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) sum += a[ii]*a[ii];

		return sqrt(sum)/(n[0]*n[1]);
#endif  // XDIM == 2


#if XDIM == 3
  long i, j, k, ii, jj, kk, dj = n[0]+2, dk = dj*(n[1]+2);
	
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk) 
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) sum += a[ii]*a[ii];

			return sqrt(sum)/(n[0]*n[1]*n[2]);
#endif
}			
#endif    // #ifdef MGFAS



#ifdef MGLIN
// Multigrid linear Poisson equation solver
void MultiGrid::mglin(double *u, double * rho, const int *n, const double *h)
{
	int j, jj, jm, jcycle, jpre, jpost;

	// restrict rho and adjust BC
	for (j = ngrid-1, nnp[j] = 1; j > 0; --j) {
		jj = j+1;
		rstrct(irho[j], (jj != ngrid ? irho[jj] : rho), nn[j], nn[jj]);
		addbc2rho((jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);
#ifdef MGPRNT
		prnt("rho ", nn[jj], hh[jj], (jj != ngrid ? irho[jj] : rho));
#endif
	}
	addbc2rho(irho[1], nn[1], hh[1]);
#ifdef MGPRNT
	prnt("rho ", nn[1], hh[1], irho[1]);
#endif
	fill0(iu[1],nn[1]);			
	slvsml(iu[1], irho[1], nn[1], hh[1]);  // initial solution on coarsest grid
#ifdef MGPRNT
	prnt("slvsml ", nn[1], hh[1], iu[1]);
#endif
	
	for (j = 2; j <= ngrid; ++j) {
		jm = j-1;
		
		// interpolate from coarse to finer grid
		interp((j != ngrid ? iu[j] : u), iu[jm], (j != ngrid ? irho[j] : rho), nn[j], nn[jm]);

		for (jcycle = 0; jcycle < ncycle; ++jcycle) {  // V-cycle loop
			for (jj = j, jm = jj-1; jj > 1; --jm, --jj) {  // downward stroke of V

				// pre-smoothing
				for (jpre = 0; jpre < npre; ++jpre) 
					relax((jj != ngrid ? iu[jj] : u), (jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);
#ifdef MGPRNT
				if (jj == j) prnt("presmoothing ", nn[jj], hh[jj], (jj != ngrid ? iu[jj] : u));
#endif

				resid(itmp[jj], (jj != ngrid ? iu[jj] : u), (jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);

				rstrct(irho[jm], itmp[jj], nn[jm], nn[jj]);
				fill0(iu[jm], nn[jm]);  // initial guess next relaxation
			}
			
			slvsml(iu[1], irho[1], nn[1], hh[1]);  // bottom of V, solve on coarsest grid
			
			for (jj = 2, jm = jj-1; jj <= j; jm = jj, ++jj) {  // upward stroke of V
				// ires used as temporary storage in addint
				addint((jj != ngrid ? iu[jj] : u), iu[jm], itmp[jj], itmp[jj], nn[jj], nn[jm]);

				// post-smoothing
				for (jpost = 0; jpost < npost; ++jpost) 
					relax((jj != ngrid ? iu[jj] : u), (jj != ngrid ? irho[jj] : rho), nn[jj], hh[jj]);
#ifdef MGPRNT
				if (jj == j) prnt("postsmoothing ", nn[jj], hh[jj], (jj != ngrid ? iu[jj] : u));				
#endif
			}
		}
	}
}
			



void MultiGrid::addint(double *u, const double *uc, double *res, 
	const double *rho, const int *n, const int *nc)
{
  interp(res, uc, rho, n, nc);
#if XDIM == 1
	int i;
	
	for (i = 1; i <= n[0]; ++i) u[i] += res[i];
#endif  // XDIM == 1

#if XDIM == 2
  long i, j, ii, dj = n[0]+2;
	
	for (j = 1; j <= n[1]; ++j)
		for (i = 1, ii = dj*j+i; i <= n[0]; ++i, ++ii) u[ii] += res[ii];
#endif  // XDIM == 2

#if XDIM == 3
  long i, j, k, ii, jj, dj = n[0]+2, dk = dj*(n[1]+2);
	
	for (k = 1; k <= n[2]; ++k)
		for (j = 1, jj = dk*k; j <= n[1]; ++j)
			for (i = 1, ii = jj+dj*j+i; i <= n[0]; ++i, ++ii) u[ii] += res[ii];
#endif  // XDIM == 3
}



void MultiGrid::resid(double *res, double *u, const double *rhs, const int *n, const double *h)
{
#if XDIM == 1
	int i, i0, il;
	double A = 1./(h[0]*h[0]), B = 2.*A;

#if BCX0+BCX1 > 3  // periodic
  i0 = 1, il = n[0]-1;
	u[0] = u[il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 2;
	res[1] = 0.;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 1;
	u[0] = u[2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
	res[n[0]] = 0.;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	u[n[0]+1] = u[n[0]-1];
#endif

  for (i = i0; i <= il; ++i)  // interior points
		res[i] = -A*(u[i-1]+u[i+1])+B*u[i]+rhs[i];
#endif  // XDIM == 1


#if XDIM == 2
	long i, i0, il, ii, j, j0, jl, dj = n[0]+2, dj2 = 2*dj;
	double A = 1./(h[0]*h[0]), B = 1./(h[1]*h[1]), C = 2.*(A+B);

#if BCX0+BCX1 > 3  // periodic
	i0 = 1, il = n[0]-1;
	for (j = 1, ii = dj*j; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 2;
	for (j = 1, ii = dj*j+1; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 1;
	for (j = 1, ii = dj*j; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
	for (j = 1, ii = dj*j+n[0]; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	for (j = 1, ii = dj*j+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii-2];
#endif
#if BCY0+BCY1 > 3  // periodic
	j0 = 1, jl = n[1]-1;
	for (i = 1, ii = i, j = jl*dj; i <= n[0]; ++i, ++ii) u[ii] = u[ii+j];
#endif
#if BCY0 == 0  // Diriclet
	j0 = 2;
	for (i = 1, ii = dj+i; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY0 == 1  // von Neumann
	j0 = 1;
	for (i = 1, ii = i; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dj2];
#endif
#if BCY1 == 0  // Diriclet
	jl = n[1]-1;
	for (i = 1, ii = dj*n[1]+i; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY1 == 1  // von Neumann
	jl = n[1];
	for (i = 1, ii = dj*n[1]+dj+i; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dj2];
#endif

  for (j = j0; j <= jl; ++j)  // interior points
		for (i = i0, ii = dj*j+i; i <= il; ++i, ++ii) 
			res[ii] = -A*(u[ii-1]+u[ii+1])-B*(u[ii-dj]+u[ii+dj])+C*u[ii]+rhs[ii];
#endif  // XDIM == 2


#if XDIM == 3
	long i, i0, il, ii, j, j0, jl, jj, k, k0, kl, kk, 
		dj = n[0]+2, dk = dj*(n[1]+2), dj2 = 2*dj, dk2 = 2*dk;
	double A = 1./(h[0]*h[0]), B = 1./(h[1]*h[1]), C = 1./(h[2]*h[2]), D = 2.*(A+B+C);

#if BCX0+BCX1 > 3  // periodic
	i0 = 1, il = n[0]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 2;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+1; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]; j <= n[1]; ++j, ii += dj) res[ii] = 0.;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii-2];
#endif
#if BCY0+BCY1 > 3  // periodic
	j0 = 1, jl = n[1]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+1, j = jl*dj; i <= n[0]; ++i, ++ii) u[ii] = u[ii+j];
#endif
#if BCY0 == 0  // Diriclet
	j0 = 2;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY0 == 1  // von Neumann
	j0 = 1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dj2];
#endif
#if BCY1 == 0  // Diriclet
	jl = n[1]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj*n[1]+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCY1 == 1  // von Neumann
	jl = n[1];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj*n[1]+dj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dj2];
#endif
#if BCZ0+BCZ1 > 3  // periodic
	k0 = 1, kl = n[2]-1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+1, k = kl*dk; i <= n[0]; ++i, ++ii) u[ii] = u[ii+k];
#endif
#if BCZ0 == 0  // Diriclet
	k0 = 2;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk+jj+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCZ0 == 1  // von Neumann
	k0 = 1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dk2];
#endif
#if BCZ1 == 0  // Diriclet
	kl = n[2]-1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk*n[2]+jj+1; i <= n[0]; ++i, ++ii) res[ii] = 0.;
#endif
#if BCZ1 == 1  // von Neumann
	kl = n[2];
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk*n[2]+dk+jj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dk2];
#endif

  for (k = k0, kk = dk*k; k <= kl; ++k, kk += dk)  // interior points
		for (j = j0, jj = kk+dj*j0; j <= jl; ++j, jj += dj) 
			for (i = i0, ii = jj+i; i <= il; ++i, ++ii) 
				res[ii] = -A*(u[ii-1]+u[ii+1])-B*(u[ii-dj]+u[ii+dj])-C*(u[ii-dk]+u[ii+dk])+D*u[ii]+rhs[ii];
#endif  // XDIM == 3
}

#endif   // #ifdef MGLIN



void MultiGrid::addbc2rho(double *rho, const int *n, const double *h)
{
  // add effects of boundary conditions to rho
#if XDIM == 1
#if BCX0 == 1  // von Neumann
	rho[1] += 2.*rho[0]/h[0];
#endif
#if BCX1 == 1  // von Neumann
	rho[n[0]] -= 2.*rho[n[0]+1]/h[0];
#endif
#endif  // XDIM == 1

#if XDIM == 2
	
#if BCX0 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, ii = dj+1; j <= n[1]; ++j, ii += dj) 
		rho[ii] += 2.*rho[ii-1]/h[0];
#endif
#if BCX1 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, ii = dj+n[0]; j <= n[1]; ++j, ii += dj)
		rho[ii] -= 2.*rho[ii+1]/h[0];
#endif
#if BCY0 == 1  // von Neumann
	for (long i = 1, dj = n[0]+2; i <= n[0]; ++i) 
		rho[i+dj] += 2.*rho[i]/h[1];
#endif
#if BCY1 == 1  // von Neumann
	for (long i = 1, dj = n[0]+2, ii = dj*n[1]+1; i <= n[0]; ++i, ++ii)
		rho[ii] -= 2.*rho[ii+dj]/h[1];
#endif
#endif  // DIM == 2


#if XDIM == 3
#if BCX0 == 1  // von Neumann
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) 
			rho[ii+1] += 2.*rho[ii]/h[0];
#endif
#if BCX1 == 1  // von Neumann
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long j = 1, ii = kk+dj+n[0]; j <= n[1]; ++j, ii += dj)
			rho[ii] -= 2.*rho[ii+1]/h[0];
#endif
#if BCY0 == 1  // von Neumann 
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long i = 1, ii = kk+1; i <= n[0]; ++i, ++ii) 
			rho[ii+dj] += 2.*rho[ii]/h[1];
#endif
#if BCY1 == 1  // von Neumann
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long i = 1, ii = kk+dj*n[1]+1; i <= n[0]; ++i, ++ii)
			rho[ii] -= 2.*rho[ii+dj]/h[1];
#endif
#if BCZ0 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, jj = dj; j <= n[1]; ++j, jj += dj)
		for (long i = 1, ii = jj+1; i <= n[0]; ++i, ++ii) 
			rho[ii+dk] += 2.*rho[ii]/h[2];
#endif
#if BCZ1 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, dk = dj*(n[1]+2), jj = dj; j <= n[1]; ++j, jj += dj)
		for (long i = 1, ii = jj+dk*n[2]+1; i <= n[0]; ++i, ++ii)
			rho[ii] -= 2.*rho[ii+dk]/h[2];
#endif
#endif  // XDIM == 3
}

	
void MultiGrid::delbc2rho(double *rho, const int *n, const double *h)
{
  // remove added boundary conditions from rho

#if XDIM == 1
#if BCX0 == 1  // von Neumann
	rho[1] -= 2.*rho[0]/h[0];
#endif
#if BCX1 == 1  // von Neumann
	rho[n[0]] += 2.*rho[n[0]+1]/h[0];
#endif
#endif  // XDIM == 1


#if XDIM == 2
#if BCX0 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, ii = dj+1; j <= n[1]; ++j, ii += dj) 
		rho[ii] -= 2.*rho[ii-1]/h[0];
#endif
#if BCX1 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, ii = dj+n[0]; j <= n[1]; ++j, ii += dj)
		rho[ii] += 2.*rho[ii+1]/h[0];
#endif
#if BCY0 == 1  // von Neumann
	for (long i = 1; i <= n[0]; ++i) 
		rho[i+dj] -= 2.*rho[i]/h[1];
#endif
#if BCY1 == 1  // von Neumann
	for (long i = 1, dj = n[0]+2, ii = dj*n[1]+1; i <= n[0]; ++i, ++ii)
		rho[ii] += 2.*rho[ii+dj]/h[1];
#endif
#endif  // DIM == 2

#if XDIM == 3
#if BCX0 == 1  // von Neumann
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) 
			rho[ii+1] -= 2.*rho[ii]/h[0];
#endif
#if BCX1 == 1  // von Neumann
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long j = 1, ii = kk+dj+n[0]; j <= n[1]; ++j, ii += dj)
			rho[ii] += 2.*rho[ii+1]/h[0];
#endif
#if BCY0 == 1  // von Neumann 
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long i = 1, ii = kk+1; i <= n[0]; ++i, ++ii) 
			rho[ii+dj] -= 2.*rho[ii]/h[1];
#endif
#if BCY1 == 1  // von Neumann
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long i = 1, ii = kk+dj*n[1]+1; i <= n[0]; ++i, ++ii)
			rho[ii] += 2.*rho[ii+dj]/h[1];
#endif
#if BCZ0 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, dk = dj*(n[1]+2), jj = dj; j <= n[1]; ++j, jj += dj)
		for (long i = 1, ii = jj+1; i <= n[0]; ++i, ++ii) 
			rho[ii+dk] -= 2.*rho[ii]/h[2];
#endif
#if BCZ1 == 1  // von Neumann
	for (long j = 1, dj = n[0]+2, dk = dj*(n[1]+2), jj = dj; j <= n[1]; ++j, jj += dj)
		for (long i = 1, ii = jj+dk*n[2]+1; i <= n[0]; ++i, ++ii)
			rho[ii] += 2.*rho[ii+dk]/h[2];
#endif
#endif  // XDIM == 3
}

 
 
	
void MultiGrid::rstrct(double *uc, const double *u, const int *nc, const int *n)
// Half-weighting restriction and copy BC to rim points
{
#if XDIM == 1
	int ic, i;
	
  for (ic = 2, i = 3; ic < nc[0]; ++ic, i += 2)  // interior points
		uc[ic] = .5*u[i]+.25*(u[i-1]+u[i+1]);
		
	uc[1] = u[1], uc[nc[0]] = u[n[0]];  // boundary points
	uc[0] = u[0], uc[nc[0]+1] = u[n[0]+1];  
#endif  // XDIM == 1

#if XDIM == 2
  long ic, jc, i, j, iic, ii, jjc, jj, djc, dj;
	
	djc = nc[0]+2, dj = n[0]+2;
	for (jc = 2, j = 3; jc < nc[1]; ++jc, j += 2) { // interior points
		jjc = djc*jc, jj = dj*j;
		for (ic = 2, ii = jj+3, iic = jjc+ic; ic < n[0]; ++ic, ++iic, ii += 2)
			uc[iic] = .5*u[ii]+.125*(u[ii-1]+u[ii+1]+u[ii-dj]+u[ii+dj]);
	}
	
	for (jc = 1, j = 1; jc <= nc[1]; ++jc, j += 2) {  // x-boundary points
		ic = djc*jc, i = dj*j;
		uc[ic] = u[i], uc[ic+1] = u[i+1];
		ic += nc[0], i += n[0];
		uc[ic+1] = u[i+1], uc[ic] = u[i];
	}
	for (ic = 1, i = 1; ic <= nc[0]; ++ic, i += 2) {  // y-boundary points
		jc = ic, j = i;
		uc[jc] = u[j], uc[jc+djc] = u[j+dj];
		jc += djc*nc[1], j += dj*n[1];
		uc[jc+djc] = u[j+dj], uc[jc] = u[j];
	}
#endif  // XDIM == 2

#if XDIM == 3
  long ic, jc, kc, i, j, k, iic, ii, jjc, jj, kkc, kk, djc, dj, dkc, dk;
	double w = 1./12.;
	
	djc = nc[0]+2, dj = n[0]+2, dkc = djc*(nc[1]+2), dk = dj*(n[1]+2);
	for (kc = 2, k = 3; kc < nc[2]; ++kc, k += 2) { // interior points
		kkc = dkc*kc, kk = dk*k;
		for (jc = 2, j = 3; jc < nc[1]; ++jc, j += 2) { 
			jjc = kkc+djc*jc, jj = kk+dj*j;
			for (ic = 2, ii = jj+3, iic = jjc+ic; ic < n[0]; ++ic, ++iic, ii += 2)
				uc[iic] = .5*u[ii]+w*(u[ii-1]+u[ii+1]+u[ii-dj]+u[ii+dj]+u[ii-dk]+u[ii+dk]);
		}
	}
	
	for (kc = 1, k = 1; kc <= nc[2]; ++kc, k += 2) {  // x-boundary points
		kkc = dkc*kc, kk = dk*k;
		for (jc = 1, j = 1; jc <= nc[1]; ++jc, j += 2) {
			ic = kkc+djc*jc, i = kk+dj*j;
			uc[ic] = u[i], uc[ic+1] = u[i+1];
			ic += nc[0], i += n[0];
			uc[ic+1] = u[i+1], uc[ic] = u[i];
		}
	}
	for (kc = 1, k = 1; kc <= nc[2]; ++kc, k += 2) {  // y-boundary points
		kkc = dkc*kc, kk = dk*k;
		for (ic = 1, i = 1; ic <= nc[0]; ++ic, i += 2) {
			jc = kkc+ic, j = kk+i;
			uc[jc] = u[j], uc[jc+djc] = u[j+dj];
			jc += djc*nc[1], j += dj*n[1];
			uc[jc+djc] = u[j+dj], uc[jc] = u[j];
		}
	}
	for (jc = 1, j = 1; jc <= nc[1]; ++jc, j += 2) {  // z-boundary points
		jjc = djc*jc, jj = dj*j;
		for (ic = 1, i = 1; ic <= nc[0]; ++ic, i += 2) {
			kc = jjc+ic, k = jj+i;
			uc[kc] = u[k], uc[kc+dkc] = u[k+dk];
			kc += dkc*nc[2], k += dk*n[2];
			uc[kc+dkc] = u[k+dk], uc[kc] = u[k];
		}
	}
#endif  // XDIM == 3
}
			



void MultiGrid::interp(double *u, const double *uc, const double *rho, const int *n, const int *nc)
{
#if XDIM == 1
	int ic, i;
	
  for (ic = 1, i = 1; ic <= nc[0]; ++ic, i += 2)  // copy identical elements
		u[i] = uc[ic];
	for (i = 2; i <= n[0]; i += 2)  // interpolate intermediate elements
		u[i] = .5*(u[i-1]+u[i+1]);
#endif  // XDIM == 1

#if XDIM == 2
  long ic, i, jc, j, ii, jjc, jj, i0, il, j0, jl, djc = nc[0]+2, dj = n[0]+2, dj2 = 2*dj;
	
	for (jc = 1, j = 1; jc <= nc[1]; ++jc, j += 2) {	// copy identical elements
		jjc = djc*jc, jj = dj*j;
		for (ic = 1, ii = jj+1; ic <= nc[0]; ++ic, ii += 2)	u[ii] = uc[jjc+ic];
	}
	
	i0 = j0 = 1, il = n[0], jl = n[1];
#if BCX0 == 0  // Diriclet
	++i0;
	for (j = 2, ii = dj*j+1; j < n[1]; j += 2, ii += dj2) u[ii] = rho[ii-1];
#endif
#if BCX1 == 0  // Diriclet
	--il;
	for (j = 2, ii = dj*j+n[0]; j < n[1]; j += 2, ii += dj2) u[ii] = rho[ii+1];
#endif
#if BCY0 == 0  // Diriclet
	j0 += 2;
	for (i = 2, ii = dj+i; i < n[0]; i += 2, ii += 2) u[ii] = rho[ii-dj];
#endif
#if BCY1 == 0  // Diriclet
	jl -= 2;
	for (i = 2, ii = dj*n[1]+i; i < n[0]; i += 2, ii += 2) u[ii] = rho[ii+dj];
#endif

	for (j = j0, jj = dj*j; j <= jl; j += 2, jj += dj2)  // interpolate odd rows horizontally
		for (i = 2, ii = jj+i; i < n[0]; i += 2, ii += 2) u[ii] = .5*(u[ii-1]+u[ii+1]);

	for (j = 2, jj = dj*j; j < n[1]; j += 2, jj += dj2)  // interpolate even rows vertically
		for (i = i0, ii = jj+i; i <= il; ++i, ++ii)	u[ii] = .5*(u[ii-dj]+u[ii+dj]);
#endif  // XDIM == 2

#if XDIM == 3
  long ic, i, jc, j, kc, k, ii, jjc, jj, kkc, kk, i0, il, j0, jl, k0, kl, sw,
		djc = nc[0]+2, dj = n[0]+2, dkc = djc*(nc[1]+2), dk = dj*(n[1]+2), dj2 = 2*dj, dk2 = 2*dk;
	
	for (kc = 1, k = 1; kc <= nc[2]; ++kc, k += 2) {	// copy identical elements
		kkc = dkc*kc, kk = dk*k;
		for (jc = 1, j = 1; jc <= nc[1]; ++jc, j += 2) {
			jjc = kkc+djc*jc, jj = kk+dj*j;
			for (ic = 1, ii = jj+1; ic <= nc[0]; ++ic, ii += 2)	u[ii] = uc[jjc+ic];
		}
	}

	i0 = j0 = k0 = 1, il = n[0], jl = n[1], kl = n[2];
#if BCX0 == 0  // Diriclet
	++i0;
	for (k = 1, kk = dk*k, sw = 2; k <= n[2]; ++k, kk += dk, sw = 3-sw)
		for (j = sw, ii = kk+dj*j+1; j <= n[1]; j += sw, ii += dj*sw) u[ii] = rho[ii-1];
#endif
#if BCX1 == 0  // Diriclet
	--il;
	for (k = 1, kk = dk*k, sw = 2; k <= n[2]; ++k, kk += dk, sw = 3-sw)
		for (j = sw, ii = kk+dj*j+n[0]; j <= n[1]; j += sw, ii += dj*sw) u[ii] = rho[ii+1];
#endif
#if BCY0 == 0  // Diriclet
	j0 += 2;
	for (k = 1, kk = dk*k, sw = 2; k <= n[2]; ++k, kk += dk, sw = 3-sw)
		for (i = sw, ii = kk+dj+i; i <= n[0]; i += sw, ii += sw) u[ii] = rho[ii-dj];
#endif
#if BCY1 == 0  // Diriclet
	jl -= 2;
	for (k = 1, kk = dk*k, sw = 2; k <= n[2]; ++k, kk += dk, sw = 3-sw)
		for (i = sw, ii = kk+dj*n[1]+i; i <= n[0]; i += sw, ii += sw) u[ii] = rho[ii+dj];
#endif
#if BCZ0 == 0  // Diriclet
	k0 += 2;
	for (j = 1, jj = dj*j, sw = 2; j <= n[1]; ++j, jj += dj, sw = 3-sw)
		for (i = sw, ii = dk+jj+i; i <= n[0]; i += sw, ii += sw) u[ii] = rho[ii-dk];
#endif
#if BCZ1 == 0  // Diriclet
	kl -= 2;
	for (j = 1, jj = dj*j, sw = 2; j <= n[1]; ++j, jj += dj, sw = 3-sw)
		for (i = sw, ii = dk*n[2]+jj+i; i <= n[0]; i += sw, ii += sw) u[ii] = rho[ii+dk];
#endif

	for (k = k0, kk = dk*k; k <= kl; k += 2, kk += dk2) { // for odd planes do
		for (j = j0, jj = kk+dj*j; j <= jl; j += 2, jj += dj2)  // interpolate odd y-rows along x
			for (i = 2, ii = jj+i; i < n[0]; i += 2, ii += 2) u[ii] = .5*(u[ii-1]+u[ii+1]);

		for (j = 2, jj = kk+dj*j; j < n[1]; j += 2, jj += dj2)  // interpolate even y-rows along y
			for (i = i0, ii = jj+i; i <= il; ++i, ++ii)	u[ii] = .5*(u[ii-dj]+u[ii+dj]);
	}

  j0 = (1+j0)/2, jl = (n[1]+jl)/2;
	for (k = 2, kk = dk2; k < n[2]; k += 2, kk += dk2)  // for even planes interpolate between planes
		for (j = j0, jj = kk+dj*j; j <= jl; ++j, jj += dj)
			for (i = i0, ii = jj+i; i <= il; ++i, ++ii)	u[ii] = .5*(u[ii-dk]+u[ii+dk]);
#endif  // XDIM == 3
}



void MultiGrid::relax(double *u, const double *rhs, const int *n, const double *h)
{
	// red-black Gauss-Seidel iteration
#ifdef NZP
	double r, dfNZP;
//	double fNZP(const double u, double &dfNZP);
#endif

#if XDIM == 1
	int i, i0, il, ipass, isw;
	double A = 1./(h[0]*h[0]), B = -2.*A;
	
#if BCX0+BCX1 > 3  // periodic
	i0 = 0, il = n[0]-1;
	u[0] = u[il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 1;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 0;
	u[0] = u[2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	u[il+1] = u[il-1];
#endif

  for (ipass = 0, isw = 1; ipass < 2; ++ipass, isw = 3-isw) // red-black sweep
		for (i = isw+i0; i <= il; i += 2) {
#ifdef NZP
			r = A*(u[i-1]+u[i+1])+B*u[i]-rhs[i];
			u[i] -= (r+fNZP(u[i], dfNZP))/(B+dfNZP);
#else
			u[i] = -(A*(u[i-1]+u[i+1])-rhs[i])/B;
#endif
#if BCX0+BCX1 > 3  //periodic
			if (ipass == 0) u[n[0]] = u[1];
#endif
		}
#endif  // XDIM == 1


#if XDIM == 2
	long i, i0, il, ii, j, j0, jl, jj, ipass, isw, jsw, dj = n[0]+2, dj2 = 2*dj;
	double A = 1./(h[0]*h[0]), B = 1./(h[1]*h[1]), C = -2.*(A+B);
  
#if BCX0+BCX1 > 3  // periodic
	i0 = 0, il = n[0]-1;
	for (j = 1, ii = dj; j <= n[1]; j += 2, ii += dj2) u[ii] = u[ii+il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 1;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 0;
	for (j = 1, ii = dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	for (j = 1, ii = dj+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii-2];
#endif
#if BCY0+BCY1 > 3  // periodic
	j0 = 0, jl = n[1]-1;
	for (i = 1, ii = 1, j = jl*dj; i <= n[0]; i += 2, ii += 2) u[ii] = u[ii+j	];
#endif
#if BCY0 == 0  // Diriclet
	j0 = 1;
#endif
#if BCY0 == 1  // von Neumann
	j0 = 0;
	for (i = 1, ii = 1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dj2];
#endif
#if BCY1 == 0  // Diriclet
	jl = n[1]-1;
#endif
#if BCY1 == 1  // von Neumann
	jl = n[1];
	for (i = 1, ii = dj*n[1]+dj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dj2];
#endif

	for (ipass = 0, jsw = 1; ipass < 2; ++ipass, jsw = 3-jsw) { // red-black sweep
		isw = jsw;
		for (j = 1+j0, jj = dj*j; j <= jl; ++j, jj += dj, isw = 3-isw)
			for (i = isw+i0, ii = jj+i; i <= il; i += 2, ii += 2) {
#ifdef NZP
				r = A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*u[ii]-rhs[ii];
				u[ii] -= (r+fNZP(u[ii], dfNZP))/(C+dfNZP);
#else
				u[ii] = -(A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])-rhs[ii])/C;
#endif				
			}
#if BCX0+BCX1 > 3
		if (ipass == 0) {
			jj = n[0]-1;
			for (j = 1, ii = dj+i; j <= n[1]; j += 2, ii += dj2) u[ii+jj] = u[ii];
			for (j = 2, ii = dj2; j <= n[1]; j += 2, ii += dj2) u[ii] = u[ii+jj];
		}
#endif
#if BCY0+BCY1 > 3
		if (ipass == 0) {
			jj = (n[1]-1)*dj;
			for (i = 1, ii = dj+i; i <= n[0]; i += 2, ii += 2) u[ii+jj] = u[ii];
			for (i = 2; i <= n[0]; i += 2) u[i] = u[i+jj];
		}
#endif
	}
#if BCX0+BCX1 > 3
	jj = n[0]-1;
  for (j = 2, ii = dj2; j <= n[1]; j += 2, ii += dj2) u[ii+jj] = u[ii];
#endif
#if BCY0+BCY1 > 3
  jj = (n[1]-1)*dj;
	for (i = 2, ii = dj+i; i <= n[0]; i += 2, ii += 2) u[ii+jj] = u[ii];
#endif
#endif  // XDIM == 2


#if XDIM == 3
	long i, i0, il, ii, j, j0, jl, jj, k, k0, kl, kk, ipass, isw, jsw, ksw, 
		dj = n[0]+2, dj2 = 2.*dj, dk = dj*(n[1]+2), dk2 = 2*dk;
	double A = 1./(h[0]*h[0]), B = 1./(h[1]*h[1]), C = 1./(h[2]*h[2]), D = -2.*(A+B+C);
  
#if BCX0+BCX1 > 3  // periodic
  i0 = 0, il = n[0]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+il];
#endif
#if BCX0 == 0  // Diriclet
	i0 = 1;
#endif
#if BCX0 == 1  // von Neumann
	i0 = 0;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)		
		for (j = 1, ii = kk+dj; j <= n[1]; ++j, ii += dj) u[ii] = u[ii+2];
#endif
#if BCX1 == 0  // Diriclet
	il = n[0]-1;
#endif
#if BCX1 == 1  // von Neumann
	il = n[0];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]+1; j <= n[1]; ++j, ii += dj) u[ii] = u[ii-2];
#endif
#if BCY0+BCY1 > 3  // periodic
	j0 = 0, jl = n[1]-1;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+1, j = jl*dj; i <= n[0]; ++i, ++ii) u[ii] = u[ii+j];
#endif
#if BCY0 == 0  // Diriclet
	j0 = 1;
#endif
#if BCY0 == 1  // von Neumann
	j0 = 0;
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dj2];
#endif
#if BCY1 == 0  // Diriclet
	jl = n[1]-1;
#endif
#if BCY1 == 1  // von Neumann
	jl = n[1];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj*n[1]+dj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dj2];
#endif
#if BCZ0+BCZ1 > 3  // periodic
	k0 = 0, kl = n[2]-1;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+1, k = kl*dk; i <= n[0]; ++i, ++ii) u[ii] = u[ii+k];
#endif
#if BCZ0 == 0  // Diriclet
	k0 = 1;
#endif
#if BCZ0 == 1  // von Neumann
	k0 = 0;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii+dk2];
#endif
#if BCZ1 == 0  // Diriclet
	kl = n[2]-1;
#endif
#if BCZ1 == 1  // von Neumann
	kl = n[2];
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk*n[2]+dk+jj+1; i <= n[0]; ++i, ++ii) u[ii] = u[ii-dk2];
#endif

	for (ipass = 0, ksw = 1; ipass < 2; ++ipass, ksw = 3-ksw) { // red-black sweep
		jsw = ksw;
		for (k = 1+k0, kk = dk*k; k <= kl; ++k, kk += dk, jsw = 3-jsw) {
			isw = jsw;
			for (j = 1+j0, jj = kk+dj*j; j <= jl; ++j, jj += dj, isw = 3-isw)
				for (i = isw+i0, ii = jj+i; i <= il; i += 2, ii += 2) {
#ifdef NZP
					r = A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*(u[ii-dk]+u[ii+dk])
						+D*u[ii]-rhs[ii];
					u[ii] -= (r+fNZP(u[ii], dfNZP))/(D+dfNZP);
#else
					u[ii] = -(A*(u[ii-1]+u[ii+1])+B*(u[ii-dj]+u[ii+dj])+C*(u[ii-dk]+u[ii+dk])
						-rhs[ii])/D;
#endif
				}
		}
	}
#endif  // XDIM == 3
}




void MultiGrid::slvsml(double *u, const double *rhs, const int *n, const double *h)
{
	
#if XDIM == 1
	double A = 1./(h[0]*h[0]), B = 2.*A;
#if BCX0 == 0
	u[1] = rhs[0];
#if BCX1 == 0
	u[3] = rhs[4];
#endif
#if BCX1 == 1
	u[3] = rhs[0]-(rhs[2]+rhs[3])/A;
#endif
#endif
#if BCX0 == 1
#if BCX1 == 0
	u[1] = rhs[4]-(rhs[1]+rhs[2])/A;
	u[3] = rhs[4];
#endif
#endif
	u[2] = (A*(u[1]+u[3])-rhs[2])/B;
#ifdef NZP
  for (int i = 0; i < 10; ++i) relax(u, rhs, n, h);
#endif
#endif  // XDIM == 1


#if XDIM == 2
#if BCX0 == 0
	for (long j = 1, dj = n[0]+2, ii = dj+1; j <= n[1]; ++j, ii += dj) 
		u[ii] = rhs[ii-1];
#endif
#if BCX1 == 0
	for (long j = 1, dj = n[0]+2, ii = dj+n[0]; j <= n[1]; ++j, ii += dj) 
		u[ii] = rhs[ii+1];
#endif
#if BCY0 == 0
	for (long i = 1, dj = n[0]+2, ii = dj+1; i <= n[0]; ++i, ++ii)
		u[ii] = rhs[ii-dj];
#endif
#if BCY1 == 0
	for (long i = 1, dj = n[0]+2, ii = dj*n[1]+1; i <= n[0]; ++i, ++ii)
		u[ii] = rhs[ii+dj];
#endif

  for (int i = 0; i < 10; ++i) relax(u, rhs, n, h);
#endif  // XDIM == 2


#if XDIM == 3
#if BCX0 == 0
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long j = 1, ii = kk+dj+1; j <= n[1]; ++j, ii += dj) u[ii] = rhs[ii-1];
#endif
#if BCX1 == 0
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long j = 1, ii = kk+dj+n[0]; j <= n[1]; ++j, ii += dj) u[ii] = rhs[ii+1];
#endif
#if BCY0 == 0
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long i = 1, ii = kk+dj+1; i <= n[0]; ++i, ++ii) u[ii] = rhs[ii-dj];
#endif
#if BCY1 == 0
	for (long k = 1, dj = n[0]+2, dk = dj*(n[1]+2), kk = dk; k <= n[2]; ++k, kk += dk)
		for (long i = 1, ii = kk+dj*n[1]+1; i <= n[0]; ++i, ++ii) u[ii] = rhs[ii+dj];
#endif
#if BCZ0 == 0
	for (long j = 1, dj = n[0]+2, dk = dj*(n[1]+2), jj = dj; j <= n[1]; ++j, jj += dj)
		for (long i = 1, ii = dk+jj+1; i <= n[0]; ++i, ++ii) u[ii] = rhs[ii-dk];
#endif
#if BCZ1 == 0
	for (long j = 1, dj = n[0]+2, dk = dj*(n[1]+2), jj = dj; j <= n[1]; ++j, jj += dj)
		for (long i = 1, ii = dk*n[2]+jj+1; i <= n[0]; ++i, ++ii) u[ii] = rhs[ii+dk];
#endif

  for (int i = 0; i < 10; ++i) relax(u, rhs, n, h);
#endif  // XDIM == 3
}




void MultiGrid::fill0(double *u, const int *n)
{
#if XDIM == 1
  long i;
	
  for (i = 1; i <= n[0]; ++i) u[i] = 0.;
#endif  // XDIM == 1


#if XDIM == 2
  long i, j, ii, jj, dj = n[0]+2;
	
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) u[ii] = 0.;
#endif  // XDIM == 2


#if XDIM == 3
  long i, j, k, ii, jj, kk, dj = n[0]+2, dk = dj*(n[1]+2);
	
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk) 
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) u[ii] = 0.;
#endif
}




void MultiGrid::TestMG()
{
  // Test of multigrid routines
	
	int i, n[XDIM];
  double x[XDIM], h[XDIM], Uf, dUf[XDIM], d2Uf, *U, *d2U, *rho, *u, D;
	
	long gmax = nnp[ngrid], dim;	
	for (i = 0; i < XDIM; ++i) { n[i] = nn[ngrid][0], h[i] = hh[ngrid][i]; }
	U = new double[gmax];
	d2U = new double[gmax];
	rho = new double[gmax];
	u = new double[gmax];
	
#if XDIM == 1
  cout<<"BC "<<BCX0<<" "<<BCX1<<endl;
	x[0] = 0.;
	Ufun(Uf, dUf, d2Uf, x);
	U[1] = Uf, d2U[1] = rho[1] = d2Uf+fNZP(Uf,D); 
#if BCX0 == 0  // Diriclet
  d2U[0] = rho[0] = Uf;
#endif
#if BCX0 == 1  // von Neumann
	d2U[0] = rho[0] = dUf[0];
#endif

	for (i = 2; i <= n[0]; ++i) { // interior points
		x[0] += h[0];
		Ufun(Uf, dUf, d2Uf, x);
		U[i] = Uf, d2U[i] = rho[i] = d2Uf+fNZP(Uf,D);
	}

#if BCX1 == 0  // Diriclet
	d2U[n[0]+1] = rho[n[0]+1] = Uf;
#endif
#if BCX1 == 1  // von Neumann 
	d2U[n[0]+1] = rho[n[0]+1] = dUf[0];
#endif
	
#endif  // XDIM == 1

#if XDIM == 2
	long j, ii, jj, dj = n[0]+2;
	
  cout<<"BC "<<BCX0<<" "<<BCX1<<" "<<BCY0<<" "<<BCY1<<endl;
	x[0] = 0.;
	for (j = 1, ii = dj+1, x[1] = 0.; j <= n[1]; ++j, ii += dj, x[1] += h[1]) {
		Ufun(Uf, dUf, d2Uf, x);
		U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCX0 == 0  // Diriclet
		d2U[ii-1] = rho[ii-1] = Uf;
#endif
#if BCX0 == 1  // von Neumann
		d2U[ii-1] = rho[ii-1] = dUf[0];
#endif
	}

	x[0] = L[0];
	for (j = 1, ii = n[0]+dj, x[1] = 0.; j <= n[1]; ++j, ii += dj, x[1] += h[1]) {
		Ufun(Uf, dUf, d2Uf, x);
		U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCX1 == 0  // Diriclet
		d2U[ii+1] = rho[ii+1] = Uf;
#endif
#if BCX1 == 1  // von Neumann 
		d2U[ii+1] = rho[ii+1] = dUf[0];
#endif
	}

	x[1] = 0.;
	for (i = 1, ii = dj+1, x[0] = 0.; i <= n[0]; ++i, ++ii, x[0] += h[0]) {
		Ufun(Uf, dUf, d2Uf, x);
		U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCY0 == 0  // Diriclet
		d2U[ii-dj] = rho[ii-dj] = Uf;
#endif
#if BCY0 == 1  // von Neumann
		d2U[ii-dj] = rho[ii-dj] = dUf[1];
#endif
	}

	x[1] = L[1];
	for (i = 1, ii = dj*n[1]+1, x[0] = 0.; i <= n[0]; ++i, ++ii, x[0] += h[0]) {
		Ufun(Uf, dUf, d2Uf, x);
		U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCY1 == 0  // Diriclet
		d2U[ii+dj] = rho[ii+dj] = Uf;
#endif
#if BCY1 == 1  // von Neumann 
		d2U[ii+dj] = rho[ii+dj] = dUf[1];
#endif
	}

	// interior points
	for (j = 2, jj = dj*j, x[1] = h[1]; j < n[1]; ++j, jj += dj, x[1] += h[1])
		for (i = 2, ii = jj+i, x[0] = h[0]; i < n[0]; ++i, ++ii, x[0] += h[0]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
		}
#endif  // XDIM == 2

#if XDIM == 3
	long j, k, ii, jj, kk, dj = n[0]+2, dk = dj*(n[1]+2);
	
  cout<<"BC "<<BCX0<<" "<<BCX1<<" "<<BCY0<<" "<<BCY1<<" "<<BCZ0<<" "<<BCZ1<<endl;
	x[0] = 0.;
	for (k = 1, kk = dk*k, x[2] = 0.; k <= n[2]; ++k, kk += dk, x[2] += h[2])
		for (j = 1, ii = kk+dj*j+1, x[1] = 0.; j <= n[1]; ++j, ii += dj, x[1] += h[1]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCX0 == 0  // Diriclet
			d2U[ii-1] = rho[ii-1] = Uf;
#endif
#if BCX0 == 1  // von Neumann
			d2U[ii-1] = rho[ii-1] = dUf[0];
#endif
		}

	x[0] = L[0];
	for (k = 1, kk = dk*k, x[2] = 0.; k <= n[2]; ++k, kk += dk, x[2] += h[2])
		for (j = 1, ii = kk+dj*j+n[0], x[1] = 0.; j <= n[1]; ++j, ii += dj, x[1] += h[1]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCX1 == 0  // Diriclet
			d2U[ii+1] = rho[ii+1] = Uf;
#endif
#if BCX1 == 1  // von Neumann 
			d2U[ii+1] = rho[ii+1] = dUf[0];
#endif
		}

	x[1] = 0.;
	for (k = 1, kk = dk*k, x[2] = 0.; k <= n[2]; ++k, kk += dk, x[2] += h[2])
		for (i = 1, ii = kk+dj+i, x[0] = 0.; i <= n[0]; ++i, ++ii, x[0] += h[0]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCY0 == 0  // Diriclet
			d2U[ii-dj] = rho[ii-dj] = Uf;
#endif
#if BCY0 == 1  // von Neumann
			d2U[ii-dj] = rho[ii-dj] = dUf[1];
#endif
		}
	
	x[1] = L[1];
	for (k = 1, kk = dk*k, x[2] = 0.; k <= n[2]; ++k, kk += dk, x[2] += h[2])
		for (i = 1, ii = kk+dj*n[1]+i, x[0] = 0.; i <= n[0]; ++i, ++ii, x[0] += h[0]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCY1 == 0  // Diriclet
			d2U[ii+dj] = rho[ii+dj] = Uf;
#endif
#if BCY1 == 1  // von Neumann 
			d2U[ii+dj] = rho[ii+dj] = dUf[1];
#endif
		}

	x[2] = 0.;
	for (j = 1, jj = dk+dj*j, x[1] = 0.; j <= n[1]; ++j, jj += dj, x[1] += h[1])
		for (i = 1, ii = jj+i, x[0] = 0.; i <= n[0]; ++i, ++ii, x[0] += h[0]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCZ0 == 0  // Diriclet
			d2U[ii-dk] = rho[ii-dk] = Uf;
#endif
#if BCZ0 == 1  // von Neumann
			d2U[ii-dk] = rho[ii-dk] = dUf[2];
#endif
		}
	
	x[2] = L[2];
	for (j = 1, jj = dk*n[2]+dj*j, x[1] = 0.; j <= n[1]; ++j, jj += dj, x[1] += h[1])
		for (i = 1, ii = jj+i, x[0] = 0.; i <= n[0]; ++i, ++ii, x[0] += h[0]) {
			Ufun(Uf, dUf, d2Uf, x);
			U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
#if BCZ1 == 0  // Diriclet
			d2U[ii+dk] = rho[ii+dk] = Uf;
#endif
#if BCZ1 == 1  // von Neumann 
			d2U[ii+dk] = rho[ii+dk] = dUf[2];
#endif
		}

  // interior points
	for (k = 2, kk = dk*k, x[2] = h[2]; k < n[2]; ++k, kk += dk, x[2] += h[2])
		for (j = 2, jj = kk+dj*j, x[1] = h[1]; j < n[1]; ++j, jj += dj, x[1] += h[1])
			for (i = 2, ii = jj+i, x[0] = h[0]; i < n[0]; ++i, ++ii, x[0] += h[0]) {
				Ufun(Uf, dUf, d2Uf, x);
				U[ii] = Uf, d2U[ii] = rho[ii] = d2Uf+fNZP(Uf,D);
			}
#endif  // XDIM == 3

#ifdef MGPRNT
  prnt("U ", n, h, U);
#endif

#ifdef MGFAS
	mgfas(u, rho, n, h);
#endif	
#ifdef MGLIN
	mglin(u, rho, n, h);
#endif

  // shift solution for all von Neumann BC
#if XDIM == 1
#if BCX0*BCX1 == 1
  double shift = 0.;
	for (i = 1; i <= n[0]; ++i) shift += u[i]-U[i];
	shift /= n[0];
	for (i = 1; i <= n[0]; ++i) u[i] -= shift;
#endif
#endif

#if XDIM == 2
#if BCX0*BCX1*BCY0*BCY1 == 1
  double shift = 0.;
	for (j = 1, jj = dj*j; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) shift += u[ii]-U[ii];
	shift /= (n[0]*n[1]);
	for (j = 1, jj = dj*j; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) u[ii] -= shift;
#endif
#endif
  

#if XDIM == 3
#if BCX0*BCX1*BCY0*BCY1*BCZ0*BCZ1 == 1
  double shift = 0.;
	for (k = 1, kk = dk*k; k <= n[2]; ++k, kk += dk)
		for (j = 1, jj = kk+dj*j; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) shift += u[ii]-U[ii];
	shift /= (n[0]*n[1]);
	for (k = 1, kk = dk*k; k <= n[2]; ++k, kk += dk)
		for (j = 1, jj = kk+dj*j; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) u[ii] -= shift;
#endif
#endif
  
  dump(n, h, u, U);
	
#if XDIM == 1
	D = 0., dim = n[0];
	for (i = 1; i <= n[0]; ++i) 
		D += (U[i]-u[i])*(U[i]-u[i]);
	cout<<"D = "<<sqrt(D)/dim<<endl;
#endif
#if XDIM == 2
	D = 0., dim = n[0]*n[1];
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) 
			D += (U[ii]-u[ii])*(U[ii]-u[ii]);
	cout<<"D = "<<sqrt(D)/dim<<endl;
#endif
#if XDIM == 3
	D = 0., dim = n[0]*n[1]*n[2];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) 
				D += (U[ii]-u[ii])*(U[ii]-u[ii]);
	cout<<"D = "<<sqrt(D)/dim<<endl;
#endif
	
	addbc2rho(d2U, n, h);
	fill0(u, n);
#if XDIM == 1
#if BCX0 == 0
  u[1] = U[1];
#endif
#if BCX1 == 0
  u[n[0]] = U[n[0]];
#endif
#endif
#if XDIM == 2
#if BCX0 == 0
	for (j = 1, ii = dj+1; j <= n[1]; ++j, ii += dj)  u[ii] = U[ii];
#endif
#if BCX1 == 0
	for (j = 1, ii = dj+n[0]; j <= n[1]; ++j, ii += dj)  u[ii] = U[ii];
#endif
#if BCY0 == 0
	for (i = 1, ii = dj+i; i <= n[0]; ++i, ++ii)  u[ii] = U[ii];
#endif
#if BCY1 == 0
	for (i = 1, ii = dj+i; i <= n[0]; ++i, ++ii)  u[ii] = U[ii];
#endif
#endif
#if XDIM == 3
#if BCX0 == 0
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+1; j <= n[1]; ++j, ii += dj)  u[ii] = U[ii];
#endif
#if BCX1 == 0
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, ii = kk+dj+n[0]; j <= n[1]; ++j, ii += dj)  u[ii] = U[ii];
#endif
#if BCY0 == 0
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (i = 1, ii = kk+dj+i; i <= n[0]; ++i, ++ii)  u[ii] = U[ii];
#endif
#if BCY1 == 0
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
	for (i = 1, ii = kk+dj+i; i <= n[0]; ++i, ++ii)  u[ii] = U[ii];
#endif
#if BCZ0 == 0
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk+jj+i; i <= n[0]; ++i, ++ii)  u[ii] = U[ii];
#endif
#if BCZ1 == 0
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = dk+jj+i; i <= n[0]; ++i, ++ii)  u[ii] = U[ii];
#endif
#endif


	for (i = 0; i < 10; ++i) relax(u, d2U, n, h);
#ifdef MGPRNT
	prnt("u relax ONLY ", n, h, u);
#endif

#if XDIM == 1
	D = 0., dim = n[0];
	for (i = 1; i <= n[0]; ++i) 
		D += (U[i]-u[i])*(U[i]-u[i]);
	cout<<"D = "<<sqrt(D)/dim<<endl;
#endif
#if XDIM == 2
	D = 0., dim = n[0]*n[1];
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj)
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) 
			D += (U[ii]-u[ii])*(U[ii]-u[ii]);
	cout<<"D = "<<sqrt(D)/dim<<endl;
#endif
#if XDIM == 3
	D = 0., dim = n[0]*n[1]*n[2];
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj)
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) 
				D += (U[ii]-u[ii])*(U[ii]-u[ii]);
	cout<<"D = "<<sqrt(D)/dim<<endl;
#endif
	
	
	delete [] U;
	delete [] d2U;
	delete [] u;
}




void Ufun(double &Uf, double *dUf, double &d2Uf, double *x)
{
	double K[3] = {2.*M_PI, 2.*M_PI, 2.*M_PI};

#if XDIM == 1
  Uf = sin(K[0]*x[0]);
	dUf[0] = K[0]*cos(K[0]*x[0]);
	d2Uf = -K[0]*K[0]*Uf;
#endif

#if XDIM == 2
  Uf = sin(K[0]*x[0])*sin(K[1]*x[1]);
	dUf[0] = K[0]*cos(K[0]*x[0])*sin(K[1]*x[1]);
	dUf[1] = K[1]*sin(K[0]*x[0])*cos(K[1]*x[1]);
	d2Uf = -(K[0]*K[0]+K[1]*K[1])*Uf;
#endif

#if XDIM == 3
  Uf = sin(K[0]*x[0])*sin(K[1]*x[1])*sin(K[2]*x[2]);
	dUf[0] = K[0]*cos(K[0]*x[0])*sin(K[1]*x[1])*sin(K[2]*x[2]);
	dUf[1] = K[1]*sin(K[0]*x[0])*cos(K[1]*x[1])*sin(K[2]*x[2]);
	dUf[2] = K[2]*sin(K[0]*x[0])*sin(K[1]*x[1])*cos(K[2]*x[2]);
	d2Uf = -(K[0]*K[0]+K[1]*K[1]+K[2]*K[2])*Uf;
#endif
}




double fNZP(const double u, double &dfNZP)
{
	dfNZP = 0.; //2.*u;  //0.;
	return 0.;  //u*u;  // 0.;
}



void prnt(const char txt[], const int *n, const double *h, const double *u)
{
#if XDIM == 1
	int i, istep;
	
	cout<<txt<<n[0]<<" "<<h[0]<<endl;
	istep = n[0]/8; if (istep == 0) istep = 1;
//	cout<<u[0]<<" "; 
	for (i = 1; i <= n[0]; i += istep) cout<<u[i]<<" "; 
//	cout<<u[n[0]+1];
	cout<<endl;
#endif

#if XDIM == 2
	long i, j, ii, jj, istep, jstep, dj = n[0]+2;
	cout<<txt; for (i = 0; i < XDIM; ++i) cout<<n[0]<<" ";
	for (i = 0; i < XDIM; ++i) cout<<h[0]<<" "; cout<<endl;
	istep = n[0]/8; if (istep == 0) istep = 1;
	jstep = n[1]/8; if (jstep == 0) jstep = 1;
	
//	cout<<u[dj*(n[1]+1)]<<" ";
//	for (i = 1, ii = dj*(n[1]+1); i <= n[0]; i += istep, ii += istep) cout<<u[ii]<<" "; 
//	cout<<u[dj*(n[1]+1)+n[0]+1]<<endl;
	for (j = n[1], jj = dj*j; j >= 1; j -= jstep, jj -= dj*jstep) {
//		cout<<u[jj]<<" ";
		for (i = 1, ii = jj+i; i <= n[0]; i += istep, ii += istep) cout<<u[ii]<<" ";
//		cout<<u[jj+n[0]+1];
		cout<<endl;
	}
//	cout<<u[0]<<" ";
//	for (i = 1; i <= n[0]; i += istep) cout<<u[i]<<" ";
//	cout<<u[n[0]+1]<<endl;
#endif

#if XDIM == 3
	long i, j, k, ii, jj, kk, istep, jstep, kstep, dj = n[0]+2, dk = dj*(n[1]+2);
	cout<<txt; for (i = 0; i < XDIM; ++i) cout<<n[0]<<" ";
	for (i = 0; i < XDIM; ++i) cout<<h[0]<<" "; cout<<endl;
	istep = n[0]/8; if (istep == 0) istep = 1;
	jstep = n[1]/8; if (jstep == 0) jstep = 1;
	kstep = n[2]/8; if (kstep == 0) kstep = 1;
	
//	for (k = 1, kk = dk; k <= n[2]; k += kstep, kk += dk*kstep) {
	for (k = n[2]/2+1, kk = dk*k; k <= n[2]/2+1; k += kstep, kk += dk*kstep) {
		cout<<"k = "<<k<<endl;
		for (j = 1, jj = kk+dj; j <= n[1]; j += jstep, jj += dj*jstep) {
			for (i = 1, ii = jj+i; i <= n[0]; i += istep, ii += istep) cout<<u[ii]<<" "; cout<<endl;
		}
		cout<<endl;
	}
#endif
}


void dump(const int *n, const double *h, const double *u, const double *U)
{
  long i;
	string utfil = "/Users/jtrulsen/pic/picMC/MG.dat";
	ofstream ut(utfil.c_str(), ios::out);

	ut<<XDIM<<endl;
	for (i = 0; i < XDIM; ++i) ut<<n[i]<<" "; ut<<endl;
	for (i = 0; i < XDIM; ++i) ut<<h[i]<<" "; ut<<endl;
	
#if XDIM == 1
  for (i = 1; i <= n[0]; ++i) ut<<u[i]<<" "; ut<<endl;
  for (i = 1; i <= n[0]; ++i) ut<<U[i]<<" "; ut<<endl;
#endif

#if XDIM == 2
  long j, ii, jj, dj = n[0]+2;
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj) {
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) ut<<u[ii]<<" "; ut<<endl;
	}
	for (j = 1, jj = dj; j <= n[1]; ++j, jj += dj) {
		for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) ut<<U[ii]<<" "; ut<<endl;
	}
#endif

#if XDIM == 3
  long j, k, ii, jj, kk, dj = n[0]+2, dk = dj*(n[1]+2);
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj) {
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) ut<<u[ii]<<" "; ut<<endl;
		}
	for (k = 1, kk = dk; k <= n[2]; ++k, kk += dk)
		for (j = 1, jj = kk+dj; j <= n[1]; ++j, jj += dj) {
			for (i = 1, ii = jj+i; i <= n[0]; ++i, ++ii) ut<<U[ii]<<" "; ut<<endl;
		}
#endif

	ut.close();
}

