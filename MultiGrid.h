/*
 *  MultiGrid.h
 *  pic3D
 *
 *  Created by Jan Trulsen on 05/01/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MULTIGRID_H
#define MULTIGRID_H

//#include "Definitions.h"

/*
#define MAXCYC 1
#define MAXIT0 100
#define MAXDIT 0
#define EPSOLVE1 1.e-1
#define EPSOLVE2 1.e-8
*/

// Maximum multigrid depth
#define NGMAX 15

// Boundary type: 0 Dirichlet, 1 von Neumann, 2 periodic
//#define BCX0 1
//#define BCX1 1
//#define BCY0 0
//#define BCY1 1
//#define BCZ0 0
//#define BCZ1 1

#define MGLIN
//#define MGFAS
//#define NZP

class MultiGrid
{

 public:
	int ngrid, ncycle, npre, npost, nn[NGMAX][XDIM];
	long nnp[NGMAX];
	double hh[NGMAX][XDIM], L[XDIM], alpha, trerr, res;
	double *irho[NGMAX-1], *iu[NGMAX-1], *itmp[NGMAX];
#ifdef MGFAS
	double *itau[NGMAX];
#endif
	
 public:
  MultiGrid(const int* n, const double* h);
	~MultiGrid();
	
#ifdef MGFAS
  void mgfas(double *u, double *rho, const int *n, const double *h);
#endif
#ifdef MGLIN
  void mglin(double *u, double *rho, const int *n, const double *h);
#endif
	void relax(double *u, const double *rhs, const int *n, const double *h);
	void addbc2rho(double *u, const int *no, const double *h);
	void delbc2rho(double *u, const int *no, const double *h);

  void TestMG();

 private:
#ifdef MGFAS
	void lop(double *res, double *u, const int *n, const double *h);
	void matsub(double *c, const double *a, const double *b, const int *n);
	void matadd(double *c, const double *a, const double *b, const int *n);
	double anorm(const double *a, const int *n);
#endif
#ifdef MGLIN
	void addint(double *u, const double *uc, double *res, const double *rho, const int *n, const int *nc);
	void resid(double *res, double *u, const double *rhs, const int *n, const double *h);
#endif
	void rstrct(double *uc, const double *u, const int *nc, const int *n);
	void interp(double *u, const double *uc, const double *rho, const int *n, const int *nc);
	void slvsml(double *u, const double *rhs, const int *n, const double *h);
	void fill0(double *u, const int *n);
};


	void Ufun(double &Uf, double *dUf, double &d2Uf, double *x);
	double fNZP(const double u, double &dfNZP);
	void prnt(const char txt[], const int *n, const double *h, const double *u);
	void dump(const int *n, const double *h, const double *u, const double *U);
#endif
