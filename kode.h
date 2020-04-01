	
	/* **************************
	 * *     HEADER START       *
	 * **************************
	 * */
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
//#include <time.h>
//#include <sys/time.h>
#include <ctime>
using namespace std;


#define pi 3.14159265359
/*
  //   PRIMEROOT RANDOM NUMBER, WOJCIECH
#define BUCKETSIZE 1000
double primerootbucket[BUCKETSIZE];
double primerootno;
#ifndef DRAND48
inline double primeroot(void){
	double r;
	double sqrt_two=sqrt(2.0);
	int draw=(int)floor((drand48()*BUCKETSIZE));
	r=primerootbucket[draw];
	primerootno=primerootno+sqrt_two-floor((primerootno+sqrt_two));
	primerootbucket[draw]=primerootno;
	return r;
}
#else
double primeroot(void){return drand48();}
#endif
void init_primeroot(double seed){
	int i;
	double x,sqrt_two;
	sqrt_two=sqrt(2);
	srand48(0);
	x=seed;
	for(i=0; i<BUCKETSIZE; i++){
		x=x+sqrt_two-(int)(x+sqrt_two);
		primerootbucket[i]=x;
	}
	primerootno=x;
}
*/

double sign(double x){return x/abs(x);}
void lowpassfilter(double f[], int N){
	for(int i=1; i<N-1; ++i) f[i]=0.25*f[i-1]+0.5*f[i]+0.25*f[i+1];
}


void ErrorMessage(string text){
	cout << text << endl;
	exit(1);
}

/*
	// Box-Mueller algorithm for normal distributed numbers 
	// N(mu = 0, sigma = 1)
double BoxMuller(){
	static int FIRST = 1;
	static double NEXT;
	double v1, v2, sqr, fac;
	
	if (FIRST) {
		do {
			v1 = 2.*drand48() - 1.;
			v2 = 2.*drand48() - 1.;
			sqr = v1*v1 + v2*v2; 
		} while (sqr >= 1. || sqr == 0.);

		fac = sqrt(-2.*log(sqr)/sqr);
		NEXT = v2*fac;
		FIRST = 0;
		return v1*fac;
	} else {
		FIRST = 1;
		return NEXT;
	}
}
*/








	/* ***********************************************
	 *            ERROR FUNCTION WITH SUBROUTINES
 	 * ************************************************/
	/*	Error function
	 * Source: Numerical Recepies in C.
	 * Changed: float to double declarations
	 * */

#include <math.h>
//Returns the errorfunction 
double erf(double x){
	double t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return  x >= 0.0 ? 1.0-ans : ans-1.0;
}


//Returns the complementary error function
double erfc(double x){
	double t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return  x >= 0.0 ? ans : 2.0-ans;
}



#include <math.h>
#define ITMAX 100 //Maximum allowed number of iterations.
#define EPS 3.0e-7 //Relative accuracy.
#define FPMIN 1.0e-30 //Number near the smallest representable

/* Numerical Recipes standard error handler */
void nrerror(char error_text[]){
        fprintf(stderr,"Numerical Recipes run-time error...\n");
        fprintf(stderr,"%s\n",error_text);
        fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}



//Metropolis algortihms. 
//Written by Jan K. Trulsen, UiO in 2007 as an aid to my masters thesis. 
//Source: Trulsen's lecturenotes on statistics in numerical analysis
//Some adjustments made by me to fit my DL simulations (master thesis).
//							Vegard L. Rekaa
void Metropolis(double n[], int N, double *x, double *f, 
				double delta, double a, double b){
	// Metropolis algorithm for a bounded domain (a, b)
  	// pdf(x) : desired (unnormalized) probability density function
	// x      : input previous draw, output present draw
	// f      : f = pdf(x) at both input and output
	// delta  : random walk step length
	// Initialize new series by inputing f = 0 and x near max of pdf()
	
	double xt, ft;
	int i;
	double w;

	xt = *x + delta*(2.*drand48() - 1.);
	if (xt < a) xt = 2.*a-xt;  // mirroring at lower boundary
	if (xt > b) xt = 2.*b-xt;  // mirroring at upper boundary
	
	w=xt*N;
	i=int(floor(w));
	w-=i;

	ft = n[i]*(1.-w) + n[i+1]*w;
	if ( ft > *f || ft >= *f*drand48()) {
		*x = xt, *f = ft;
	}
}

void Metropolis(double (*pdf)(double,double,double), double u, double s,
		double *x, double *f, double delta, double a, double b){
	// Metropolis algorithm for a bounded domain (a, b)
  	// pdf(x) : desired (unnormalized) probability density function
	// x      : input previous draw, output present draw
	// f      : f = pdf(x) at both input and output
	// delta  : random walk step length
	// Initialize new series by inputing f = 0 and x near max of pdf()
	
	double xt, ft;

	xt = *x + delta*(2.*drand48() - 1.);
	if (xt < a) xt = 2.*a-xt;  // mirroring at lower boundary
	if (xt > b) xt = 2.*b-xt;  // mirroring at upper boundary
		
	ft = pdf(xt, u, s);
	if ( ft > *f || ft >= *f*drand48()) {
		*x = xt, *f = ft;
	}
}

void Metropolis(double (*pdf)(double, double), double u, double *x, 
		double *f, double delta, double a){
	// Metropolis algorithm for a bounded domain (a, b)
  	// pdf(x) : desired (unnormalized) probability density function
	// x      : input previous draw, output present draw
	// f      : f = pdf(x) at both input and output
	// delta  : random walk step length
	// Initialize new series by inputing f = 0 and x near max of pdf()
	
	double xt, ft;

	xt = *x + delta*(2.*drand48() - 1.);
	if (xt < a) xt = 2.*a-xt;  // mirroring at lower boundary
		
	ft = pdf(xt, u);
	if ( ft > *f || ft >= *f*drand48()) {
		*x = xt, *f = ft;
	}
}
