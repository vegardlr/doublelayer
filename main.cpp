	/*		DOUBLE LAYER PLASMA SIMULATION
	 *			Vegard L. Rekaa 
	 *			UiO, 2007/2008
	 *
	 *	This file contains globally decleared arrays for position
	 *	and velocities, and parameters used in most functions. 
	 *	I am not sure wether it is smart of me to have so many 
	 *	variables and arrays as global. Hopefully I will learn 
	 *	from this work.
	 *
	 *	OTHER FILES IN THIS FOLDER:
	 *	kode.h: 
	 *		standard headerlines, random number generators, 
	 *		basic mathematical functions (errorfuncion, 
	 *		gammafunction, etc)
	 *	Utils.cpp:
	 *		code from Jan Trulsen. Includes the Inversion 
	 *		method.
	 *	initializeDL.cpp:
	 *		from a known finite potenial leap and 
	 *		distrobutionfunctions at the boundaries 
	 *		particles are initialixed to correspond with 
	 *		a double layer.
	 *	phi_integrals.cpp:
	 *		contains the function phi(x), distrobution 
	 *		functions f(x,v) and densityintegrals n(x)
	 *	1D_plasma_simulations.cpp:
	 *		functions and procedure needed for simulating 
	 *		overall charge neutral H-plasma in equilibrium
	 *	flux.cpp:
	 *		handles particles that exit the domain and removes 
	 *		them. From integration of distrobution functions, 
	 *		the influx of new particles are found and the 
	 *		particles are injected.
	 *	diagnostics.cpp:
	 *		prints the different values, fields, etc. for 
	 *		every timestep.
	 *
	 *	FOLDERS NEEDED FOR SIMULATIONS:
	 *	Datafiles are written to specific folders. For the 
	 *	simulation to work the folders "log" and "diag".
	 *	If you also add the folder "gnuplot" and "plot" with 
	 *	"gnuplot"s original contents, you may create complete 
	 *	plotfiles in the folder "plot" by rinning the command 
	 *	"load 'all.p'" from within gnuplot. 
	 * */


#define XDIM 1
#ifdef BC_DD
	#define BCX0 0
	#define BCX1 0
#endif
#ifdef BC_vND
	#define BCX0 1
	#define BCX1 0
#endif
#ifdef BC_vNvN
	#define BCX0 1
	#define BCX1 1
#endif
#ifndef BC_vNvN
	#ifndef BC_vND
		#ifndef BC_DD
	#define BCX0 0
	#define BCX1 1
		#endif
	#endif
#endif


double me, mp, Qp, Qe, QMe, QMp, vthe, vthp, drifte, driftp;
long Ne, Np, particles;
double L, lx, ux, dt; 
double *xe, *xp, *ve, *vp, *dve, *dvp, phi0, phiDL, BCphi0, BCphiDL;
double flux_ea, flux_erd, flux_ia, flux_ird;
int N_time;

//Headerfile
#include "kode.h"
//Jans contribution
#include "Utils.cpp"
#include "MultiGrid.cpp"
//My files
#include "phi_integrals.cpp"
#include "waterbag.cpp"
#include "initializeDL.cpp"
#include "TVhist.cpp"
#include "1D_plasma_simulations.cpp"
//#include "influxhist.cpp"
#include "flux.cpp"
#include "diagnostics.cpp"
//#include "tracing.cpp"
#include "fourier.cpp"
#include "presentation.cpp"


void set_globals(int *N){

	double TionTe, mach, space, P, Adx, Avdtdx, time, num_part;

	//Variable paramters
	L=250;			//Spatial range
	phiDL=200.0;		//Potential leap (must be >0)
	mp=10.;		//Ion mass
	mach=5;			//Driftvelocity

	//For the scriptfile BCcopy
	double temp_L, temp_phiDL, temp_mp, temp_mach;
	temp_L=0; 
	temp_phiDL=0; 
	temp_mp=0;
	temp_mach=0;

	particles=long(1e5);
	time=10.;		//Simulation time in units of dt

	num_part=1e5;	//Temporary number for no. of particles


	//"SIZES" OF SIMULATION	
	space=25;		//Size of plasma added
	Adx=0.5;		//Accuracy x
	Avdtdx=0.6;		//Accuracy vdt/dx

	//PHYSICAL QUANTITIES
	me=1.0;			//Electron mass
	TionTe=1.0;		//Temperature ratio
	phi0=0.0;



	//	DO NOT CHANGE BEYOND THIS LINE
	//Total simulation time
	//N_time=int(pow(2.,P))-1;//Must be >=10 !
	//Position boundaries
	lx=0.-space; 		ux=L+space;
	//Charge
	//Code is hardcoded on the assumption of negative 
	//electrons and positive ions. Do not change!
	Qp=1./num_part;		Qe=-1./num_part;	
	//Mass
	mp/=num_part;		me/=num_part;
	//Mass charge ratio
	QMe=Qe/me;		QMp=Qp/mp;
	//Thermal velocity
	vthe=1.0;		vthp=sqrt(TionTe*me/mp);
	//Driftvelocity
	//Ion drift must be negative, electron drift must be positive
	drifte=mach*vthe;	driftp=-mach*vthp;
	//If Dirichlet boundary conditions, these are needed
	BCphi0=phi0;		BCphiDL=phiDL;
	//Automatic set spatial resolution
	double dx=1.;
	for(double p=7.; dx>Adx ; ++p){
		*N=int(pow(2.,p))+1;
		dx=(ux-lx)/((*N)-1);
		//cout<<"dx="<<dx<<" p="<<p<<" N="<<*N<<endl;
	}
	//Automatic set temporal resolution
	double vth=vthe; if(vthp>vth) vth=vthp;
	dt=Avdtdx*dx/((mach+2.)*vth);
	double vdtdx=vth*(mach+2.)*dt/dx;
	
	N_time=int(time/dt);
/*
	double t=0.;
	for(double p=8.; t<time && p<14.; ++p){
		N_time=int(pow(2.,p))-1;
		t=N_time*dt;
	//	cout<<" time: "<<t<<" "<<N_time<<" "<<p<<" "<<dt<<endl;
	}
	N_time=int(N_time);
*/
	void print_simdetails(double, int&);
	void initialstabilitycheck(double, double, double, int&);
	print_simdetails(vdtdx,*N);
	initialstabilitycheck(vdtdx,mach, TionTe, *N);

}

void finalize(bool end){
	
	finalizemother();
	printTVhist(N_time);

	if(end) exit(1);
}

int main (int argc, const char* argv[] ) {

	/*
	 *	INITIALIZATION
	 */
	system("rm diag/*.dat");	//Removing old data
	srand48(2);			//Initializing drand48
	int N;				//Size of spatial mesh
	set_globals(&N);		//Set global variables

	double RHO[N], PHI[N], E[N+1];	//Fields
	long temp=initDL(false);	//Calc. auto part numbers
	reweight(temp);			//Adjusting part. number
	temp=initDL(true);		//Initializing Double Layer


	initialdiagnostics(PHI,N_time, N);		//Diagnostics;

	//MG poisson solver initialization
	int n[XDIM]; double h[XDIM];
	initMGpoisson(n,h,N);
	MultiGrid *mg = new MultiGrid(n,h);


	//Div. test procedures
#ifdef TANHPHI
	void tanhphi();	tanhphi();	//in waterbag.cpp
#endif
#ifdef TESTMOTHER
	void testmother(); testmother();
#endif
#ifdef INIT
	cout << "Checked initialization, exiting."<<endl<<endl;
	exit(0);
#endif

	/*
	 *	TEMPORAL SIMULATIONS
	 */
	cout<<endl<<"          ---STARTING TEMPORAL SIMULATIONS---"<<endl;
	//Zeroth timestep
	charge(RHO,N);			//Assign particles to spacegrid
	poissonMG(PHI,RHO,N,mg,n,h);	//Find ES potential
	Efield(E,PHI,N);		//Find Electric field
	deltaV(E,N);			//Find change in velocity
	init_leapfrog();		//Backstepping velocities to t=-1/2
	
	diagnostics(RHO,PHI,E,N, 0);		//Diagnostics
	double phidl[N_time+1]; phidl[0]=phiDL;	//Diagnostics
	
	//Timeloop start
	for(int t=1; t<N_time+1; ++t){

		deltaV(E,N);			//Find change in vel of part
		leapfrog(t);			//Stepping pos and vel
		flux(t);			//Flux of through boundary

		charge(RHO,N);			//Assign part to spacegrid
		poissonMG(PHI,RHO,N,mg,n,h);	//Find ES potential
		Efield(E,PHI,N);		//Find Electric field

		phidl[t]=PHI[N-1]-PHI[0];	//Diagnostics
	//	trace(xp,vp,t);			//Diagnostics
		diagnostics(RHO,PHI,E,N,t);	//Diagnostics
		
	}
	
	FFTphidl(phidl,N_time+1,dt);
	
//	mg.~MultiGrid();
	finalize(false);
	
	cout<<"               -----SIMULATIONS DONE--------" << endl;
  	return 0;
}


