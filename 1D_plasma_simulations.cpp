

void setDLpotential(double pot[], int N){
	double dx=(ux-lx)/(N-1);
	for(int i=0; i<N; ++i){
		pot[i]=phi(i*dx+lx);
	}
//	PHIDL=phiDL;

}







void charge(double Qcell[], int N){
	//Applying the particle-in-cell (PIC) scheme to the particles 
	//position. All continuos positions are related to a grid point
	//through the PIC scheme Cloud-in-cell. Normalization is included 
	//the factor (dx)^2 wich is not a part of the charge density. I am 
	//doing this to save simulation time by not adding them later in 
	//my poisson solver.
	
	int index;
	double x;
	double dx=(ux-lx)/(N-1);
	
	for(int i=0; i<N; ++i) Qcell[i]=0.0;

	//Electrons
	for(int i=0; i<Ne; ++i){
                x=(xe[i]-lx)/dx;	
		index = int(x);	//index=grid point adress
		//x is now the distance to the nearest lower gridpoint, 
		//a number between 0 and 1 where 1 is the distance 
		//between gridpoints
		x-=index;	

		//Assigning a particles charge as it was a cload to points
	    	Qcell[index]+=(1.0-x)*Qe;
		Qcell[index+1]+=x*Qe;
	}
	//Protons (see comments above)
	for(int i=0; i<Np; ++i){
        	x=(xp[i]-lx)/dx;
		index = int(x);
		x-=index;

		Qcell[index]+=(1.0-x)*Qp;	
		Qcell[index+1]+=x*Qp;
	}
	//Boundary conditions
        Qcell[0]*=2.0;		Qcell[N-1]*=2.0;

	for(int i=0; i<N; ++i) Qcell[i]/=dx;

	lowpassfilter(Qcell,N);
}


void initMGpoisson(int n[], double h[], int N){
	double dx=(ux-lx)/(N-1);
	for(int i=0;i<XDIM;++i)	n[i]=N,	h[i]=dx;
	//Assuming XDIM == 1

}


void poissonMG(double PHI[], double RHO[], int N, MultiGrid *mg, 
	int n[], double h[]){

	//Assuming XDIM == 1
	double rhoMG[n[0]+2], phiMG[n[0]+2];
	//Copying RHO to array of n[0]+2 dim and reversing sign
	for(int i=1; i<n[0]+1; ++i) rhoMG[i]=-RHO[i-1];
	
	//Boundary conditions
#if BCX0 == 0	//If phi(x=0) : Dirichlet BC 
	rhoMG[0]=BCphi0;	//phi=0
#endif
#if BCX0 == 1	//If phi(x=0) : von Neuman 
	rhoMG[0]=0.;		//dphi/dx=0
#endif
#if BCX1 == 0	//If phi(x=L) : Dirichlet BC 
	rhoMG[n[0]+1]=BCphiDL;	//phi=phiDL
#endif
#if BCX1 == 1	//If phi(x=L) : von Neuman 
	rhoMG[n[0]+1]=0.;	//dphi/dx=0
#endif

	(*mg).mglin(phiMG,rhoMG,n,h);
	(*mg).delbc2rho(rhoMG,n,h);
	
	//Copy back to PHI-array
#ifdef BC_vNvN
	double phimid=0.0;
	phimid=phiMG[n[0]/2]-BCphiDL/2.;
	for(int i=1; i<n[0]+1; ++i) PHI[i-1]=phiMG[i]-phimid;
#else
	for(int i=1; i<n[0]+1; ++i) PHI[i-1]=phiMG[i];
#endif
	phi0=PHI[0];		phiDL=PHI[N-1];
	if(isnan(phi0)) cout << "phi0 is NaN" << endl;
	if(isnan(phiDL)) cout << "phiDL is NaN" << endl;

#ifdef PTEST
	ofstream pout("diag/poissonMG.dat");
	double x;	int i;
	for(i=1, x=lx;i<=n[0];++i,x+=dx)
		pout<<x<<" "<<-rhoMG[i]<<" "<<phiMG[i]<<endl;
	pout.close();
	cout<<"Poissontest forcing exit!"<<endl; exit(0);
#endif
	
}





void Efield(double E[], double PHI[], int N){
	//Find the electric field with the first derivative of the 
	//potential. E is calculated as a midpoint in i+1/2 not in 
	//i compared to the grid of Qcell and phi. Notice that E is 
	//a matrix with two more elements than needed to fill the 
	//midpoints mentioned above. It is so to have the boundaries 
	//in the arrays first and last element.
	
	double dx=(ux-lx)/(N-1);
	
	for(int i=1; i<N; ++i)	E[i]=-(PHI[i]-PHI[i-1])/dx;

	//Boundary conditions
	//Extrapolating E's outside of the domain
	E[0]=2*E[1]-E[2];
	E[N]=2*E[N-1]-E[N-2];
}




void deltaV(double E[], int N){
	//From the electric field I find the acceleration multiplied with 
	//the timestep (wich means this is the change in velocity pr 
	//timestep. Electric field is mapped back to particle position 
	//in the same fassion the particles positions where mapped onto 
	//the grid earlier. The procedures are still different since 
	//E is calculated in midpoints (i+1/2).
	
	int index;
	double x; 
	double dx=(ux-lx)/N;
	double dv_factor_p=dt*QMp;	//Velocitystep, protons
	double dv_factor_e=dt*QMe;		//Velocitystep electrons
	
	//Electrons (for comments on procedure, see "void PIC()"
	for(int i=0; i<Ne; ++i){

		x=(xe[i]-lx)/dx-0.5;	
		index = int(floor(x))+1;		
		x-=index;

		//CIC weighting of Efields in gridpoints to particles pos.
		dve[i]=(E[index]*(1.0-x)+E[index+1]*x)*dv_factor_e;
	}		
	
	//Protons
	for(int i=0; i<Np; ++i){
		
		x=(xp[i]-lx)/dx-0.5;
		index = int(floor(x))+1;
		x-=index;
		
		//CIC weighting of Efields in gridpoints to particles pos.
		dvp[i]=(E[index]*(1.0-x)+E[index+1]*x)*dv_factor_p;
	}
}





void init_leapfrog(){
	//Backstepping one half timestep (Leapfrog)
	
	for(int i=0; i<Ne; ++i)	
		ve[i] -= dve[i]/2.0; 

	for(int i=0; i<Np; ++i)	
		vp[i] -= dvp[i]/2.0;
}











void leapfrog(int time){
	//Calculating v(t+dt/2) and x(t+dt) (leapfrog)
	//Leapfrog: velocities found in midpoints between 
	//those of positions and accelerations
#ifdef TVSELECT	
	double xold,x,v[1];
#endif
	for(long i=0; i<Ne; ++i){
#ifdef TVSELECT
		xold=xe[i];
#endif
		ve[i]+=dve[i]; 	
		xe[i]+=ve[i]*dt;
#ifdef TVSELECT
		x=xe[i];
		//In TVhist particle register
		if((xold<xlowlimit && xlowlimit<x) || 
			(x<xuplimit && xuplimit<xold)){
				v[0]=ve[i];
				inTVhist(time,v,vine, tine, i, 1);
		}

		//Out TVhist particle register
		if((xold>xlowlimit && xlowlimit>x) || 
			(x>xuplimit && xuplimit>xold))
				outTVhist(time,i,Ne,vine,tine,htve,lve,uve);
#endif
	}
	
	for(long i=0; i<Np; ++i){
#ifdef TVSELECT
		xold=xp[i];
#endif
		vp[i]+=dvp[i];	
		xp[i]+=vp[i]*dt;
#ifdef TVSELECT
		x=xe[i];
		//In TVhist particle register
		if((xold<xlowlimit && xlowlimit<x) || 
			(x<xuplimit && xuplimit<xold)){
				v[0]=vp[i];
				inTVhist(time,v,vinp, tinp, i, 1);
		}

		//Out TVhist particle register
		if((xold>xlowlimit && xlowlimit>x) || 
			(x>xuplimit && xuplimit>xold))
				outTVhist(time,i,Np,vinp,tinp,htvp,lvp,uvp);
#endif
	}
}



