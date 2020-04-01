
void reweight(long temp){

	double f=double(temp)/(double(particles));
	Qp*=f;	Qe*=f;	me*=f;	mp*=f;
}




void rhoofphi(double rho[], double alpha[], int N){
	//Returns the charge density as a function of potential

	double nea(double);
	double ner(double);
	double ned(double);
	double nia(double);
	double nir(double);
	double nid(double);
	
	double p, rhosum=0.0;
	double d_phi=phiDL/(N-1);
	for(int i=0; i<N; ++i){
		p=i*d_phi;
		rho[i]=Qp*(alpha[2]*nia(p) + alpha[3]*nir(p) + 
				alpha[3]*nid(p)) 
			+ Qe*(alpha[0]*nea(p) + alpha[1]*ner(p) + 
				alpha[1]*ned(p));
		rhosum+=rho[i];
	}
	rhosum*=d_phi;





	//Diagnostics	
	double max=0.0, r, absint=0.0;	//Find maximum rho value!
	for(int i=0; i<N; ++i){
		r=abs(rho[i]); 
		if(r>max)max=r;
		absint+=r;
	}
	absint*=d_phi;

	//bc0: Charge density at x=0 taken to max charge density in DL
	//bcL: ---------"------- x=L  ------------"----------------
	//bc: Integrated charge density -----------------"------------
	double bc0=rho[0]/max*100.;
	double bcL=rho[N-1]/max*100.; 
	double bc=rhosum/absint*100.;
	
	//Print results to screen
	cout <<" Charge density" << endl;
	cout <<"   Overall neutrality= " 
				<< bc << "%" << endl;
	cout <<"   Neutrality at boundaries= " 
		<< bc0 << "%  " << bcL << "%"<<endl<< endl;
	//If essential boundary conditions are violated, print warning!	
	if(bc0*bc0>20. || bcL*bcL>20. || bc*bc>20.){
		cout 
		<<"      *****************************"<<endl
		<<"      CHARGE DENSITY BC VIOLATED!!!"<<endl
		<<"      *****************************"<<endl;
	}

}











void sagdeev(double alpha[], double phi[], 
		double (*phifunction)(double), int N){
	//Through charge density and the Sagdeev potential, 
	//calculate the last alpha so that the amplitude of the DL 
	//charge density mathches the DL potential leap phiDL. 
	//Also sets the final potential profile phi(x) wich is stored in 
	//global phi_array accessible through function "double phi(x)".
	//Returns the array alpha after adjusted to its final values.

	//Finding rho(phi)
	double rho[N];	rhoofphi(rho,alpha,N);

	//A fix to avoid V(phi)<0 due to big number subtractions
//	for(int i=0; i<5; ++i) if(rho[i]>0.) rho[i]=-rho[i];
//	for(int i=N-1; i>N-6; --i) if(rho[i]<0.) rho[i]=-rho[i];
//	ofstream r("diag/rho.dat");
//	double dp=phiDL/(N-1);
//	for(int i=0; i<N; ++i) r << i*dp << " " << rho[i]<< endl;
//	r.close();

	//Finding the Sagdeev potential V(phi)
	//from the charge density
	double V[N];
	double d_phi=phiDL/(N-1);
	V[0]=0.;
	for(int i=1; i<N;++i){
		V[i] = V[i-1] - 0.5*(rho[i-1]+rho[i])*d_phi;
		if(V[i]*V[i]<1e-14 && V[i]<0.0) V[i]=-V[i];
	}
	for(int i=1; i<N; ++i) 	//Test if V[i]<0. If so, print 
		if(V[i]<=0.0){ 	//warning and exit program
			cout <<"  **ERROR**: Wrong Sagdeev potential: "; 
			cout <<V[i] << " " << rho[0] << " "<< i << endl;
			cout <<"Forcing exit!!"<<endl<<endl;
			exit(0);
		}
	

	
	//Calculate alpha raising the magnitude of the DL charge density
	//to a level matching phiDL.
	double I[N];
	I[0]=0.0; 	//Special treatment of singularity
	I[1]=2.*sqrt(-d_phi/.5/rho[1]);	
	if(isnan(I[1])){
		cout << "     I[1] is Nan!!!!!"<<endl; 
		exit(0);
	}
	for(int i=2;i<N-1;++i) I[i]=I[i-1] + d_phi/sqrt(V[i]);	//Loop
	I[N-1]=I[N-2] + 2.*sqrt(d_phi/.5/rho[N-2]);		//Boundary
	double alpha_s=0.5*I[N-1]*I[N-1]/L/L;			//Result


	//find x(phi)
	double x[N];
	double a=1./sqrt(2.*alpha_s);
	x[0]=0.;
	for(int i=0; i<N-1; ++i) x[i]=a*I[i];
	x[N-1]=L;
	

	//inverting to phi(x)
	phi[0]=0.;
	int j;
	double w, xi;
	double dx=L/(N-1);
	j=1;
	for(int i=1;i<N;++i){
		xi=i*dx;
		//j=0;
		while(xi>x[j])++j;
		if(j>N-1)	cout << "j too large: " << j << endl;
		
		w=(xi - x[j-1])/(x[j] - x[j-1]);
		phi[i] = (j-1)*d_phi + w*d_phi;

		if(phi[i]<0.0 || phi[i]>phiDL){
			cout << "INCORRECT phi[i]="<<phi[i]<<endl;
			cout << "Forcing exit!"<<endl;
			exit(0);
		}
	}
	
	//Returnvalues
	for(int i=0; i<4; ++i) alpha[i]*=alpha_s;


	//For diagnostics: find rho(x)
	double rhox[N], p;
	for(int i=0; i<N; ++i){
		p=(*phifunction)(i*dx);
		rhox[i]=Qp*(alpha[2]*nia(p)+alpha[3]*nir(p)+
				alpha[3]*nid(p))
			+Qe*(alpha[0]*nea(p)+alpha[1]*ner(p)+
				alpha[1]*ned(p));
	}

//	double rhoxx[N], dxdx=dx*dx*100;
//	rhoxx[0]=0;
//	rhoxx[N-1]=0;
//	for(int i=1; i<N-1; i+=10){
//		rhoxx[i]=-(phi[i-10]+phi[i+10]-2.*phi[i])/dxdx;
//		
//	}

	//Diagnostics
	ofstream sag("diag/supermethod.dat");
	for(int i=0;i<N;++i)sag << i*dx << " " 
				<< i*d_phi << " " 
				<< V[i] << " " 
				<< x[i] << " "
				<< phi[i] << " "
				<< rho[i] << " "
				<< rhox[i] << " "
//				<< rhoxx[i] << " "
				<< I[i] << " " 
				<< endl;
	sag.close();

}
















void BC_ChargeNeutrality(double alpha[]){
	double nea(double);
	double ner(double);
	double ned(double);
	double nia(double);
	double nir(double);
	double nid(double);
	
	double alpha_ir=1.0, alpha_er=1.0;
	double alpha_ea=1.0, alpha_ia=1.0;

	int N=200;
	double d_phi=phiDL/(N-1);

	//Integrated density over space for all six sorts
	double Nea=0.0, Ner=0.0, Ned=0.0, Nia=0.0, Nir=0.0, Nid=0.0;
	for(double p=0.; p<phiDL; p+=d_phi){
		Nea+=nea(p)*d_phi;	
		Ner+=ner(p)*d_phi;	
		Ned+=ned(p)*d_phi;	
		Nia+=nia(p)*d_phi;	
		Nir+=nir(p)*d_phi;	
		Nid+=nid(p)*d_phi;	
	}
	
	//Density of the six particle sorts taken at boundaries
	double nea_0=nea(0);
	double nea_L=nea(phiDL);
	//ner_0=0
	double ner_L=ner(phiDL);
	double ned_0=ned(0);
	double ned_L=ned(phiDL);
	double nia_0=nia(0);
	double nia_L=nia(phiDL);
	double nir_0=nir(0);
	//nir_L=0
	double nid_0=nid(0);
	double nid_L=nid(phiDL);
	


	double aea=1.0, aer, aia, air;
	//alpha of reflected electrons
	aer=aea*(   	(Nea*nia_L/Nia-nea_L) + 
			(nea_0-Nea*nia_0/Nia)*(nid_L-(Nir+Nid)*nia_L/Nia)
			/(nir_0+nid_0-(Nir+Nid)*nia_0/Nia)	)/(   
			(ner_L+ned_L-(Ner+Ned)*nia_L/Nia) - 
			(ned_0-(Ner+Ned)*nia_0/Nia)
			*(nid_L-(Nir+Nid)*nia_L/Nia)
			/(nir_0+nid_0-(Nir+Nid)*nia_L/Nia)	);
	//alpha of reflected ions
	air=(aer*(ned_0-(Ner+Ned)*nia_0/Nia)+aea*(nea_0-Nea*nia_0/Nia))
		/(nir_0 + nid_0-(Nir+Nid)*nia_0/Nia);
	//alpha of accelerated ions
	aia= (aea*Nea + aer*(Ner+Ned) - air*(Nir+Nid))/Nia;

	//Storing alpha results in return vector alpha
	alpha[0]=aea;	alpha[1]=aer;	alpha[2]=aia;	alpha[3]=air;



	//Diagnostics
	for(int i=0; i<4; ++i) if(alpha[i]<0.){
		cout 	<< "ERROR: NEGATIVE PARTICLE NUMBERS!!"
		<<" Forcing exit!"<<endl<<endl;
		exit(0);
	}

}















void fill(long *Npart, double x[], double v[], 
	long Nfill, double temp_x[], double temp_v[]){ 
	
	for(long i=0; i<Nfill; ++i){
		x[*Npart+i] = temp_x[i];
		v[*Npart+i] = temp_v[i];
	}
	*Npart+=Nfill;
}




























void particlenumbers(double alpha[], long N[]){
	//Calculating the number of particles of specie 's' to be
	//generated from N_s=\int_0^phiDL alpha_s n_s(phi) d_phi

	double nea(double);
	double ner(double);
	double ned(double);
	double nia(double);
	double nir(double);
	double nid(double);

	double ea=0., er=0., ed=0, ia=0., ir=0., id=0.;
	double dx=(ux-lx)/300;
	double p;
	//Integral of density functions for the six species
	for(double x=lx;x<=ux;x+=dx){
		p=phi(x);
		ea+=nea(p)*dx;	
		er+=ner(p)*dx;
		ed+=ned(p)*dx;
		ia+=nia(p)*dx;
		ir+=nir(p)*dx;
		id+=nid(p)*dx;
	}
	
	//Multiplying with alpha to get particle numbers
	N[0]=long(alpha[0]*ea);	//Accelerated electrons
	N[1]=long(alpha[1]*er);	//Reflected electrons
	N[2]=long(alpha[1]*ed);	//Decelerated electrons

	N[3]=long(alpha[2]*ia);	//Accelerated ions
	N[4]=long(alpha[3]*ir);	//Reflected ions
	N[5]=long(alpha[3]*id);	//Decelerated ions
	
	//Global variables telling the number of particles generated 
	//at each moment. Naturally they are equal to 0 now.
	Ne=0; 	Np=0; 

	//Diagnostics output
	cout 	<<" Particle numbers"<< endl
	 	<<"   Nea="<<N[0]<<"  Ner="<<N[1]<<"  Ned="<<N[2] 
		<<"  Nia="<<N[3]<<"  Nir="<<N[4]<<"  Nid="<<N[5] 
		<<endl 
		<<"   Ne="<<N[0]+N[1]+N[2]<<"   Ni="<<N[3]+N[4]+N[5] 
		<<"   Ni-Ne="<<(-(N[0]+N[1]+N[2])+(N[3]+N[4]+N[5]))
		<<"   (Ni-Ne)/(avg(N))=" 
		<<double(-(N[0]+N[1]+N[2])+(N[3]+N[4]+N[5]))
			/double(+(N[0]+N[1]+N[2])+(N[3]+N[4]+N[5]))*50.
		<< "%"<< endl<< endl;
}















void allocate_xvarrays(long ne, long ni){
	//Allocating arrays to global pointers *xe,*ve,*xp,*vp,*dve,*dvp
	//Stores positions, velocities and change in velocities for 
	//electrons and protons

	//Setting size of arrays larger than the number of particles I'm 
	//starting with, to allow a positive net flux of particles 
	long alloc_e=long(ne*1.5);
	long alloc_i=long(ni*1.5);

	//Allocating
	xe = new double[alloc_e];		
	ve = new double[alloc_e];		
	xp = new double[alloc_i];		
	vp = new double[alloc_i];		
	dve = new double[alloc_e];		
	dvp = new double[alloc_i];		
	
	//Zeroing all elements of all arrays before use.
	for(long i=0; i<alloc_e; ++i){ xe[i]=0.0;  ve[i]=0.0; dve[i]=0.0;}
	for(long i=0; i<alloc_i; ++i){ xp[i]=0.0;  vp[i]=0.0; dvp[i]=0.0;}
}
















void generate(long K[], double alpha[]){
	//From known particle distrobution functions f(x,v) use the 
	//inversionmethod to decide how particles are distributed in 
	//phasespace, then draw them after this pattern. Last step stores 
	//all drawn particles in global arrays for position/velocity for 
	//both particle species.

	double *x, *y, fnorm;

	long N=0; for(int i=0; i<6; ++i) if(K[i]>N) N=K[i];
	N=long(N+1);
	
	x = new double[N];	//Position temp array
	y = new double[N];	//Velocity temp array

	void ps_plot(string, double [], double [], long, int);
	const int res=300;	//Space/velocity gridsize in inversion
	long start, stop;
/*
	void histogram(long [], double (*n)(double), double , 
		double (*f)(double,double), double [], 
		double [], long , double );
	void histogram(long [],long [],long [],long [],int);
	long histea[200]; 
	long hister[200]; 
	long histed[200]; 
	long histia[200]; 
	long histir[200]; 
	long histid[200]; 
*/	
	//For all six particle species, particles are drawn from 
	//their respective PDF's, put into temporary arrays x and y, 
	//then placed in common arrays for electrons and ions 
	//respectively.
	cout<<" Generating particles" << endl;
	//ACCELERATED ELECTRONS
	cout<<"   Electrons: ";
	cout << "Accelerated " ;
	N=K[0];	
	Inversion ea(lx,ux,res,lea(lx),uea(ux),res,lea,uea,fea,fnorm);
	ea.draw(N, x, y);
//	ps_plot("genea-x-v",x,y,N,0);
//	histogram(histea, nea, alpha[0], fea, x, y, N, vthe);
	fill(&Ne,xe,ve,N,x,y);
	
	//REFLECTED ELECTRONS
	cout << "- Reflected " ;
	N=K[1];	
	Inversion er(lx,ux,res,ler(ux),uer(ux),res,ler,uer,fer,fnorm);
	er.draw(N, x, y);
//	ps_plot("gener-x-v",x,y,N,0);
//	histogram(hister, ner, alpha[1], fer, x, y, N, vthe);
	fill(&Ne,xe,ve,N,x,y);

	//DECELERATED ELECTRONS
	if(K[2]!=0){
		cout << "- Decelerated" ;
		N=K[2];	start=stop; 	stop+=N;
		Inversion ed(lx,ux,res,led(ux),ued(lx),res,led,ued,fed,fnorm);
		ed.draw(N, x, y);
//		ps_plot("gened-x-v",x,y,N,0);
//		histogram(histed, ned, alpha[1], fed, x, y, N, vthe);
		fill(&Ne,xe,ve,N,x,y);
	}
	cout<<endl<<"   Ions: ";
	//ACCELERATED IONS
	cout << "Accelerated ";
	N=K[3];	
	Inversion ia(lx,ux,res,lia(lx),uia(ux),res,lia,uia,fia,fnorm);
	ia.draw(N, x, y);
//	ps_plot("genia-x-v",x,y,N,0);
//	histogram(histia, nia, alpha[2], fia, x, y, N, vthp);
	fill(&Np,xp,vp,N,x,y);
	
	//REFELECTED IONS
	cout << "- Reflected " ;
	N=K[4];	
	Inversion ir(lx,ux,res,lir(lx),uir(lx),res,lir,uir,fir,fnorm);
	ir.draw(N, x, y);
//	ps_plot("genir-x-v",x,y,N,0);
//	histogram(histir, nir, alpha[3], fir, x, y, N, vthp);
	fill(&Np,xp,vp,N,x,y);

	//DECELERATED IONS
	if(K[5]!=0){
		cout << "- Decelerated";
		N=K[5];	
		Inversion id(lx,ux,res,lid(ux),uid(lx),res,lid,uid,fid,fnorm);
		id.draw(N, x, y);
//		ps_plot("genid-x-v",x,y,N,0);
//		histogram(histir, nid, alpha[3], fid, x, y, N, vthp);
		fill(&Np,xp,vp,N,x,y);
	}
	cout<<endl;
//	histogram(histea, hister, histia, histir, 200);
	delete x,y;
}










long initDL(bool gen){
	cout 	<< endl<< "                    "
		<< "--- INITIALIZING DL ---" << endl;
	//Array N: number of particles	
	//	[0] : accelerated electrons
	//	[1] : reflected electrons
	//	[2] : decelerated electrons
	//	[3] : accelerated ions
	//	[4] : reflected ions
	//	[5] : decelerated ions
	long N[6];

	//Array alpha: weighting of PDF's
	//	[0] : accelerated electrons
	//	[1] : reflected and decelerated electrons
	//	[2] : accelerated ions
	//	[4] : reflected and decelerated ions
	double alpha[4];
	

//	densitycheck(alpha);	exit(0);
	
	//The boundarycondition of having overall charge neutrality
	//and charge neutrality at the boundaries gives 3 equations 
	//for solving the alphas, wich is done here
	BC_ChargeNeutrality(alpha);

	//Through a calculation of the Sagdeev potential, 
	//the final criterium for alpha is laid and 
	//the potential profile consistent with a 
	//DL is found. After this, phi_array is accessible through
	//the function "double phi(double x)".
	sagdeev(alpha, phi_array, phi, phi_array_length);
//	a_simple_test(alpha);	

	//Calculating the expected influx number at both boundaries for 
	//both particle species. Set the value of the global variables 
	//flux_ea, flux_erd, flux_ia and flux_ird.
	void initflux(double []);
	initflux(alpha);

	//Integrating the particle densities with the weightings alpha
	//to find the number of each particle specie. Also setting global 
	//variables Ne and Np (total electron and ion number)
	particlenumbers(alpha,N);
	
	if(gen){
		//Allocate the global arrays for 
		//electron/ion position/velocity
		allocate_xvarrays(N[0]+N[1]+N[2],N[3]+N[4]+N[5]);
	
		//Draw particles from PDF's, 
		//in the amount described by N, 
		//and place them in global position/velocity arrays.
		generate(N,alpha);

		return Ne;
	}else{ 	
		cout <<" Adjusting particle numbers" << endl;
		return N[0]+N[1]+N[2]; 
	}
}



























	/*
	 *	HISTOGRAMS
	 */
/*
void histogram(long hist[], double (*n)(double), double alpha, double (*f)(double,double), double x[], double v[], long Nxv, double vth)
{

	cout << " Nxv=" << Nxv << endl;
	//Maximum velocity
	double vmax=0.0;
	for(long i=0; i<Nxv; ++i) if(abs(v[i])>vmax) vmax=abs(v[i]);
7
	//Position histogram initializations
	int xbins=200;
	double xbin=(ux-lx)/(xbins-1);
	long xhist[xbins];
	for(int m=0; m<xbins; ++m)	xhist[m]=0;
	
	//Velocity histogram initializations
	int vbins=xbins;
	double vbin=2.*vmax/(vbins-1);
	long vhist[vbins];
	for(int m=0; m<vbins; ++m)	vhist[m]=0;
	
	//Indexvalues
	int ix,iv;

	//Filling histograms
	for(long i=0;i<Nxv;++i){
		ix=int((x[i]-lx)/xbin);
		iv=int( (v[i]+vmax)/(vmax+vmax)*(vbins-1)+0.5);
		if(ix<0 || iv<0 || ix>xbins || iv>vbins)
			cout	<<"  HI ERROR!!  "
				<< ix << " " << iv << "   X   ";
		xhist[ix]+=1;
		vhist[iv]+=1;
	}

	//Check particle numbers
	int sumx=0, sumv=0;
	for(int i=0; i<xbins; ++i){sumv+=vhist[i]; sumx+=xhist[i];}
	double dp=phiDL/(xbins-1);
	double dx=(ux-lx)/(xbins-1);
	double sumn=0., sumnx=0.;
	for(double p=0.;p<=phiDL;p+=dp) sumn+=alpha*(*n)(p)*dp;
	for(double X=lx; X<=ux; X+=xbin){ 
		if(X<=ux)dp=phi(X+xbin) - phi(X);
		else	dp=phi(X) - phi(X-xbin);
		sumnx+=alpha*(*n)(phi(X))*dx; 
	}
	cout << "Histogram paticlenumber comparison:" << endl;
	cout 	
		<< "  in density functions(x)   : " << sumnx << endl
		<< "  in density functions(phi) : " << sumn << endl
		<< "  in particle number N      : " << Nxv << endl
		<< endl;


	double X, V, DP; double P=0.;
	ofstream ofile("diag/histogram.dat");
	for(int i=0; i<xbins; ++i){
		X=i*xbin+lx;
		V=i*vbin-vmax;
		if(i<xbins)DP=phi(X+xbin)-phi(X);
		else	DP=phi(L)-phi(L-xbin);

		ofile 	<< X << " " 
			<< xhist[i]/dx << " "
			<< alpha*(*n)(phi(X)) << " "
			<< V << " " 
			<< vhist[i] << " " 
			<< endl;
	P+=dp;
	}	
	ofile.close();

	for(int i=0; i<200; ++i) hist[i] = xhist[i];
	cout << "Pause at histogram end. \n Do your plots,"
		<<" type something here and press enter to continue." 
		<< endl << "  :";
	string tull; cin >> tull;
}

void histogram(long histea[], long hister[], long histia[], 
		long histir[], int bins){

	ofstream fout("diag/histcomp.dat");
	for(int i=0; i<bins-1; ++i)
	fout 	<< i << " "
		<< histea[i] << " " 
		<< hister[i] << " " 
		<< histia[i] << " " 
		<< histir[i] << " " 
		<< Qp*(histia[i] + histir[i]) 
			+ Qe*(histea[i] + hister[i]) << " "
		<< endl;
	fout.close();
}
*/




/*
void densitycheck(double a[]){
	double p=phiDL/1.1;
	nea(p); ner(p); ned(p); nia(p); nir(p); nid(p);
}
*/




/*
void a_simple_test(double alpha[]){
	ofstream f("diag/initdensity.dat");
	double dx=(ux-lx)/127, p;
	for(double x=0.;x<=L;x+=dx){
		p=phi(x);
		f<< x << " "
		<< alpha[0]*nea(p) << " " 	
		<< alpha[1]*ner(p) << " "
		<< alpha[1]*ned(p) << " "
		<< alpha[2]*nia(p) << " " 	
		<< alpha[3]*nir(p) << " "
		<< alpha[3]*nid(p) << " "
		<< Qp*(alpha[0]*nia(p)+alpha[1]*nir(p)+alpha[1]*nid(p))
		   +Qe*(alpha[2]*nea(p)+alpha[3]*ner(p)+alpha[3]*ned(p))
		<<endl;
	}
	f.close();
	double dp=phiDL/499;
	double ea=0., er=0., ed=0., ia=0., ir=0., id=0.;
	for(double p=0.; p<=phiDL; p+=dp){
		ea+=alpha[0]*nea(p)*dp;
		er+=alpha[1]*ner(p)*dp;
		ed+=alpha[1]*ned(p)*dp;
		ia+=alpha[2]*nia(p)*dp;
		ir+=alpha[3]*nir(p)*dp;
		id+=alpha[3]*nid(p)*dp;
	}
	cout 	<< (long)ea << " " << (long)er << " " << (long)ed 
		<< " " << (long)ia << " " << (long)ir << " " 
		<< (long)id << endl;
}
*/

