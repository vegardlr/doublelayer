void finalize(bool);


int diagcounter=0;
int carpetcounter=0;

void initialdiagnostics(double PHI[], int N_time, int N){

	void densityplot(double [], int, int);
	void initTVhist(long,long);
	void resetTVhist();
	//void reset_influxhist();
	//void selecttrace(double [], double [], long, double);
	void gnuplotpartnum();

	densityplot(PHI,N,2);
	initTVhist(Ne, Np);
	resetTVhist();
	//reset_influxhist();
	//selecttrace(xp,vp,Np,QMp*phiDL);
	gnuplotpartnum();

	void printf(double (*fa)(double, double), double (*la)(double), 
		double (*ua)(double), double (*fr)(double,double),
		double (*lr)(double), double (*ur)(double), 
		double (*fd)(double, double), double (*ld)(double), 
		double (*ud)(double), string specie);

	cout<<endl<<" Average no. of particles pr. cell = "<<double(Ne+Np)/(2.*N)<<endl<<endl;
	cout<<"Printf takes time";
	printf(fea, lea, uea, fer, ler, uer, fed, led, ued, "e");
	printf(fia, lia, uia, fir, lir, uir, fid, lid, uid, "i");
	cout<<", done."<<endl;


	ofstream phiout("diag/initialphi.dat");
	double x,dx=(ux-lx)/(N-1);
	for(x=lx;  x<=ux; x+=dx) phiout<<x<<" "<<phi(x)<<endl;
	phiout.close();
}





void gnuplotpartnum(){
	ofstream fout("diag/partnum.dat");
	fout<<" Ne="<<Ne<<" Ni="<<Np<<endl;
}



void diagnostics(double RHO[], double PHI[], double E[], int N, int time){

	void checkfornan(double [], double [], double [], int);
	void phidlprintout(int,double);
	void fields(double [], double [], double [], int, int);
	void densityplot(double [], int, int);
//	void idlpsplot(string, int , double [], double [], int);
	void ps_plot(string, double [], double, double, 
		double [], double, double, long, int);
	void boundaryvdf(double x[], double v[], long Nxv, 
		double (*fa)(double, double), double la, double ua, 
		double (*frd)(double, double), double lrd, double urd, 
		int, string);
	void fieldcarpet(int,double [], double [], double [], int);
	void gnuplottimestamps(double);
	void gnuplotmark(string,string);
	
	void mother(double*, int, int);
	mother(E,N,time);


	phidlprintout(time, (PHI[N-1]-PHI[0]) );	


	//Data written 200 times during simulations
	if(time*100/N_time==carpetcounter){
		fieldcarpet(time,RHO,PHI,E,N);
		++carpetcounter;
	}


	//Data written 10 times during simulations
	//plottong routines based on the fact that the number 10 is 10!
	if(time*10/N_time==diagcounter){
		cout << " Diagnostics at t="<<time;

		gnuplottimestamps(time*dt);

		//print_influxhist(0.,uid(lx),time,"ir");
		checkfornan(RHO,PHI,E,N);//Check ALL value for NaN
		fields(RHO, PHI, E, N, time);
		densityplot(PHI, N, time);
		ps_plot("xe-ve",xe,lx,ux,ve,led(ux),uea(ux),Ne,time);
		ps_plot("xp-vp",xp,lx,ux,vp,lia(lx),uid(lx),Np,time);

//		stabilitycheck(RHO,PHI,E,N);
//		idlpsplot("electron",time,xe,ve,Ne);
		boundaryvdf(xe,ve,Ne,fea,lea(0.),
			uea(0.),fer,led(L),uer(L),time,"e");
		boundaryvdf(xp,vp,Np,fia,lia(L),
			uia(L),fir,lir(0.),uid(0.),time,"p");

		diagcounter++;
		cout << " ...Finished"<<endl;
	}

}

void gnuplottimestamps(double time){

	char snu[13]; sprintf(snu, "%d", diagcounter);
	string str="gnuplot/marks/time";
	str = str+snu;
	ofstream fout(str.c_str());
	fout.precision(3);
	fout<<time<<endl;
	fout.close();

}


void gnuplotmark(string name, string file){

	char snu[13]; sprintf(snu, "%d", diagcounter);
	string str="gnuplot/marks/";
	str = str+name+snu;
	ofstream fout(str.c_str());
	fout<<"../"<<file;
	fout.close();

}




void fieldcarpet(int time, double RHO[], double PHI[], 
			double E[], int N){
	ofstream carpet("diag/carpet.dat", ios::app);

	double x, dx=(ux-lx)/(N-1);
	int i; 
	for(i=0, x=lx; i<N; ++i,x+=dx)
		carpet<<time*dt<<" "<<x<<" "
			<<RHO[i]<<" "<<PHI[i]<<" "<<E[i]<<endl;
	
	carpet<<endl;
	carpet.close();
}


void phidlprintout(int time, double dl){
	//phiDL printout
	ofstream pout("diag/phidl.dat",ios::app);
	pout<<time*dt<<" "<<dl<<endl;
	pout.close();
	if(phiDL<0){ 
		cout<<endl<<" ERROR: phiDL<0"<<endl; 
		finalize(true);
	}
}


void fields(double RHO[], double PHI[], double E[], int N, int time){

	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/fields.";
	str = str+snu+".dat";
	ofstream fout(str.c_str());

	gnuplotmark("fields",str);

	double dx=(ux-lx)/(N-1);
	double x;
	for(int i=0; i<N; ++i)
		fout 	<<  i*dx+lx << " " 
			<< RHO[i] << " " 
			<< PHI[i] << " " 
			<< E[i] << " " 
			<< endl;

	fout.close();


}







void ps_plot(string specie, double x[], double xl, double xu, 
	double v[], double vl, double vu, long N, int time){

	int xbins=200, vbins=xbins;
	double dx=(xu-xl)/(xbins-1), dv=(vu-vl)/(vbins-1);
	double hist[xbins][vbins];
	int ix,iv;
	double X,V;
	//Zeroing histogram
	for(ix=0; ix<xbins; ++ix) for(iv=0;iv<vbins;++iv) hist[ix][iv]=0.;
	//Filling histogram
	for(long i=0;i<N;++i)
		++hist[ int((x[i]-xl)/dx) ][ int((v[i]-vl)/dv) ];
	
	
	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/ps.";
	str = str+specie+"."+snu+".dat";
	ofstream ps(str.c_str());

	gnuplotmark("ps"+specie,str);

	for(ix=0, X=xl; ix<xbins; ++ix, X=ix*dx+xl){ 
		for(iv=0, V=vl; iv<vbins; ++iv, V=iv*dv+vl){
			ps << X << " " << V << " "
				<< hist[ix][iv] << endl;;
		}
		ps << endl;
	}
}





/*
void ps_plot(string specie, double x[], double v[], long N, int time){
	//Phase space plottable datafiles	
	//Set filename as string	
	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/ps.";
	str = str+specie+"."+snu+".dat";

	//Open file, write to file
	ofstream ps(str.c_str());
	for(long i=0;i<N;++i) ps << x[i] << " " << v[i] << endl;
	ps.close();
}
*/











/*
void idlpsplot(string specie, int time, double x[], double v[], int Nxv){
	
	string str="data/faserom.";
	str=str+specie+".dat";
	ofstream fout(str.c_str(), ios::app);
	fout 	<< specie << endl 
		<< (ux-lx) << endl
		<< time << endl 
		<< Nxv << endl;
	for(long i=0; i<Nxv; ++i) fout << x[i] << " " << v[i] << endl;
}
*/























double min(double a[], long N){
	double m=1e300;
	for(long i=0; i<N; ++i) if(a[i]<m) m=a[i];
	return m;
}
double max(double a[], long N){
	double m=-1e300;
	for(long i=0; i<N; ++i) if(a[i]>m) m=a[i];
	return m;
}


void boundaryvdf(double x[], double v[], long Nxv, 
	double (*fa)(double, double), double la, double ua, 
	double (*frd)(double, double), double lrd, double urd, 
	int time, string specie){
	//VDF AT BOUNDARIES
	
	double vl=min(v,Nxv);
	double vu=max(v,Nxv);
	//vl=-0.5;
	//vu=0.5;
	//cout << "Histogram vdf at boundary: "<< specie << endl;
	//cout << "vl and vu: "<< vl << " " << vu << endl;
	
	
	int bins=100;
	double bin=(vu-vl)/(bins-1);
	long vhu[bins];		//Velocity histogram lower
	long vhl[bins];		//Velocity histogram upper
	for(int i=0; i<bins; ++i){
		vhl[i]=0;
		vhu[i]=0;
	}

	double f=0.002;
	int iv;
	for(long i=0; i<Nxv; ++i){
		iv=int((v[i]-vl)/(vu-vl)*(bins-1));
		if(iv>=0 && iv<bins){
		if(x[i]<lx+f*L){
			if(x[i]<lx) cout << "vdfplot: x<0: "<<x[i]<<" "<<v[i]<<endl;
			vhl[iv]+=1;
		}else if(x[i]>ux-f*L){
			if(x[i]>ux) cout << "vdfplot: x>L: "<<x[i]<<" "<<v[i]<<endl;
			vhu[iv]+=1;
		}
		}
	}

	//cout << "vhl=";
	//for(int i=0; i<bins; ++i) cout << vhl[i] << " ";cout << endl;
	//cout << "vhu=";
	//for(int i=0; i<bins; ++i) cout << vhu[i] << " ";cout << endl;


	const int t=time;
	char snu[13]; sprintf(snu, "%d", t);
	string str="diag/vh";
	str=str+specie+"."+snu+".dat";
	
	gnuplotmark("vdf"+specie,str);

	double V;
	ofstream fvh(str.c_str(),ios::app);
	for(int i=0; i<bins; ++i) {
		V=(vu-vl)*(double(i)+0.5)/double(bins-1)+vl;
		fvh 	<<  V  << " " 
			<< vhl[i] << " " 
			<< (*fa)(0.,V) << " "
			<< (*frd)(0.,V) << " " 
			<< vhu[i] << " "
			<< (*fa)(L,V) << " "
			<< (*frd)(L,V) << " " 
			<< endl;
	}
	fvh.close();
}


double CIC(double x, double l, double u, double array[], int N){
	//Returns a linear interpolated value from array[], given at a
	//position(etc.) x in domain [l,u] (lower, upper). 
	//Size of array is N

	double w=(x-l)/(u-l)*(N-1);
	int m=int(floor(w));
	w-=m;

	if(w>1. || w<0.){
		cout<<" CIC(): w out of range: "<<w<<endl;
		cout<<"      "<<m<<" "<<l<<" "<<u<<" "<<x<<" "<<endl;
	}
	if(m<0 || m>=N-1){
		cout<<" CIC(): m out of range: "<<m<<endl;
		cout<<"      "<<m<<" "<<l<<" "<<u<<" "<<x<<" "<<endl;
	}
	
	return (1.-w)*array[m] + w*array[m+1];
}



void densityplot(double PHI[], int N, int time){
	int m;
	int bins=200;
	double bin=(ux-lx)/(bins-1);
	long ea[bins];
	long er[bins];
	long ia[bins];
	long ir[bins];

	for(int i=0; i<bins; ++i){
		ea[i]=0;
		er[i]=0;
		ia[i]=0;
		ir[i]=0;
	}

	double x, vv, vs2;
	double f=0.50;
	for(long i=0; i<Ne; ++i){
		x=xe[i];
		vv=ve[i]*ve[i];
		vs2=-2.*QMe*f*CIC(x,lx,ux,PHI,N);
		if(vv<vs2 && x>L/2.)	er[int((x-lx)/bin)]+=1;
		else	 		ea[int((x-lx)/bin)]+=1;
	}
	for(long i=0; i<Np; ++i){
		x=xp[i];
		vv=vp[i]*vp[i];
		vs2=2.*QMp*f*(phiDL-CIC(x,lx,ux,PHI,N));
		if(vv<vs2 && x<L/2.)	ir[int((x-lx)/bin)]+=1;
		else	 		ia[int((x-lx)/bin)]+=1;
	}

	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/density.";
	str=str+snu+".dat";
	gnuplotmark("density",str);
	
	ofstream fout(str.c_str());
	for(int i=0; i<bins-1; ++i)
		fout 	<<i*bin+lx<<" "
			<<ea[i]/bin<<" "
			<<er[i]/bin<<" "
			<<ia[i]/bin<<" "
			<<ir[i]/bin<<" "
			<<ea[i]/bin+er[i]/bin<<" "
			<<ia[i]/bin+ir[i]/bin<<" "
			<<endl;
	fout.close();
}



bool checknan(string type, double x, int i){
	if(isnan(x)){ 
		cout << type <<" isnan: " << i << " "<<x<<endl; 
		finalize(true);
		return true;
	}else 	return false;
}


void checkfornan(double RHO[], double PHI[], double E[], int N){
	cout << " - Checking for nans ";
	bool stop=false;
	for(int i=0; i<N; ++i){
		stop=checknan("rho",RHO[i],i);
		stop=checknan("phi",PHI[i],i);
		stop=checknan("E",E[i],i);
	}
	for(long i=0; i<Ne; ++i){ 
		if(isnan(xe[i])){ 
			cout<<"xe : "<<i<<" "<<xe[i]<<endl; break;}
		if(isnan(ve[i])){ 
			cout<<"ve : "<<i<<" "<<ve[i]<<endl; break;}
		if(isnan(dve[i])){ 
			cout<<"dve : "<<i<<" "<<dve[i]<<endl; break;}
	}
	for(long i=0; i<Np; ++i){
		if(isnan(xp[i])){ 
			cout<<"xp : "<<i<<" "<<xp[i]<< endl; break;}
		if(isnan(vp[i])){ 
			cout<<"vp : "<<i<<" "<<vp[i] << endl; break;}
		if(isnan(dvp[i])){ 
			cout<<"dvp : "<<i<<" "<<dvp[i]<<endl; break;}
	}

	if(stop){
		cout << "Forcing exit due to NaN's"<<endl; 
		finalize(true);
	}
}

















void print_simdetails(double vdtdx, int &N){
	//Stores total simulationtime as a system variable
	//char snu[13]; sprintf(snu, "%d", Ntime);
	//string str="export CRUNSIMTIME=";
	//str=str+snu;
	//system(str.c_str());

	double mach=abs(drifte/vthe);
	double dx=(ux-lx)/(N-1);
	ofstream params("diag/params.dat");
	params 	<<"  |  u="<<mach<<"   p="<<phiDL<<"   m="<<mp/me<<"   |  "
		<<"   t="<<N_time<<"    t*dt=";
	params 	<< setprecision(3) << N_time*dt << "{/Symbol w}_{pe}^{-1}"; 
	params	<<"   Nc="<<N<<"   |   "; 
	params	<< setprecision(2) 	<<" dt="<< dt <<" dx="<< dx
					<<" vdt/dx="<< vdtdx <<"   |  ";

	cout	<<endl;
	cout 	<<" Simulation parameters: "	<<endl
		<<"       x = ("<<lx<<","<<0<<","<<L<<","<<ux<<")"<<endl
		//<<"       L = "<<L		<<endl
		//<<"      lx = "<<lx		<<endl 
		//<<"      ux = "<<ux		<<endl 
		<<"      dx = "<<dx		<<endl
		<<"      dt = "<<dt		<<endl
		<<"  vdt/dx = "<<vdtdx		<<endl
		<<"    time = "<<N_time		<<endl
		<<" time*dt = "<<N_time*dt	<<endl
		<<endl;
	cout 	<<" Plasma and DL parameters: "	<<endl
		<<"   mi/me = "<<mp/me		<<endl
		<<"   Ti/Te = "<<vthp*vthp/me*mp<<endl
		<<"    mach = "<<mach		<<endl
		<<"   phiDL = "<<phiDL		<<endl
		<<endl;
#ifdef BC_vND
	cout<<" Boundary condition: vN - D"<<endl;
	params<<" vN-D";
#endif
#ifdef BC_DD
	cout<<" Boundary condition: D - D"<<endl;
	params<<" D-D";
#endif
#ifdef BC_vNvN
	cout<<" Boundary condition: vN - vN"<<endl;
	params<<" vN-vN";
#endif
#ifndef BC_DD
	#ifndef BC_vND
		#ifndef BC_vNvN
	cout<<" Boundary condition: D - vN"<<endl;
	params<<" D-vN";
		#endif
	#endif
#endif
	params.close();
#ifdef PTEST
	cout<<" Running test on Poisson solver... "<<endl;
#endif
}
















void initialstabilitycheck(double vdtdx, double mach, double TionTe, int &N){
	bool stop=false;

	cout << " Initial stabilitycheck... ";	

	//Test-something-something
//	if(vdtdx >= 0.6) cout << " dt too large! ", stop=true;
	
	//Bohm criterium
	if(mach<=(3.0+TionTe)) 
		cout<<" Bohm criterium violated! ";
	else	cout<<"";
//		cout<<" Bohm criterium OK!" <<endl;
	
	//Strong double layer
	if(phiDL<10.) 
		cout << " Weak double layer! ", stop=true;	
	
	//f(v=upper limit / lower limit)=0
	double zero=1e-10;
	double ea0=fea(0.,uea(0.));
	double eaL=fea(L,uea(L)) ;
	double ed0=fed(0.,led(0.));
	double edL=fed(L,led(L))  ;
	double ia0=fia(0.,lia(0.));
	double iaL=fia(L,lia(L))  ;
	double id0=fid(0.,uid(0.));
	double idL=fid(L,uid(L));
	if(ea0>zero) cout<<" fea(0,inf)=" <<ea0<<" ", stop=true;
	if(eaL>zero) cout<<" fea(L,inf)=" <<eaL<<" ", stop=true;
	if(ed0>zero) cout<<" fed(0,-inf)="<<ed0<<" ", stop=true;
	if(edL>zero) cout<<" fed(L,-inf)="<<edL<<" ", stop=true;
	if(ia0>zero) cout<<" fia(0,-inf)="<<ia0<<" ", stop=true;
	if(iaL>zero) cout<<" fia(L,-inf)="<<iaL<<" ", stop=true;
	if(id0>zero) cout<<" fid(0,inf)=" <<id0<<" ", stop=true;
	if(idL>zero) cout<<" fid(L,inf)=" <<idL<<" ", stop=true;
	
	//If important stabilitytest violated, stop simulations
	if(stop){
		cout	<<" forcing exit!"<<endl<<endl;
		exit(1);
	}else{
		cout << "  Finished" << endl;
	}
}
















void stabilitycheck(double phi[], double rho[], double E[], int N){
	
	cout << " Stabilitycheck:";
	//phi(x)>0
	for(int i=0; i<N; ++i)	
		if(phi[i]<0.){cout<<" phi(x)<0"; break;}
	
	//E<0
	for(int i=0; i<N; ++i)	
		if(E[i]>0.){ cout<<" E(x)>0"; break;}
	
	//rho
	double rhosum=0., dx=(ux-lx)/(N-1), m=0.0;
	for(int i=0; i<N; ++i){ 
		rhosum+=rho[i];			//Total charge density
		if(m>abs(rho[i])) m=abs(rho[i]);//Maximum value
	}
	double total=rhosum*dx/m*100.;		//Total rho to max in %
	double lb=abs(rho[0])*dx/m*100.;	//Lower boundary ---"---
	double ub=abs(rho[N-1])*dx/m*100.;	//Upper boundary ---"---
	double acceptance=5.;			//Unit percentage
	if(lb>acceptance) cout << " rho(0)>>0";
	if(ub>acceptance) cout << " rho(L)>>0";
	if(total>acceptance) cout << " sum(rho)>>0";

	cout << " ...Finished" << endl;
}






void printf(	double (*fa)(double, double), double (*la)(double), 
		double (*ua)(double), double (*fr)(double,double),
		double (*lr)(double), double (*ur)(double), 
		double (*fd)(double, double), double (*ld)(double), 
		double (*ud)(double), string specie){

	int Nx=200, Nv=200;

	double dx=(ux-lx)/(Nx-1);
	double vmax=0.0, vmin=0.0, tmp;
	for(double x=lx; x<=ux; x+=dx){
		tmp=(*ua)(x); if(vmax < tmp) vmax=tmp;
		tmp=(*ud)(x); if(vmax < tmp) vmax=tmp;
		tmp=(*la)(x); if(vmin > tmp) vmin=tmp;
		tmp=(*ld)(x); if(vmin > tmp) vmin=tmp;
	}
	double dv=(vmax-vmin)/(Nv-1);


	double x,v,f;
	string str="diag/printf.";
	str=str+specie+".dat";
	ofstream fout(str.c_str());
	
	for(double x=lx; x<=ux; x+=dx){
		for(double v=vmin; v<=vmax; v+=dv){
			if((*la)(x)<v && v<(*ua)(x)) 
					f=(*fa)(x,v);
			else if((*lr)(x)<=v && v<=(*ur)(x)) 
					f=(*fr)(x,v);
			else if((*ld)(x)<v && v<(*ud)(x)) 
					f=(*fd)(x,v);
			else		f=0.0;

			//fout<<x<<" "<<v<<" "<<f<<endl; 
			fout<<x<<" "<<v<<" "<<f<<endl; 
		}
		fout << endl;
	}
	fout.close();


}










void vdftest(){
	double dx=(ux-lx)/100;
	for(double x=lx; x<=ux; x+=dx){
		cout << fea(x,lea(x)) << " "
			<< fia(x,uia(x)) << " " 
		<< endl;
	}
}








