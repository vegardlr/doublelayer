	/*
	 *	INITIALIZE FLUX
	 *
	 */
double NRintegration(double (*f)(double),double l, double u, int steps){
	//Newton Raphsons integration scheme

	double dx=(u-l)/(steps-1);
	double I=.5*(*f)(l)*dx;
	for(double x=l; x<u; x+=dx) I += (*f)(x)*dx;
	I+=.5*(*f)(u)*dx;
	return I;
}

void initflux(double alpha[]){
	//Finds the expected number of particles to enter at each
	//timestep. Experiencing a flux in 100 times lower than 
	//what leaves the domain. 
	//This is independent of what I choose L and dt to be.

	flux_ea  = alpha[0]*dt*NRintegration(vfea, 0., uea(lx), 200);
	flux_erd = alpha[1]*dt*NRintegration(vfer, led(ux), 0., 200);
	flux_ia  = alpha[2]*dt*NRintegration(vfia, lia(ux), 0., 200);
	flux_ird = alpha[3]*dt*NRintegration(vfir, 0., uid(lx), 200);

	
	cout << " Influx " << endl;
	cout << "   Pea="<<flux_ea<<" Perd="<<flux_erd<<" Pia="
		<<flux_ia<<" Pird="<<flux_ird<<endl<<endl;

}


























	/*
	 *	INFLUX OF PARTICLES
	 */
long fluxnumber(double P){
	//From a double telling the flux, return a long containing total 
	//number of particles to be drawn. Can handle non integer P's.
	
	if(drand48()<=(P-floor(P))) return long(P)+1;
	else return long(P);
}

double max(double x1, double x2, double x3, double x4){
	//Return the maximum of the 4 double variables.
	
	double m=0.;
	if(x1>m)m=x1;
	if(x2>m)m=x2;
	if(x3>m)m=x3;
	if(x4>m)m=x4;
	return m;
}
	
void position(long N, double v[], double x[], double x0){
	//From known velocities, find their positions	
	double drand;
	for(long i=0; i<N; ++i){ 
		x[i] = x0 + v[i]*dt*(drand=drand48());
		if(x[i]>ux || x[i]<lx){
			cout
			<<"    flux.cpp: position(): x out of bounds (x=" 
			<< x[i] << ") " << v[i] << " " << drand << endl;
		}
	}
}




void influx(int time){
	//Generate particles according to their distrobution functions 
	//vfea, etc. Use the Inversion method supplied in Utils.cpp.

	void histogram(int time, string specie, double x[], double v[], long Nxv, 
			double (*f)(double), double, double);
	double *x, *v, fnorm;
	long N;

	double m=max(flux_ea, flux_erd, flux_ia, flux_ird);
	x = new double[long(m+1)];
	v = new double[long(m+1)];
	
//	cout << "influx()"<<endl;

	N=fluxnumber(flux_ea);
	Inversion ea(0., uea(lx), 200, vfea, fnorm);
	ea.draw(N,v,false);
#ifndef TVSELECT
	inTVhist(time, v, vine, tine, Ne, N);
#endif
	position(N,v,x,lx);
	fill(&Ne,xe,ve,N,x,v);
//	histogram(time, "ea", x, v, N, vfea, 0., uea(lx));

	N=fluxnumber(flux_erd);
	Inversion erd(led(ux), 0., 200, vfer, fnorm);
	erd.draw(N,v,false);
#ifndef TVSELECT
	inTVhist(time, v, vine, tine, Ne, N);
#endif
	position(N,v,x,ux);
	fill(&Ne,xe,ve,N,x,v);
//	histogram(time, "erd",x, v, N, vfer, led(ux), 0.);

	N=fluxnumber(flux_ia);
	Inversion ia(lia(ux), 0., 200, vfia, fnorm);
	ia.draw(N,v,false);
#ifndef TVSELECT
	inTVhist(time, v, vinp, tinp, Np, N);
#endif
	position(N,v,x,ux);
	fill(&Np,xp,vp,N,x,v);
//	histogram(time,"ia",x, v, N, vfia, lia(L), 0.);

	N=fluxnumber(flux_ird); 
	Inversion ird(0., uid(lx), 200, vfir, fnorm);
	ird.draw(N,v,false);
#ifndef TVSELECT
	inTVhist(time, v, vinp, tinp, Np, N);
#endif
	position(N,v,x,lx);
	fill(&Np,xp,vp,N,x,v);
//	histogram(time,"ird",x, v, N, vfir, 0., uid(0.));
//	fill_influxhist(v,0.,uid(ux),N);
	

	
	delete x,v;
}



















	/*
	 *	OUTFLUX OF PARTICLES
	 */
void outflux(int time){
	for(int i=0;i<Np;++i) while(xp[i]>ux||xp[i]<lx){
		--Np;
		xp[i]=xp[Np];
		vp[i]=vp[Np];
	};
	for(int i=0;i<Ne;++i) while(xe[i]>ux||xe[i]<lx){
		--Ne;
		xe[i]=xe[Ne];
		ve[i]=ve[Ne];
	};
}


void diagoutflux(int time){
	double f=0.5;
	double vse2=-2.*QMe*phiDL*f;
	double vsp2=2.*QMp*phiDL*f;
	int out_e0a=0, out_e0rd=0, out_eLa=0, out_eLrd=0;
	int out_p0a=0, out_p0rd=0, out_pLa=0, out_pLrd=0;



	//OUT
	for(int i=0;i<Np;++i) while(xp[i]>ux||xp[i]<lx){
#ifndef TVSELECT	
		outTVhist(time,i,Np,vinp,tinp,htvp,lvp,uvp);
#endif
		if(xp[i]>ux){
			if(vp[i]*vp[i]>vsp2){
				out_pLa+=1;	
			}else{
				out_pLrd+=1;
			}	
		}else{
			if(vp[i]*vp[i]>vsp2){
				out_p0a+=1;
			}else{
				out_p0rd+=1;
			}	
		}
		--Np;
		xp[i]=xp[Np];
		vp[i]=vp[Np];
		tinp[i]=tinp[Np];
		vinp[i]=vinp[Np];

	};
	for(int i=0;i<Ne;++i) while(xe[i]>ux||xe[i]<lx){
#ifndef TVSELECT
		outTVhist(time,i,Ne,vine,tine,htve,lve,uve);
#endif		
		if(xe[i]>ux){
			if(ve[i]*ve[i]>vse2){
				out_eLa+=1;	
			}else{
				out_eLrd+=1;
			}	
		}else{
			if(ve[i]*ve[i]>vsp2){
				out_e0a+=1;
			}else{
				out_e0rd+=1;
			}	
		}
		--Ne;
		xe[i]=xe[Ne];
		ve[i]=ve[Ne];
		tine[i]=tine[Ne];
		vine[i]=vine[Ne];
	};
	
	ofstream fout("diag/flux2.dat",ios::app);
	fout	<< time*dt << " "
		<< flux_ea << " " 
		<< out_e0a << " " 
		<< out_eLa << " "
		<< flux_erd << " "
		<< out_e0rd << " "
		<< out_eLrd << " "
		<< flux_ia << " " 
		<< out_p0a << " " 
		<< out_pLa << " "
		<< flux_ird << " "
		<< out_p0rd << " "
		<< out_pLrd << " "
		<< " - " 
		<< (double)out_eLa/flux_ea << " "
		<< (double)out_eLrd/flux_erd << " "
		<< (double)out_p0a/flux_ia << " "
		<< (double)out_p0rd/flux_ird << " "
		<< endl;
	fout.close();

	
}
















	/*
	 *	FLUX "MAIN" SCHEME
	 */
void flux(int time){
	//printTVhist(time);
	int oute=Ne, outp=Np;	//Count particles that leave/enter domain
	//outflux(time);	//Remove particles that left the domain
	diagoutflux(time);
	oute-=Ne; outp-=Np;;	//Count particles that leave/enter domain
	int ine=Ne, inp=Np;;	//Count particles that leave/enter domain
	influx(time);		//Insert new particles
	ine=Ne-ine; inp=Np-inp;;//Count particles that leave/enter domain

	//Write in-/outflux to file
	ofstream fout("diag/flux.dat",ios::app);
	fout	<< time*dt << " "
		<< Ne << " "
		<< Np << " "
		<< ine << " "
		<< oute << " "
		<< inp << " "
		<< outp << " " 
		<< -ine+oute-outp+inp << " "
		<< endl;
	fout.close();
}






















	/*
	 *	DIAGNOSTICS
	 */ 
/*
void testfluxnumber(){
	long ea=0, erd=0, ia=0, ird=0;
	int N=9000;
	for(int i=0; i<N; ++i){
		ea+=fluxnumber(flux_ea);
		erd+=fluxnumber(flux_erd);
		ia+=fluxnumber(flux_ia);
		ird+=fluxnumber(flux_ird);
	}
	cout<<flux_ea<<" "<<flux_erd<<" "<<flux_ia<<" "<<flux_ird<<endl;
	cout	<<double(ea)/double(N)<<" "
		<<double(erd)/double(N)<<" "
		<<double(ia)/double(N)<<" "
		<<double(ird)/double(N)<<" "
	<<endl;
}


void testflux(){
	//initflux(alpha);
	//testfluxnumber();
	double high=2000.;
	flux_ea=high; 
	flux_erd=high;
	flux_ia=high;
	flux_ird=high;
	influx(0);
}
*/

/*
void histogram(int time, string specie, double x[], double v[], long Nxv, double (*f)(double), double vl, double vu){
	
	//Position histogram initializations
	int xbins=200;
	double xbin=(ux-lx)/(xbins-1);
	double xhist[xbins];
	for(int m=0; m<xbins; ++m)	xhist[m]=0.;
	
	//Velocity histogram initializations
	int vbins=xbins;
	double vbin=(vu - vl)/(vbins-1);
	double vhist[vbins];
	for(int m=0; m<vbins; ++m)	vhist[m]=0.;
	
	//Indexvalues
	int ix,iv;
	//Filling histograms
	for(long i=0;i<Nxv;++i){
		ix=int((x[i]-lx)/xbin);
		iv=int( (v[i]-vl)/vbin+0.5);
		if(ix<0 || iv<0 || ix>xbins || iv>vbins)
			cout<<"  HI ERROR!!  "<<ix<<" "<<iv<<"   X   ";
		xhist[ix]+=1.;
		vhist[iv]+=1.;
	}

	//Forcing normalizing of histograms and functions.
	double vnorm=0., k;
	for(int i=0; i<xbins; ++i){
		k=(vu-vl)*double(i)/double(vbins-1)+vl;
	
		vnorm+=(*f)(k);
	}
	vnorm=double(Nxv)/vnorm;

	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/influxhist.";
	str=str+specie+"."+snu+".dat";
	ofstream fout(str.c_str());

	double V;
	for(int i=0; i<xbins; ++i){ 
		V=(vu-vl)*double(i)/double(vbins-1)+vl;
		fout << i*xbin+lx << " "
				<< xhist[i] << " "
				<< V << " " 
				<< vhist[i] << " "
				<< vnorm*(*f)(V) << " " 
				<< endl; 
	}
	fout.close();
//	cout << "Flux histogram written, do plots." << endl;
//	double t; cout << "Enter: "; cin >> t;

}
*/
