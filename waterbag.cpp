void tanhphi(){
	cout << " Doing the phi(x)~tanh(x) thing! "<<endl;

	//Set phi
	int N=200;
	double phidl=phiDL, length=100, width=10;
	double x, dx=length/(N-1);
	double phi[N];
	int i;
	for(i=0, x=0;i<N; ++i, x+=dx)
		phi[i]=(tanh((x-length/2.)/width)+1.)/2.*phidl;
	
	
	//Set densities
	double ea[N], er[N], ed[N], ia[N], ir[N], id[N];
	for(i=0,x=0;i<N;++i,x+=dx){
		ea[i]=nea(phi[i]);
		er[i]=ner(phi[i]);
		ed[i]=ned(phi[i]);
		ia[i]=nia(phi[i]);
		ir[i]=nir(phi[i]);
		id[i]=nid(phi[i]);
	}

	ofstream fout("diag/tanhphi.dat");
	for(i=0, x=0;i<N; ++i, x+=dx)
		fout	<<x<<" "<<phi[i]<<" "
			<<ea[i]<<" "
			<<er[i]<<" "
			<<ed[i]<<" "
			<<ia[i]<<" "
			<<ir[i]<<" "
			<<id[i]<<" "
			<<endl;
	fout.close();

	cout<< "tanh phi tests done, exiting."<<endl;
	exit(0);
}



//Waterbag
const double w1e=0.2, w4e=2.0;
const double w1i=2.0, w4i=0.2;
const double wbf0=1.0;


double wbuea(double x){	return sqrt(w4e-2.*QMe*phi(x));	}
double wblea(double x){ return sqrt(-2.*QMe*phi(x));	}
double wbfea(double x, double y){
	//if(y>wblea(x) && y<=wbuea(x)) return wbf0;
	//else return 0.;
	return wbf0;
}

double wbuer(double x){	return sqrt(-2.*QMe*phi(x));	}
double wbler(double x){	return -sqrt(-2.*QMe*phi(x));	}	
double wbfer(double x, double y){
	//if(y>=wbler(x) && y<=wbuer(x)) return wbf0;
	//else return 0.;
	return wbf0;
}

double wbued(double x){	return -sqrt(-2.*QMe*phi(x));	}
double wbled(double x){	return -sqrt(w1e-2.*QMe*phiDL);	}
double wbfed(double x, double y){
	//if(y>=wbled(x) && y<wbued(x)) return wbf0;
	//else return 0.;
	return wbf0;
}


//	ION HELP FUNCTIONS
double wbuia(double x){	return -sqrt(2.*QMp*(phiDL-phi(x)));	}
double wblia(double x){	return -sqrt(w1i+2.*QMp*(phiDL-phi(x)));}
double wbfia(double x, double y){
	//if(y>=wblia(x) && y<wbuia(x)) return wbf0;
	//else return 0.;
	return wbf0;
}

double wbuir(double x){	return sqrt(2.*QMp*(phiDL-phi(x)));	}
double wblir(double x){	return -sqrt(2.*QMp*(phiDL-phi(x)));	}
double wbfir(double x, double y){	
	//if(y>=wblir(x) && y<=wbuir(x)) return wbf0;
	//else return 0.;
	return wbf0;
}

double wbuid(double x){	return sqrt(w1i+2.*QMp*phi(x));		}
double wblid(double x){	return sqrt(2.*QMp*(phiDL-phi(x)));	}
double wbfid(double x, double y){
	//if(y>wblid(x) && y<=wbuid(x)) return wbf0;
	//else return 0.;
	return wbf0;
}

void waterbag(){
	void printf(double (*fa)(double, double), double (*la)(double), 
		double (*ua)(double), double (*fr)(double,double),
		double (*lr)(double), double (*ur)(double), 
		double (*fd)(double, double), double (*ld)(double), 
		double (*ud)(double), string specie);
	cout << "   Electrons"<<endl;
	printf(wbfea, wblea, wbuea, wbfer, wbler, wbuer, wbfed, wbled, wbued, "e");
	cout << "    Ions"<<endl;
	printf(wbfia, wblia, wbuia, wbfir, wblir, wbuir, wbfid, wblid, wbuid, "i");
}

