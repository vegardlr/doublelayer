
int phi_array_length=10000;
double *phi_array = new double[phi_array_length];
double phi(double x){	
	
	double w;
	int m;
	//Linear interpolation
	w=x*phi_array_length/L;
	m=int(floor(w));
	w-=m;
	
	if(m <= 0){ 
		return 0;
	}
	else if(m >= phi_array_length-1)	return phiDL;
	else	return (1.-w)*(phi_array[m]) + w*(phi_array[m+1]);
}




//ELECTRON HELP FUNCTIONS
double fea(double x, double y){
	double a=y*y+2.*QMe*phi(x);
	if(a*a<1e-8 && a < 0.0)a=0.0;
	double u=sqrt(a);
	
	//return 1./sqrt(2.*pi)/vthe*exp(-.5*(u-drifte)*(u-drifte)/vthe/vthe);
	return exp(-.5*(u-drifte)*(u-drifte)/vthe/vthe);
}
double vfea(double v){
	return v*fea(0., v);
}
double uea(double x){
	double veinf2=pow(15.*vthe+drifte,2);		//INFINTIY
	return sqrt(veinf2-2.*QMe*phiDL);
}
double lea(double x){
	return sqrt(-2.*QMe*phi(x));
}


double fer(double x, double y){
	//return 1./sqrt(2*pi)/vthe
	//	*exp( (-.5*y*y + QMe*(phiDL-phi(x)))/vthe/vthe  );
	return exp( (-.5*y*y + QMe*(phiDL-phi(x)))/vthe/vthe  );
}
double vfer(double v){
	return abs(v)*fer(L, v);
}
double uer(double x){
	return sqrt(-2.*QMe*phi(x));
}
double ler(double x){
	return -sqrt(-2.*QMe*phi(x));
}



double fed(double x, double y){
	//return fer(x,y);
	return fer(x,y);
}
double ued(double x){
	return -sqrt(-2.*QMe*phi(x));
}
double led(double x){
	double veinf2=pow(15.*vthe+drifte,2);		//INFINTIY
	return -sqrt(veinf2-2.*QMe*phiDL);
}




//	ION HELP FUNCTIONS
double fia(double x, double y){

	double a=y*y-2.*QMp*(phiDL-phi(x));
	if(a*a<1e-8 && a < 0.0)	a=0.0;
	double u=-sqrt(a);
	//return 1./sqrt(2*pi)/vthp
	//	*exp(  -.5*(u-driftp)*(u-driftp)/vthp/vthp  );
	return exp(  -.5*(u-driftp)*(u-driftp)/vthp/vthp  );
}
double vfia(double v){
	return abs(v)*fia(L, v);
}
double uia(double x){
	return -sqrt(2.*QMp*(phiDL-phi(x)));
}
double lia(double x){
	double vpinf2=pow(15.*vthp-driftp,2);		//INFINITY
	return -sqrt(vpinf2+2.*QMp*phiDL);
}



double fir(double x, double y){	
	//return 1./sqrt(2*pi)/vthp
	//	*exp((-.5*y*y - QMp*phi(x))/vthp/vthp  );
	return exp((-.5*y*y - QMp*phi(x))/vthp/vthp  );
}
double vfir(
	double v){return v*fir(0., v);
}
double uir(double x){
	return sqrt(2.*QMp*(phiDL-phi(x)));
}
double lir(double x){
	return -sqrt(2.*QMp*(phiDL-phi(x)));
}



double fid(double x, double y){
	//return fir(x,y);
	return fir(x,y);
}
double uid(double x){
	double vpinf2=pow(15.*vthp+driftp,2);		//INFINTIY
	return sqrt(vpinf2+2.*QMp*phiDL);
}
double lid(double x){
	return sqrt(2.*QMp*(phiDL-phi(x)));
}

























//INTEGRALS OF DENSITY
double nea(double phi){

	int it=0;

	double a=-2*phi*QMe;
	double dv=0.05*vthe;
	double n=0.0;
	double d=0.0;
	double v=dv;

	//ofstream fil("diag/nea.dat");
	//fil << n << " " << d << " " << v << endl;
	do{
		d=v*fea(0.,v)/sqrt(v*v + a);
		n+=d*dv;
		v+=dv;
		++it;

		if(isnan(d)){
			cout << "nea isnan: ";
			cout 	<< fea(0.,v) << " "
				<< v*v+a << " " 
				<< v*v << " " 
				<< a << " " 
				<< endl;
		}
		//fil << n << " " << d << " " << v << endl;
		
	}while(d>0.01*dv || it < 100);
	//fil.close();

//	cout << "   nea it=" << it << " n=" << n <<endl;
	return n;
}


double ned(double phi){	
	
	int it=0;

	double v=sqrt(-2.*phiDL*QMe);
	double vs2=2*(phi-phiDL)*QMe;
	double dv=0.01*v;
	double d=0;
	double n=0.0;
	

	n=(3./2.)*sqrt(v*dv/2)*fed(L,-v);
	v+=dv;





	//ofstream fil("diag/ned.dat");
	//fil << n << " " << d << " " << v << endl;
	do{
		d=v*fed(L,-v)/sqrt(v*v - vs2);
		n+=d*dv;
		v+=dv;
		++it;

		if(isnan(d)){
			cout << "ned isnan ";
			cout 	<< fed(L,-v) << " "
				<< v << " " << drifte << " " << vthe << " " 
				<< v*v-vs2 << " " 
				<< v*v << " " 
				<< vs2 << " " 
				<< endl;
		}
	//fil << n << " " << d << " " << v << endl;

	}while(d>0.01*dv || it < 100);
	//fil.close();

//	cout << "   ned it=" << it << " n=" << n <<endl;
	return n;

}

double ner(double phi){
	double fer(double, double);


	int it=0;

	double vs2=2*(phi-phiDL)*QMe;
	double vs0=sqrt(-2*phiDL*QMe);
	double v=sqrt(abs(vs2));
	double dv=0.01*(vs0-v);
	double d=0.0;
	double n=0.0;


	if(v==0.) n=fer(L,-v)*dv;
	else{ 
		n=3.*sqrt(v*dv/2.)*fer(L,-v);
		if(dv<1e-12) return n;
	}
	
	v+=dv;



	if(isnan(n)) 
		cout 	<< "ner (1) gives nan, phi=" << phi
		<< " dv=" << dv << endl;

	
	//ofstream fil("diag/ner.dat");
	//fil << n << " " << d*dv << " " << v << endl;	
	while(v<vs0){
		d=2.*v*fer(L,-v)/sqrt(v*v - vs2);
		n+=d*dv;
		v+=dv;
		++it;

		if(isnan(n)){
			cout << "ner isnan ";
			cout 	<< fer(L,-v) << " "
				<< v*v-vs2 << " " 
				<< v*v << " " 
				<< vs2 << " " 
				<< dv << " "
				<< phi << " "
				<< endl;
		}
	
	//fil << n << " " << d*dv << " " << v << endl;
	}
	if(isnan(n)) 
		cout << "ner (2) gives nan, phi=" << phi
		<< " dv=" << dv << endl;
	
	n+=v*fer(L,-v)*dv/sqrt(v*v-vs2);
	
	//fil << n << " " << d*dv << " " << v << endl;
	//fil.close();

	if(isnan(n)) 
		cout 	<< "ner (3) gives nan, phi=" 
			<< phi<< " dv=" << dv << endl;
//	if(n > 100.) 
//		cout 	<< "ERROR r=" << d*dv/n 
//			<< "     n="<< n << " d="<< d*dv 
//			<< " dv="<< dv << endl;
//	cout << "   ner it=" << it << " n=" << n <<endl;
	return n;

}






double nia(double phi){
	double fia(double, double);

	int it=0;

	double vs2=2*(phiDL-phi)*QMp;
	double dv=0.05*vthp;
	double n=0.0;
	double d=0.0;
	double v=0.0;
	
	//ofstream fil("diag/nia.dat");
	//fil << n << " " << d << " " << v << endl;
	do{
		v+=dv;
		d=v*fia(L,-v)/sqrt(v*v + vs2);
		n+=d*dv;
		++it;

		if(isnan(d)){
			cout << "nia isnan ";
			cout 	<< fia(L,v) << " "
				<< v*v+vs2 << " " 
				<< v*v << " " 
				<< vs2 << " " 
				<< endl;
		}
	//fil << n << " " << d << " " << v << endl;
		
	}while(d>0.01*dv || it < 100);
	//fil.close();

//	cout << "   nia it=" << it <<" n=" << n << endl;
	return n;
}


double nid(double phi){	
	double fid(double,double);
	
	int it=0;

	double a=2*phi*QMp;
	double v=sqrt(2.*phiDL*QMp);
	double dv=0.01*v;
	double d=0.0;
	double n=0.0;
	
	n=(3./2.)*sqrt(v*dv/2.)*fid(0.,v);
	v+=dv;

	//ofstream fil("diag/nid.dat");
	//fil << n << " " << d << " " << v << endl;
	do{
		d=v*fid(0.,v)/sqrt(v*v - a);
		n+=d*dv;
		v+=dv;
		++it;

		if(isnan(d)){
			cout << "nid isnan ";
			cout 	<< d << " " << n << " "
				<< fid(0,v) << " "
				<< v*v-a << " " 
				<< v*v << " " 
				<< a << " " 
				<< endl;
		}
		//fil << n << " " << d << " " << v << endl;
	}while(d>0.1*dv || it < 100);
	//fil.close();

//	cout << "   nid it=" << it << " n=" << n <<endl;
	return n;

}

double nir(double phi){
	double fir(double,double);
	

	int it=0;

	double a=2.*phi*QMp;
	double vs0=sqrt(2*phiDL*QMp);
	double v=sqrt(2.*(phi)*QMp);
	double dv=0.01*(vs0-v);
	double d=0.0;
	double n=0.0; 
	

	if(v==0.) n=fir(0.,v)*dv;
	else{ 
		n=3.*sqrt(v*dv/2.)*fir(0.,v);
		if(dv<1e-12) return n;
	}


	v+=dv;



	if(isnan(n)) 
		cout 	<< "nir (1) gives nan, phi=" 
			<< phi<< " dv=" << dv << endl;

	//ofstream fil("diag/nir.dat");
	//fil << n << " " << d*dv << " " << v << endl;
	while(v<vs0){
		d=2.*v*fir(0.,v)/sqrt(v*v - a);
		n+=d*dv;
		v+=dv;
		++it;

		if(isnan(n))
			cout << "nir isnan  "
			 	<< fir(0.,v) << " "
				<< v*v-a << " " 
				<< v*v << " " 
				<< a << " " 
				<< dv << " "
				<< phi << " " 
				<< endl;
		
	//fil << n << " " << d*dv << " " << v << endl;
	}
	
	if(isnan(n)) cout 	<< "nir (2) gives nan, phi=" 
				<< phi<< " dv=" << dv << endl;
	n+=v*fir(0.,v)*dv/sqrt(v*v-a);
	
	//fil.close();
	
	if(isnan(n)) cout 	<< "nir (3) gives nan, phi=" 
				<< phi<< " dv=" << dv << endl;
//	cout << "   nir it=" << it << " n="<< n << endl;
	return n;

}
