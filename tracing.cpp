const int Ntrace=20;
long traceindex[Ntrace];
double lastx[Ntrace], lastv[Ntrace];

void selecttrace(double x[], double v[], long N, double vs2){
	int m=0;
	long i;
	while(m<Ntrace){
		i=long(drand48()*(N-1));
		// If reflected particle, 
		// and ions at x approx lx
		// for m=0-m=9: v>0
		// for m=10-m=19: v<0
		if(v[i]*v[i]<vs2){
			if(m<Ntrace/2 && v[i]>0 && x[i]<(lx+L*0.1)){
				traceindex[m]=i;
				++m;
			}else if(m>=Ntrace/2 && v[i]<0 && x[i]>(lx+L*0.01) && x[i]<(lx+L*0.05)){
				traceindex[m]=i;
				++m;
			}
		} 

	}//End while
}



void trace(double x[], double v[], int time){
	ofstream fout("diag/trace.dat",ios::app);

	double X,V; 
	fout<<time<<" ";
	for(int m=0;m<Ntrace;++m){

		X=x[traceindex[m]];
		V=v[traceindex[m]];
		
		if(X<0){ 
			fout<<X<<" "<<V<<" ";
			lastx[m]=X; lastv[m]=V;
		}else{
			fout<<lastx[m]<<" "<<lastv[m]<<" ";
		}
	}
	fout<<endl;
		
	fout.close();
}
