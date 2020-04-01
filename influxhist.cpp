const int Nhf=200;
double influxhist[Nhf];


void reset_influxhist(){
	for(int i=0; i<Nhf; ++i) influxhist[i]=0.0;
}

void fill_influxhist(double v[], double lb, double ub, long N){
	
	double bin=(ub-lb)/(Nhf-1);
	int iv;
	for(long i=0; i<N; ++i){ 
		iv=int((v[i]-lb)/bin);
		++influxhist[iv];
	}
}

void print_influxhist(double lb, double ub, int time, string specie){

	
	double dv=(ub-lb)/(Nhf-1);
	double v;
	int i;

	ofstream fout("diag/influxhist.dat", ios::app);
	for(i=0, v=lb; i<Nhf; ++i,v+=dv) 
		fout<<time<<" "<<v<<" "<<influxhist[i]/dv<<endl;
	fout<<endl;


	fout.close();


	reset_influxhist();
}


