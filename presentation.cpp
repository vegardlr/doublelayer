int Pdiag=0;


void Pdiagagnostics(double RHO[], double PHI[], double E[], int N, 
		int time){		
	
	void Pphidlprintout(int,double);
	void Pfieldcarpet(int,double [], double [], double [], int);


	if(time*20/N_time == Pdiag){
		cout<<" Presentation data at t="<<time<<endl;
		Pphidlprintout(time, (PHI[N-1]-PHI[0]) );	
		fields(RHO, PHI, E, N, time);
		Pfieldcarpet(time,RHO,PHI,E,N);
		ps_plot("xe-ve",xe,lx,ux,ve,led(ux),uea(ux),Ne,time);
		ps_plot("xp-vp",xp,lx,ux,vp,lia(lx),uid(lx),Np,time);

		++Pdiag;
	}

}




void Pphidlprintout(int time, double dl){
	//phiDL printout
	ofstream pout("diag/phidl.dat",ios::app);
	pout<<time*dt<<" "<<dl<<endl;
	pout.close();
	if(phiDL<0) cout<<endl<<" ERROR: phiDL<0"<<endl, exit(0);

	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/phidl.";
	str = str+snu+".dat";
	str = "cp diag/phidl.dat "+str;
	system(str.c_str());


}

void Pfieldcarpet(int time, double RHO[], double PHI[], 
			double E[], int N){
	ofstream carpet("diag/carpet.dat", ios::app);

	double x, dx=(ux-lx)/(N-1);
	int i; 
	for(i=0, x=lx; i<N; ++i,x+=dx)
		carpet<<time*dt<<" "<<x<<" "
			<<RHO[i]<<" "<<PHI[i]<<" "<<E[i]<<endl;
	
	carpet<<endl;
	carpet.close();

	char snu[13]; sprintf(snu, "%d", time);
	string str="diag/carpet.";
	str = str+snu+".dat";
	str = "cp diag/carpet.dat "+str;
	system(str.c_str());
}

