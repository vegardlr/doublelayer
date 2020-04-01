//#define index(x,y,c)

//int *hist[dim];
//for(i=0;i<dim;++i) new hist[i] = new int[dim2];
//hist[j][k];

const int histtime=N_time; 
int TVhistCount=1;
const int binst=200, binsv=200;
long TVhistNe, TVhistNp;
double *vine, *vinp;		//In-velocity
int *tine, *tinp;		//In-time
double lve, uve, lvp, uvp; 	//Upper/lower electron/proton velocity
double htve[binsv][binst];	//Histograms
double htvp[binsv][binst];	//Histograms
double vzoom=5.;

#ifdef TVSELECT
	double xlowlimit,xuplimit;
#endif

//Signed power of 2 of x
double spow2(double x){
	return x*x*x/abs(x);
}


void initTVhist(long ne, long np){
	
	double uea(double);
	double ued(double);
	double uia(double);
	double uid(double);

#ifdef TVSELECT
	xlowlimit=0.,xuplimit=L;
#endif
	TVhistNe=long(ne*1.5);
	TVhistNp=long(np*1.5);

	uve=spow2(uea(L));
	lve=spow2(led(L));
	uvp=spow2(uid(0.));
	lvp=spow2(lia(0.));
	
	vine = new double[TVhistNe];	
	tine = new int[TVhistNe];	
	vinp = new double[TVhistNp];	
	tinp = new int[TVhistNp];	
}


void resetTVhist(){
		for(long i=0; i<TVhistNe; ++i)	tine[i]=0;
		for(long i=0; i<TVhistNp; ++i)	tinp[i]=0;
	
		for(int v=0;v<binsv;++v)
			for(int t=0;t<binst;++t){
				htve[v][t]=0.;
				htvp[v][t]=0.;
		}

}

void inTVhist(int time, double v[], double vin[], 
			int tin[], long old, long gen){
	long i,ii;
	for(i=0, ii=old; i<gen; ++i, ++ii){
		vin[ii]=spow2(v[i]);
		tin[ii]=time;
	}
}

void outTVhist(int time, long index, long N, double vin[], int tin[], 
	double hist[binsv][binst], double lv, double uv){
	int iv, it;

	//Zoom in on velocity scale, that is look only at the inner 1/f
	double l=lv/vzoom, u=uv/vzoom;
	
	if(tin[index]>0){
		it=int(double(time-tin[index])/double(N_time)*(binst-1));
		iv=int((vin[index]-l)/(u-l)*(binsv-1));

		if(iv>0 && iv<(binsv-1) && it>0 && it<(binst-1)){
			++hist[iv][it];
		}
	}
//	THE FOLLOWING LINES ARE PLACED IN FLUX.CPP
//	tin[index]=tin[N-1];
//	vin[index]=vin[N-1];
}


void printTVhist(int time){
//	if(time/TVhistCount==histtime){
		cout << " Printing 'transit time'-'entry velocity'"
			<<" histogram at t="<<time;
	
		char snu[13]; sprintf(snu, "%d", time);
		string str="diag/TVhist.";
		string stre=str+"e."+snu+".dat";
		string strp=str+"p."+snu+".dat";
		ofstream eout(stre.c_str());
		ofstream pout(strp.c_str());

		void gnuplotmark(string,string);
		gnuplotmark("TVhiste",stre);
		gnuplotmark("TVhistp",strp);

		double Ve, Vp, T;
		for(int v=0; v<binsv; ++v){
			Ve=(double(v)*(uve-lve)/double(binsv-1)+lve)*vzoom;
			Vp=(double(v)*(uvp-lvp)/double(binsv-1)+lvp)*vzoom;
			for(int t=0;t<binst;++t){
				T=t*dt;
				eout<<spow2(Ve)<<" "<<T<<" "<<htve[v][t]<<endl;
				pout<<spow2(Vp)<<" "<<T<<" "<<htvp[v][t]<<endl;
			}	
			eout<<endl;	pout<<endl;	
		}
		eout.close();	pout.close();
	
	//	++TVhistCount;
	//	resetTVhist();
		cout << " ...Finished" << endl;
//	}
}
