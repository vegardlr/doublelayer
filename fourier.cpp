#ifdef FFTW
	#include "fftw.cpp"
#else
	#include "FFT.cpp"
#endif
void powerspectrum(double data[], int N){
#ifdef FFTW	
	void fftw(double [], double [],int,int);
	
	//Converting real data to complex array
	double re[N], im[N];
	for(int i=0; i<N; ++i){
		re[i]=data[i];	//Real part
		im[i]=0.;	//Imaginary part
	}
	//Fourier transform of complex array
	fftw(re,im,N,-1); //-1 : forward transform
	
	//Calculating powerspectrum
	for(int i=0; i<N; ++i) data[i]=re[i]*re[i] + im[i]*im[i];
#else
	FFT fft;
	double W[N];
	fft.forward(W,data,N);
	for(int i=0; i<N/2; ++i) data[i]=W[i]*W[i] + W[N-i-1]*W[N-i-1];
#endif
	
	//Checking for errors in results
	for(int i=0; i<N; ++i) if(isnan(data[i])) 
		cout<<" ps is NaN: "<<i<<endl; 

#ifdef PSTEST	
	//Frequency
	int N21=N/2-1;
	double f[N], df=1.;

	for(int i=0; i<=N21; ++i){
		f[i]=i*df;
		f[N-i-1]=-(i+1)*df;
	}

	ofstream pow("diag/ps.dat");
	for(int i=0; i<N; ++i){
		pow << f[i]<< " "<<data[i]<<endl;
	}
	pow.close();
#endif
}


int findmodes(double data[], int N, int modes[], 
		int modessize, string type){
	static int count=0;
	//Frequency
	int N21=N/2-1;
	double f[N], df=1.;

	for(int i=0; i<=N21; ++i){
		f[i]=i*df;
		f[N-i-1]=-(i+1)*df;
	}

	N/=2;		//Last half of powerspectrum is negative 
			//frequencies. Not interesting!

	int finds=0;		//Count no of finds 
	double mean=0.0;	//Find mean
	double tolerance=8.;	//Set tolerance
	for(int i=0; i<N; ++i)	mean+=data[i];	mean/=N;
	tolerance*=mean;
	
	//Finding derivative of powerspectrum
	double der[N-1];
	for(int i=0; i<N-1; ++i)der[i]=data[i+1]-data[i];

	//If local maxima > tolerance >> mean value = peak found
	int j;
	int k;
	for(int i=1; i<N-1; ++i)
		if((der[i]<0 && der[i-1]>0) && data[i] > tolerance){
			j=0;
			//Run through stored modes
			//If i matches j'th element of modes, i is
			//allready found, go to next found i.
			//If i doens't match, j'th element of modes is 
			//either empty or occupied by another mode.
			//If either of these, enter while loop.
			while(modes[j] != i && j<modessize){

				//If m[j] is empty, store i in m[j]
				if(modes[j] == -1){ 
					modes[j]=i;
					++finds;			
				}
				//If m[j] isn't empty, check m[j+1] 
				//(by exiting this loop with
				//incrementation of j
				else{	
					++j;
				}
			}

			k=0;
			if(modes[k]!=-1)
				cout<<"  Found powerspectrum modes, "
						+type+": ";
			while(modes[k]!=-1){ 
				cout << "m["<<k<<"]="<< f[modes[k]]<<" ";
				++k;
			}
			cout << endl;

		}
		
#ifdef MOTHERDIAG	
	//For diagnostics, print mean and tolerance
	char snu[13]; sprintf(snu, "%d", count);
	string str="diag/findmodes.";
	str=str+type+"."+snu+".dat";
	ofstream m(str.c_str());
	for(int i=0; i<N; ++i) m<<f[i]<<" "<<data[i]<<" "
				<<mean<<" "<<tolerance<<endl;
	m.close();
	++count;
#endif
	return finds;
		
}

void windowing(double data[], int N){
	double x, f, dx=10./N; 
	int i;

	for(i=0, x=0.; i<N/2; ++i, x+=dx){
		f=tanh(x);
		data[i]*=f;
		data[N-1-i]*=f;
	}
		
}


void writemodes(double data[], int N, 
		int modes[], double time, string type, double dx){
/*	
	//Frequency
	int N21=N/2-1;
	double f[N], df=1.;

	for(int i=0; i<=N21; ++i){
		f[i]=i*df;
		f[N-i-1]=-(i+1)*df;
	}
*/
	string tagfile="diag/powertags."+type+".dat";
	string powfile="diag/power."+type+".dat";
	string gnuplotfile="diag/power."+type+".p";
	ofstream tag(tagfile.c_str());
	ofstream pow(powfile.c_str(), ios::app);
	ofstream gnu(gnuplotfile.c_str());
	
	//If any modes found at all
	if(modes[0]!=-1){
		tag<<"time ";
		pow<<time<<" ";
		gnu<<"data='../"<<tagfile<<"'"<<endl<<"plot ";
	}

	
	double power, wavenumber;
	int j=0;
	while(modes[j]!=-1){
		power=data[modes[j]];
	//	wavenumber=f[modes[j]];
		wavenumber=2.*pi*modes[j]/(dx*N);
		tag<< setprecision(4) <<wavenumber<<" ";
		pow<<power<<" ";
		gnu	<<" data using 1:"<<j+2
			<<" title "<<j+2<<" with lines ";
		
		if(modes[j+1]!=-1)gnu<<",";
		else gnu<<endl;
		++j;
	}
	pow<<endl;
}


void finalizemother(){
	ofstream left("diag/powertags.left.dat", ios::app);
	ofstream right("diag/powertags.right.dat", ios::app);
	
	left<<endl;	right<<endl;
	left.close();	right.close();
	
	string cmd1="cat diag/power.left.dat >> diag/powertags.left.dat";
	string cmd2="cat diag/power.right.dat >> diag/powertags.right.dat";
	system(cmd1.c_str());
	system(cmd2.c_str());

}


//static double *accum;

void mother(double rho[], int N, int time){
	
	//Split rho into left of- and right of DL
	//Assuming equal amount of space on each side of DL
	double dx=(ux-lx)/(N-1);
	int size=int(-lx/dx);	//Size of added plasma in array inidces
	int M=N; while(M>size) M/=2;	//Finds biggest M=2^n<size
	int j=N-1-M;			//Muched used index
	
	if(M>2048){
		cout<<"fourier.cpp: mother(): Size, M too big"<<endl;
		cout<<"Forcing exit"<<endl;
		exit(1);
	}

	static double rps[2048];	//Arrays for accumulation power-
	static double lps[2048];	//  spectrum over time
	double left[M], right[M];
	for(int i=0; i<M; ++i){
		left[i]=rho[i];		//Temporary arrays
		right[i]=rho[i+j];
	}
	
	//Find array index of modes
	const int modessize=50;		//Max modes to be found
	static int lmodes[modessize];	//Storing of mode wavenumbers
	static int rmodes[modessize];	//Storing of mode wavenumbers
	static int found;		//Count no. of found modes
	static int w=1;			//Count no. write to datafile
	static bool initial=true;	//When first time run, do initial

	if(initial){ 
		for(int i=0; i<modessize; ++i){ 
			rmodes[i]=-1;
			lmodes[i]=-1;
		}
		for(int i=0; i<M; ++i){ 
			rps[i]=0.0;
			lps[i]=0.0;
		}
		initial=false;
	}
	
	//Diagnostics: taking copies of left,right for plottable datafile
	int k;
	double Lsignal[M], Rsignal[M];	//Signal
	double Lwindow[M], Rwindow[M];	//Signal with window
	//double Lpowers[M], Rpowers[M];	//Powerspectrum
	//Saving original signal
	for(k=0;k<M;++k){ Lsignal[k]=left[k]; Rsignal[k]=right[k];}

	//Windowing
	windowing(left,M);
	windowing(right,M);
	for(k=0;k<M;++k){ Lwindow[k]=left[k]; Rwindow[k]=right[k];}
	//Powerspectrum	
	powerspectrum(left,M);
	powerspectrum(right,M);
	//for(k=0;k<M;++k) Lpowers[i]=left[k], Rpowers=right[k];
	

	//Accumulating over time
	for(int i=0; i<M; ++i){ 
		rps[i] += right[i];
		lps[i] += left[i];
	
	}

	if(10*time/N_time == w){
#ifdef MOTHERDIAG
		//Diagnostics
		char snu[13]; sprintf(snu, "%d", time);
		string file="diag/mother.signal.";
		file=file+snu+".dat";
		ofstream fout(file.c_str());

		double l,r;
		for(k=0,l=lx,r=lx+j*dx ; k<M ; ++k,l+=dx,r+=dx){
			fout	<< l << " "	//1.Left positions
				<< r << " "	//2.Right positions
				<< Lsignal[k] << " "
				<< Rsignal[k] << " "
				<< Lwindow[k] << " "
				<< Rwindow[k] << " "
				<<endl;
		}
		fout.close();
#endif	
		//Search for new modes
		found+=findmodes(rps,M,rmodes,modessize, "right");	
		found+=findmodes(lps,M,lmodes,modessize, "left");	
		//Write power of modes to file
		writemodes(rps,M,rmodes,time*dt,"right",dx);
		writemodes(lps,M,lmodes,time*dt,"left",dx);
		//Reset accumulated arrays
		for(int i=0; i<M; ++i){
			rps[i]=0.0;
			lps[i]=0.0;
		}

		++w;
	}

	
}






void testmother(){
	int N=1024;
	double RHO[N];

	double a=1, b=1, c=.0, d=0.000;
	void mother(double *, int, int);
	double dx=8.*pi/N;
	for(int k=0; k<40; ++k){
		for(int i=0; i<N; ++i)
			RHO[i]=	(a+=d*drand48())*sin(i*dx)	
	//		RHO[i]=	sin(3.*i*dx)+cos(7.*i*dx)	;
				+(b+=d*drand48())*cos(10*i*dx)
				+(c+=d*drand48())*cos(34*i*dx)
				+0.5*(2.*drand48()-1.);
		mother(RHO, N, k);
		d=0.0002;
	}

	cout << "Fourier(mother) routines test finished, exiting."<<endl;
	exit(0);

}






void FFTphidl(double data[], int N, double delta){
	void windowing(double*, int);
	void powerspectrum(double*, int);

	//Finding highest M=2^p < N
	int M=int(pow(2.,20.));
	for(double p=20.; M>N; --p)	M=int(pow(2.,p));
	
	//Declaring data array of M
	double f[M];	int i,j;
	//Copying last M elements of data
	for(i=0, j=N-M; i<M; ++i, ++j) 	f[i]=data[j];
	
	//Powerspectrum with window
	windowing(f,M);
	powerspectrum(f,M);
	
	//Write results to file
	ofstream pout("diag/phidl.powerspectrum.dat");
	for(i=0;i<M;++i)	pout<<2.*pi*i/(delta*M)<<" "<<f[i]<<endl;
	pout.close();
}


