//gcc natestm.cpp -o test -lGL -lGLU -lglut -lglui
//gcc -g natestm.cpp -o test -lGL -lGLU -lglut -lglui

#include <stdio.h>
#include <string>
#include "natestm.h"
#include "graphics.h"
#include <GL/glui.h>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include "define_constants.h"
//#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include "lmcurve.h"

/* removing cuda functionality for the moment */
//#include <cuda.h>
//#include <cuComplex.h>
//#include "g_kernel.h"

using namespace std;

void fit(void);
template <class T>
int calcKWeights(T* weights, int rows){
	int cols = rows;
	int row, col;
	int rowMiddle = rows/2;
	int colMiddle = cols/2;
	int rowStart = 0;
	int colStart = 1;
	int length=0;
	weights[0] = 1;
	//for(row = 0; row<rowMiddle; row++){ 
	//The first row is special and we calculate it independently. 
	row = 0;
	for(col = 1; col<colMiddle; col++){
		weights[col] = 4;
		
	}
	weights[colMiddle] = 2;

	for(row = 1; row<rowMiddle; row++){
		weights[row*cols + row] = 4;
		for(col = row+1; col<colMiddle; col++){
			weights[row*cols + col] = 8;
		}
		weights[row*cols + colMiddle] = 4;
	}
	weights[rowMiddle*cols + colMiddle] = 1;

	for(row = 0; row<rows; row++){	
		for(col = 0; col<cols; col++){
			length++;
		}
	}
	return length;
}

template <class T, class U>
void dWaveGap(U gap0, T* kx, T* ky, T* gap, int rows){
	int cols = rows;
	int row, col;
	for(row = 0; row<rows; row++){
		for(col = 0; col<cols; col++){
			gap[row*cols+col]= gap0*(cos(kx[col])-cos(ky[row]))/2.0;
		}
	}
}

template <class T, class U>
void calcBandEnergy(T* kx, T* ky, U t1, U t2, U t3, U u, T* band, int rows){
	int cols = rows;
	int row, col;
	for(row = 0; row<rows; row++){
		for(col = 0; col<cols; col++){
			band[row*cols+col]= 2.0*t1*(cos(kx[col])+cos(ky[row]))+
								4.0*t2*(cos(kx[col])*cos(ky[row]))+
								2.0*t3*(cos(2*kx[col])+cos(2*ky[row]))
								- u;
		}
	}
}

template <class T>
void calcQuasiEnergy(T* band, T* gap, T* quasi, int rows){
	int cols = rows;
	int row, col;
	for(row = 0; row<rows; row++){
		for(col = 0; col<cols; col++){
			quasi[row*cols+col]= sqrt(	band[row*cols+col]*band[row*cols+col]+
										gap[row*cols+col]*gap[row*cols+col]);
		}
	}
}

/*this functions calculates G11 for one energy, the next function for a range of energies*/
template <class T>
complex<T> calcG011(T w, complex<T> gamma, T* quasi_k, T* delta_k, T* k_weights, int rows){
	int cols = rows;
	int row, col, index;
	complex<T> t1, g011(0.0,0.0), wg;
	T weight;
	wg = w + gamma;
	for(row = 0; row<rows/2+1; row++){
		for(col = 0; col<cols/2+1; col++){
			index = row*cols + col;
			weight = k_weights[index];
			if(weight == 0){
				continue;
			}
			t1 = pow(delta_k[index], 2)/(wg + quasi_k[index]);
			
			//g011 += 1.0/(w + gamma - quasi_k[row*cols+col] - t1);
			g011 += weight/(wg - quasi_k[index] - t1);
		}
	}
	return g011;
}

template <class T>
complex<T> calcG011_with_lorentzian(T w, complex<T> gamma, T* quasi_k, T* delta_k, T* k_weights, T lorentz_amplitude, T lorentz_energy, complex<T> lorentz_gamma, int rows){
	int cols = rows;
	int row, col, index;

	complex<T> t1, g011(0.0,0.0), wg, sigma_w;
	T weight;
	wg = w + gamma;
	for(row = 0; row<rows/2+1; row++){
		for(col = 0; col<cols/2+1; col++){
			index = row*cols + col;
			weight = k_weights[index];
			if(weight == 0){
				continue;
			}
			sigma_w = lorentzian(lorentz_amplitude, lorentz_energy, w, lorentz_gamma);
			
			t1 = pow(delta_k[index], 2)/(wg + quasi_k[index] -sigma_w);
			
			
			//g011 += 1.0/(w + gamma - quasi_k[row*cols+col] - t1);
			g011 += weight/(wg - quasi_k[index] - sigma_w - t1);
			
		}
	}
	return g011;
}

template <class T>
void calcG11(T wMax, int numSpecVoltages, complex<T> gamma, T* quasi_k, T* delta_k, T* k_weights,  T lorentz_amplitude, T lorentz_energy, complex<T> lorentz_gamma, complex<T>* G11, int rows){
	printf("STARTING CALCULATION\n");
	boost::timer t;
	int i;
	complex<T> t1;
	// complex<T> lorentz_gamma;
	// T lorentz_energy=1;
	// T lorentz_amplitude=10000000;
	T w, stepSize;
	if(numSpecVoltages > 1)
		stepSize = 2*wMax/(numSpecVoltages-1);
	else 
		stepSize = 0;

	//printf("vMax: %.16lf, vstep: %.16lf, maxV: %.16lf\n", wMax, stepSize, -wMax + (numSpecVoltages-1)*stepSize);
	#pragma omp parallel for
	for(i=0;i<numSpecVoltages;i++){
		w = -wMax + i*stepSize;
		//cout << "voltage is " << w << endl;
		G11[i] = calcG011_with_lorentzian(w, gamma, quasi_k, delta_k, k_weights, lorentz_amplitude, lorentz_energy, lorentz_gamma, rows);
	}
	
	double elapsed_time = t.elapsed();
	cout << "G11 elapsed time: " << elapsed_time << "\n";
	printf("ENDING CALCULATION\n");
}

template <class T>
void calcDOS(complex<T>* G11, T* DOS, int numSpecVoltages){
	int n_squared = numSpecVoltages*numSpecVoltages;
	for(int i=0; i<numSpecVoltages; i++){
		DOS[i] = -(1/(n_squared*PI))*G11[i].imag();
	}	
}

template <class T>
T fermi(T energy, T fermiEnergy, T beta){
	return 1/(pow(M_E, (energy-fermiEnergy)*beta) + 1.0);		
}

template <class T>
complex<T> lorentzian(T Amplitude, T energy, T w, complex<T> gamma){
	complex<T> first_term, second_term;
	first_term = Amplitude/(w - energy + gamma);
	second_term = Amplitude/(-w - energy + gamma);
	return first_term + second_term;
}

template <class T>
void arrayRange(T* array, T min, T max, int numSteps){
	T stepSize = (max-min)/(numSteps-1);
	//printf("step: %lf max: %lf\n", stepSize, min + stepSize*(numSteps-1));
	for(int i=0; i<numSteps; i++){
		array[i] = min + i*stepSize;
	}
}


template <class T>
T* free_and_reallocate(T* array, int length){
	if(array){
		delete[] array;
	}

	array = new T[length];
	initArray(array, 0, length);
	return array;
}



ScanUserData scanUserData;

ScanUserData* getScanUserData(void){
	return &scanUserData;
}



int main(int argc, char* argv[])
{
	
	scanUserData.gaps = scanUserData.bandEnergy = scanUserData.quasiEnergy = 0;


	scanUserData.quasiMin = scanUserData.bandMin = -1.0;
	scanUserData.quasiMax = scanUserData.bandMax = 1.0;

	scanUserData.t1 = -882;
	scanUserData.t2 = 239;
	scanUserData.t3 = -14;
	scanUserData.Em = -26;
	scanUserData.u  = -4.0*(scanUserData.t2 - scanUserData.t3) - scanUserData.Em;
	scanUserData.vMax  = 100.0;
	scanUserData.gamma = 1.0;

	scanUserData.numSpecVoltages = 256;
	scanUserData.nx = 256;
	scanUserData.gap0 = 50;

	scanUserData.lorentz_amplitude = 100;
	scanUserData.lorentz_energy = 100;
	scanUserData.lorentz_gamma = 1.0;

	init_glui(argc, argv);

	return 0;
}

void idle(){
	usleep(1000);
}

void allocateMemory(ScanUserData* scanUserData){
	printf("allocating scan memory\n");
	int nx_squared = scanUserData->nx*scanUserData->nx;
	int nx = scanUserData->nx;

	scanUserData->gaps = free_and_reallocate(scanUserData->gaps, nx_squared);
	scanUserData->bandEnergy = free_and_reallocate(scanUserData->bandEnergy, nx_squared);
	scanUserData->quasiEnergy = free_and_reallocate(scanUserData->quasiEnergy, nx_squared);
	scanUserData->kx = free_and_reallocate(scanUserData->kx, nx);
	scanUserData->ky = free_and_reallocate(scanUserData->ky, nx);
	scanUserData->k_weights = free_and_reallocate(scanUserData->k_weights, nx_squared);
	scanUserData->G11 = free_and_reallocate(scanUserData->G11, scanUserData->numSpecVoltages);
	scanUserData->spec = free_and_reallocate(scanUserData->spec, scanUserData->numSpecVoltages);

}

void calc_dos_for_fit( const double *par, int m_dat,
        const void *data, double *fvec, int *userbreak ){
	double gap0, t1, t2, t3, u, gam;
	double lorentz_amplitude, lorentz_energy, l_gamma;

	int i = 0;
	gap0 	= par[i++];
	t1   	= par[i++];
	t2   	= par[i++];
	t3   	= par[i++];
	u   	= par[i++];
	gam   		= par[i++];
	complex<double> gamma(0.0,gam);
	lorentz_amplitude = par[i++];
	lorentz_energy   	= par[i++];
	l_gamma   	= par[i++];
	complex<double> lorentz_gamma(0.0,l_gamma);

	u  = -4.0*(t2 - t3) - scanUserData.Em;

	dWaveGap(gap0, scanUserData.kx, scanUserData.ky, scanUserData.gaps, scanUserData.nx);
	calcBandEnergy(	scanUserData.kx, scanUserData.ky, 
									t1, t2, t3, u, scanUserData.bandEnergy, m_dat);

	calcQuasiEnergy(scanUserData.bandEnergy, scanUserData.gaps, scanUserData.quasiEnergy, scanUserData.nx);

	
	calcG11((double)scanUserData.vMax, scanUserData.numSpecVoltages, gamma, 
					scanUserData.bandEnergy, scanUserData.gaps, scanUserData.k_weights, 
					lorentz_amplitude, lorentz_energy, lorentz_gamma, 
					scanUserData.G11,  scanUserData.nx);
	cout << "finished G11\n";
	calcDOS(scanUserData.G11, scanUserData.spec, scanUserData.numSpecVoltages);
}

void evaluate_dos( const double *par, int m_dat,
        const void *data, double *fvec, int *userbreak ){
        /* for readability, explicit type conversion */
	ScanUserData *D;
	D = (ScanUserData*)data;

	calc_dos_for_fit( par, m_dat, data, fvec, userbreak );

	// int i;
	// // for ( i = 0; i < m_dat; i++ ){
	// // 	fvec[i] = D->y[i] - D->f( D->tx[i], D->tz[i], par );
	// // }
}

void onCalculateG011(int id){

	/* ------------------------- SETUP -------------------------------- */
	allocateMemory(&scanUserData);
	calcKWeights(scanUserData.k_weights, scanUserData.nx);
	arrayRange(scanUserData.kx, -PI, PI - (2*PI)/scanUserData.nx, scanUserData.nx);
	arrayRange(scanUserData.ky, -PI, PI - (2*PI)/scanUserData.nx, scanUserData.nx);

	/* ------------------------- CALCULATION --------------------------- */

	// dWaveGap(scanUserData.gap0, scanUserData.kx, scanUserData.ky, scanUserData.gaps, scanUserData.nx);
	// scanUserData.u  = -4.0*(scanUserData.t2 - scanUserData.t3) - scanUserData.Em;
	// calcBandEnergy(	scanUserData.kx, scanUserData.ky, 
	// 								scanUserData.t1, scanUserData.t2, scanUserData.t3, 
	// 								scanUserData.u, scanUserData.bandEnergy, scanUserData.nx);

	// calcQuasiEnergy(scanUserData.bandEnergy, scanUserData.gaps, scanUserData.quasiEnergy, scanUserData.nx);

	// complex<double> gamma(0.0,scanUserData.gamma);
	// double w = 0.0;
	// complex<double> lorentz_gamma(0.0,scanUserData.lorentz_gamma);
	// calcG11((double)scanUserData.vMax, scanUserData.numSpecVoltages, gamma, scanUserData.bandEnergy, scanUserData.gaps, scanUserData.k_weights, (double)scanUserData.lorentz_amplitude, (double)scanUserData.lorentz_energy, lorentz_gamma, scanUserData.G11,  scanUserData.nx);
	// cout << "finished G11\n";
	// calcDOS(scanUserData.G11, scanUserData.spec, scanUserData.numSpecVoltages);
	int m_dat = scanUserData.numSpecVoltages;
	double* fvec;
	int* userbreak;

	

	int i = 0;
	double par[9];
	par[i++] = scanUserData.gap0;
	par[i++] = scanUserData.t1;
	par[i++] = scanUserData.t2;
	par[i++] = scanUserData.t3;
	par[i++] = scanUserData.u; 
	par[i++] = scanUserData.gamma;
	par[i++] = scanUserData.lorentz_amplitude;
	par[i++] = scanUserData.lorentz_energy;
	par[i++] = scanUserData.lorentz_gamma;

	calc_dos_for_fit( par, m_dat,
        (void*)&scanUserData, fvec, userbreak );


	/* -------------------------- POST PROCESSING --------------------------------------------- */
	cout << "finished spec values\n";
	scanUserData.specMin = getMin(scanUserData.spec, scanUserData.numSpecVoltages);
	scanUserData.specMax = getMax(scanUserData.spec, scanUserData.numSpecVoltages);
	cout << "finished max min\n";
	cout << "starting fitting\n";
	fit();
	cout << "finished fitting\n";
	post_redisplay();
}
	
void onSave(int id){
	printf("onsave\n");
	printf("%s\n", scanUserData.filename);
	double v;
	double stepSize = 2*scanUserData.vMax/(scanUserData.numSpecVoltages-1);
  	ofstream file;
  	file.open(scanUserData.filename);
	
	double n_squared = scanUserData.nx*scanUserData.nx;
	for(int i=0; i<scanUserData.numSpecVoltages; i++){
		v = -scanUserData.vMax + i*stepSize;
		file << v << "    " << scanUserData.spec[i] << "\n";
	}	
	file.close();	
	return;
}


void onCalcWeights(int id){
	calcKWeights(scanUserData.k_weights, scanUserData.nx);
	printArray(scanUserData.k_weights, scanUserData.nx, scanUserData.nx);
}


double f( double t, const double *p )
{
		getScanUserData();
    return p[0] + p[1]*t + p[2]*t*t;
}

void fit(void){
		
    int n = 3; /* number of parameters in model function f */
    double par[3] = { 100, 0, -10 }; /* really bad starting value */
    
    /* data points: a slightly distorted standard parabola */
    int m = 9;
    int i;
    double t[9] = { -4., -3., -2., -1.,  0., 1.,  2.,  3.,  4. };
    double y[9] = { 16.6, 9.9, 4.4, 1.1, 0., 1.1, 4.2, 9.3, 16.4 };

    lm_control_struct control = lm_control_double;
    lm_status_struct status;
    control.verbosity = 9;

    printf( "Fitting ...\n" );
    /* now the call to lmfit */
    lmcurve( n, par, m, t, y, f, &control, &status );
        
    printf( "Results:\n" );
    printf( "status after %d function evaluations:\n  %s\n",
            status.nfev, lm_infmsg[status.outcome] );

    printf("obtained parameters:\n");
    for ( i = 0; i < n; ++i)
        printf("  par[%i] = %12g\n", i, par[i]);
    printf("obtained norm:\n  %12g\n", status.fnorm );
    
    printf("fitting data as follows:\n");
    for ( i = 0; i < m; ++i)
        printf( "  t[%2d]=%4g y=%6g fit=%10g residue=%12g\n",
                i, t[i], y[i], f(t[i],par), y[i] - f(t[i],par) );

    return;
  }
