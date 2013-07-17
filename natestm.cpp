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

using namespace std;


template <class T>
void calcKWeights(T* weights, int rows){
	int cols = rows;
	int row, col;
	int rowMiddle = rows/2;
	int colMiddle = cols/2;
	int rowStart = 0;
	int colStart = 1;
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
void calcG11(T wMax, int numSpecVoltages, complex<T> gamma, T* quasi_k, T* delta_k, T* k_weights, complex<T>* G11, int rows){
	printf("STARTING CALCULATION\n");
	boost::timer t;
	int i;
	complex<T> t1;
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
		G11[i] = calcG011(w, gamma, quasi_k, delta_k, k_weights, rows);
	}
	
	double elapsed_time = t.elapsed();
	cout << "G11 elapsed time: " << elapsed_time << "\n";
	printf("ENDING CALCULATION\n");
}

template <class T>
void calcG11_linear(T wMax, int numSpecVoltages, complex<T> gamma, T* quasi_k, T* delta_k, T* k_weights, complex<T>* G11, int rows){
	printf("STARTING CALCULATION\n");
	boost::timer t;
	int i;
	complex<T> t1;
	T w, stepSize;
	if(numSpecVoltages > 1)
		stepSize = 2*wMax/(numSpecVoltages-1);
	else 
		stepSize = 0;
	//printf("vMax: %.16lf, vstep: %.16lf, maxV: %.16lf\n", wMax, stepSize, -wMax + (numSpecVoltages-1)*stepSize);
	for(i=0;i<numSpecVoltages;i++){
		w = -wMax + i*stepSize;
		//cout << "voltage is " << w << endl;
		G11[i] = calcG011(w, gamma, quasi_k, delta_k, k_weights, rows);
	}
	double elapsed_time = t.elapsed();
	cout << "G11 elapsed time: " << elapsed_time << "\n";
	printf("ENDING CALCULATION\n");
}

template <class T>
T fermi(T energy, T fermiEnergy, T beta){
	return 1/(pow(M_E, (energy-fermiEnergy)*beta) + 1.0);	
	
}




template <class T>
void arrayRange(T* array, T min, T max, int numSteps){
	T stepSize = (max-min)/(numSteps-1);
	//printf("step: %lf max: %lf\n", stepSize, min + stepSize*(numSteps-1));
	for(int i=0; i<numSteps; i++){
		array[i] = min + i*stepSize;
	}
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

	init_glui(argc, argv);

	return 0;
}

void idle(){
	usleep(1000);
}

void allocateMemory(ScanUserData* scanUserData){
	printf("allocating scan memory\n");

	// if(scanUserData->filename){
	// 	delete[] scanUserData->filename;
	// 	scanUserData->filename = 0;
		
	// }

	if(scanUserData->gaps){
		delete[] scanUserData->gaps;
		scanUserData->gaps = 0;
		
	}

	if(scanUserData->bandEnergy){
		delete[] scanUserData->bandEnergy;
		scanUserData->bandEnergy = 0;
		
	}

	if(scanUserData->quasiEnergy){
		delete[] scanUserData->quasiEnergy;
		scanUserData->quasiEnergy = 0;
		
	}

	if(scanUserData->kx){
		delete[] scanUserData->kx;
		scanUserData->kx = 0;
		
	}

	if(scanUserData->ky){
		delete[] scanUserData->ky;
		scanUserData->ky = 0;
		
	}

	if(scanUserData->k_weights){
		delete[] scanUserData->k_weights;
		scanUserData->k_weights = 0;
		
	}

	if(scanUserData->G11){
		delete[] scanUserData->G11;
		scanUserData->G11 = 0;
		
	}

	if(scanUserData->spec){
		delete[] scanUserData->spec;
		scanUserData->spec = 0;
		
	}

	// scanUserData->filename = new char[100];

	scanUserData->gaps = new double[scanUserData->nx*scanUserData->nx];
	scanUserData->bandEnergy = new double[scanUserData->nx*scanUserData->nx];
	scanUserData->quasiEnergy = new double[scanUserData->nx*scanUserData->nx];
	scanUserData->kx = new double[scanUserData->nx];
	scanUserData->ky = new double[scanUserData->nx];

	scanUserData->k_weights = new double[scanUserData->nx*scanUserData->nx];

	scanUserData->G11 = new complex<double>[scanUserData->numSpecVoltages];
	scanUserData->spec = new double[scanUserData->numSpecVoltages];


	initArray(scanUserData->gaps, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->bandEnergy, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->quasiEnergy, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->k_weights, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->kx, 0, scanUserData->nx);
	initArray(scanUserData->ky, 0, scanUserData->nx);
	initArray(scanUserData->G11, 0, scanUserData->numSpecVoltages);
	initArray(scanUserData->spec, 0, scanUserData->numSpecVoltages);

}


void onCalculateG011(int id){
	allocateMemory(&scanUserData);

	arrayRange(scanUserData.kx, -PI, PI - (2*PI)/scanUserData.nx, scanUserData.nx);
	arrayRange(scanUserData.ky, -PI, PI - (2*PI)/scanUserData.nx, scanUserData.nx);
	
	dWaveGap(scanUserData.gap0, scanUserData.kx, scanUserData.ky, scanUserData.gaps, scanUserData.nx);
	scanUserData.u  = -4.0*(scanUserData.t2 - scanUserData.t3) - scanUserData.Em;
	calcBandEnergy(	scanUserData.kx, scanUserData.ky, 
					scanUserData.t1, scanUserData.t2, scanUserData.t3, 
					scanUserData.u, scanUserData.bandEnergy, scanUserData.nx);

	scanUserData.bandMin = getMin(scanUserData.bandEnergy, scanUserData.nx*scanUserData.nx);
	scanUserData.bandMax = getMax(scanUserData.bandEnergy, scanUserData.nx*scanUserData.nx);

	calcQuasiEnergy(scanUserData.bandEnergy, scanUserData.gaps, scanUserData.quasiEnergy, scanUserData.nx);
	scanUserData.quasiMin = getMin(scanUserData.quasiEnergy, scanUserData.nx*scanUserData.nx);
	scanUserData.quasiMax = getMax(scanUserData.quasiEnergy, scanUserData.nx*scanUserData.nx);

	complex<double> gamma(0.0,scanUserData.gamma);
	double w = 0.0;

	calcKWeights(scanUserData.k_weights, scanUserData.nx);
	/*actually uses bare dispersion, not quasi dispersion*/
	calcG11((double)scanUserData.vMax, scanUserData.numSpecVoltages, gamma, scanUserData.bandEnergy, scanUserData.gaps, scanUserData.k_weights, scanUserData.G11, scanUserData.nx);
	//calcG11_linear((double)scanUserData.vMax, scanUserData.numSpecVoltages, gamma, scanUserData.bandEnergy, scanUserData.gaps, scanUserData.k_weights, scanUserData.G11, scanUserData.nx);
	
	
	cout << "finished G11\n";
	double v;
	double stepSize = 2*scanUserData.vMax/(scanUserData.numSpecVoltages-1);

	
	double n_squared = scanUserData.nx*scanUserData.nx;
	cout << "calculating spec values\n";
	for(int i=0; i<scanUserData.numSpecVoltages; i++){
		v = -scanUserData.vMax + i*stepSize;
		scanUserData.spec[i] = -(1/(n_squared*PI))*scanUserData.G11[i].imag();
		//printf("%lf %lf\n", v, scanUserData.spec[i]);

	}	
	cout << "finished spec values\n";

	cout << "calculating max min\n";

	scanUserData.specMin = getMin(scanUserData.spec, scanUserData.numSpecVoltages);
	scanUserData.specMax = getMax(scanUserData.spec, scanUserData.numSpecVoltages);
	cout << "finished max min\n";
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

