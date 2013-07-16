#ifndef NATESTM_H
#define NATESTM_H

/*disable an annoying truncation warning*/
#pragma warning(disable: 4786)
#pragma warning(disable: 4275)
#pragma warning(disable: 4251)
#include <GL/glui.h>
#include <complex>




#define V_MIN -1000.0
#define V_MAX 1000.0

#define NUM_VOLTAGES_MIN 0
#define NUM_VOLTAGES_MAX 1000
#define NUM_XY_MIN 0
#define NUM_XY_MAX 10000

#define WAIT_SPINNER_SPEED 0.1
#define SPINNER_SPEED 0.1

#define SMOOTH_SPEC_MIN 0
#define SMOOTH_SPEC_MAX 1000

#define GAMMA_MAX 1000.0
#define GAMMA_MIN 0.0001




#define MOVE_ACTION 0
#define TOPO_ACTION 1
#define TOPO_AND_SPEC_ACTION 3
#define SPEC_ACTION 2

#define HEADER_SIZE 1000

#define T_MIN -10000.0f
#define T_MAX 10000.0f
#define U_MIN -10000.0f
#define U_MAX 10000.0f

/*dap time is in microseconds, this is the multiplier we will use to convert between seconds*/

#define PI 3.14159265

struct imageScale{
	double min;
	double max;
	double mean;
	double dev;
};

std::complex<double> a;

struct ScanUserData{

	int numSpecVoltages, numSpecSamples, numSpec, nx, nv;

	double vMin, vMax;
	int	action;
	
	double* kx;
	double* ky;
	double* gaps;
	double* bandEnergy;
	double* quasiEnergy;
	double* spec;
	std::complex<double>* G11;


	double bandMin;
	double bandMax;
	double quasiMin;
	double quasiMax;
	double specMin;
	double specMax;
	float gamma;

	double GMin;
	double GMax;

	float t1, t2, t3, u, Em;

	float gap0;

	int smoothSpecWidth;


	int specPixelSize;
	int topoSize;
	double specMapSize;


	char* dapFileName;


	int boolAutoSave;

};

/*graphics data holds values for the graphics, updateAfterAverages tells us how many times
  to run the simulation before displaying on the screen, mainWindow is an ID for the
  graphics library to know which window to draw in*/
struct graphicsData{
   int mainWindow;
   GLUI *glui;

	GLUI_Spinner *vMinSpinner;
	GLUI_Spinner *vMaxSpinner;
	GLUI_Spinner *numSpecVoltagesSpinner;

	GLUI_Spinner *nxSpinner;
	GLUI_Spinner *nySpinner;
	GLUI_Spinner *gapSpinner;
	GLUI_Spinner *t1Spinner;
	GLUI_Spinner *t2Spinner;
	GLUI_Spinner *t3Spinner;
	GLUI_Spinner *uSpinner;
	GLUI_Spinner *EmSpinner;
	GLUI_Spinner *vSpinner;
	GLUI_Spinner *gammaSpinner;
	GLUI_Spinner *nvSpinner;

	GLUI_Spinner *smoothSpecWidthSpinner;


	GLUI_Button *moveButton;

	GLUI_Rollout* topoRollout;
	GLUI_Rollout* specRollout;

	



	int width;
	int height;
	int autoScaleTopo;
	int autoScaleSpec;

	int boolAutoScaleSpec;
	int boolAutoScaleDeriv;
	int boolAutoScaleTopo;
	

};

template <class T>
T getMax(T* array, int length){
	if(length == 0)
		return 0;
	T max = array[0];
	for(int i=0; i<length; i++){
		if(array[i]>max)
			max = array[i];
	}
	return max;
}

template <class T>
T getMin(T* array, int length){
	if(length == 0)
		return 0;
	T min = array[0];
	for(int i=0; i<length; i++){
		if(array[i]<min)
			min = array[i];
	}
	return min;
}

template <class T>
double getMean(T* array, int length){
	if(length == 0)
		return 0;
	double mean = 0;
	for(int i=0; i<length; i++){
		mean += array[i];
	}
	return 1.0*mean/length;
}

template <class T>
double getVariance(T* array, int length){
	if(length == 0)
		return 0;
	double mean = getMean(array, length);
	double var=0;
	for(int i=0; i<length; i++){
		var += (array[i]-mean)*(array[i]-mean);
	}
	return var/length;


}

template <class T>
void deriv3(T* spec, double* derivSpec, int length){
	for(int i=1; i<(length-1); i++){
		derivSpec[i] = (spec[i+1]-spec[i-1])/2.0;
	}
	derivSpec[0] = spec[1]-spec[0];
	derivSpec[length-1] = spec[length-1]-spec[length-2];
}

template <class T>
void linsmooth(T* spec, double* smoothSpec, int halfWidth, int length){
	int i,j,k;
	double avg;

	avg = 0;
	for(j=0; j<=halfWidth; j++){
		avg += spec[j];
		
	}
	smoothSpec[0] = avg/j;

	for(i=1; i<=halfWidth; i++){
		avg += spec[j+i-1];
		smoothSpec[i] = avg/(j+i);
	}
	for(i=halfWidth+1; i <(length-halfWidth); i++){
	   avg += spec[i+halfWidth];
	   avg -= spec[i-halfWidth-1];
	   smoothSpec[i] = avg/(2*halfWidth+1);
	}
	k=1;
	for(i=(length-halfWidth); i<length; i++){
	   avg -= spec[i-halfWidth-1];
      smoothSpec[i] = avg/(2*halfWidth+1-k);
      k++;
	}     
}

template <class T>
void initArray(T* array, int initVal, int length){
	if(length == 0)
		return;
	if(array == 0)
		return;
	for(int i=0; i<length; i++){
		array[i]=0;
	}
}

template <class T>
void copyArray(T* array, T* arrayCopy, int length){
	if(length == 0)
		return;
	for(int i=0; i<length; i++){
		arrayCopy[i]=array[i];
	}
}


#endif