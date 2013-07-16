//gcc natestm.cpp -o test -lGL -lGLU -lglut -lglui
//gcc -g natestm.cpp -o test -lGL -lGLU -lglut -lglui

#include <stdio.h>
#include <string>
#include "natestm.h"
#include <GL/glui.h>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include "define_constants.h"

using namespace std;

void draw(void);
void drawTopo(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);
void drawSpectro(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);
void drawEcut(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);

void gluiCallback(int id);

void idle(void);
void allocateMemory(ScanUserData* scanUserData);
void onResize(int width, int height);
void onCalculateGap(int id);
void onCalculateBand(int id);
void onCalculateQuasi(int id);
void onCalculateG011(int id);
void onSave(int id);

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
complex<T> calcG011(T w, complex<T> gamma, T* quasi_k, T* delta_k, int rows){
	int cols = rows;
	int row, col, index;
	complex<T> t1, g011(0.0,0.0), wg;
	wg = w + gamma;
	for(row = 0; row<rows; row++){
		for(col = 0; col<cols; col++){
			index = row*cols + col;
			t1 = pow(delta_k[index], 2)/(wg + quasi_k[index]);
			
			//g011 += 1.0/(w + gamma - quasi_k[row*cols+col] - t1);
			g011 += 1.0/(wg - quasi_k[index] - t1);
		}
	}
	return g011;
}

template <class T>
void calcG11(T wMax, int numSpecVoltages, complex<T> gamma, T* quasi_k, T* delta_k, complex<T>* G11, int rows){
	
	int i;
	complex<T> t1;
	T w, stepSize;
	if(numSpecVoltages > 1)
		stepSize = 2*wMax/(numSpecVoltages-1);
	else 
		stepSize = 0;
	printf("vMax: %.16lf, vstep: %.16lf, maxV: %.16lf", wMax, stepSize, -wMax + (numSpecVoltages-1)*stepSize);
	for(i=0;i<numSpecVoltages;i++){
		w = -wMax + i*stepSize;
		//cout << "voltage is " << w << endl;
		G11[i] = calcG011(w, gamma, quasi_k, delta_k, rows);
	}
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



graphicsData myGraphicsData;


void initGL(void){
	glShadeModel(GL_SMOOTH);
	float glcolors[4] = {-1,-1,-1,-1};
	glClearColor(0,0.5, 0.6, 0);
   
   glGetFloatv(GL_COLOR_CLEAR_VALUE, glcolors);
   //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glClearDepth(2.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	
	
	
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
	cout << "U = " << scanUserData.u << "\n";
	scanUserData.vMax  = 100.0;
	scanUserData.gamma = 1.0;

	myGraphicsData.width = myGraphicsData.height = 700;
	glutInit(&argc, argv);
	glutInitWindowSize(myGraphicsData.width,myGraphicsData.height);
	glutInitDisplayMode(GLUT_DOUBLE);
	myGraphicsData.mainWindow = glutCreateWindow("natemodel");
	glutDisplayFunc(draw);
	glutReshapeFunc(onResize);
	
	initGL();
	
	

	scanUserData.numSpecVoltages = 256;

	scanUserData.nx = 256;
	scanUserData.gap0 = 50;

	scanUserData.action = TOPO_ACTION;

	scanUserData.smoothSpecWidth = 3;

	/****************************** USER INTERFACE *********************************************/

	myGraphicsData.glui = GLUI_Master.create_glui( "GLUI", 0,myGraphicsData.width+20,0 );	

	myGraphicsData.topoRollout = myGraphicsData.glui->add_rollout("Model Parameters", true);
	/****************************** TOPOGRAPHY *********************************************/

	myGraphicsData.nxSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "nx:", GLUI_SPINNER_INT, &(scanUserData.nx) );
	myGraphicsData.nxSpinner->set_int_limits( NUM_XY_MIN, NUM_XY_MAX );

	myGraphicsData.gapSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "gap0:", GLUI_SPINNER_FLOAT, &(scanUserData.gap0) );
	myGraphicsData.gapSpinner->set_float_limits( 0.0, 100.0 );

	myGraphicsData.t1Spinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "t1:", GLUI_SPINNER_FLOAT, &(scanUserData.t1) );
	myGraphicsData.t1Spinner->set_float_limits( T_MIN, T_MAX );

	myGraphicsData.t2Spinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "t2:", GLUI_SPINNER_FLOAT, &(scanUserData.t2) );
	myGraphicsData.t2Spinner->set_float_limits( T_MIN, T_MAX );

	myGraphicsData.t3Spinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "t3:", GLUI_SPINNER_FLOAT, &(scanUserData.t3) );
	myGraphicsData.t3Spinner->set_float_limits( T_MIN, T_MAX );

	myGraphicsData.EmSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "Em:", GLUI_SPINNER_FLOAT, &(scanUserData.Em) );
	myGraphicsData.EmSpinner->set_float_limits( T_MIN, T_MAX );

	myGraphicsData.uSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "u:", GLUI_SPINNER_FLOAT, &(scanUserData.u) );
	myGraphicsData.uSpinner->set_float_limits( U_MIN, U_MAX );

	myGraphicsData.vSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "vMax:", GLUI_SPINNER_FLOAT, &(scanUserData.vMax) );
	myGraphicsData.vSpinner->set_float_limits( 0, V_MAX );

	myGraphicsData.gammaSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "gamma:", GLUI_SPINNER_FLOAT, &(scanUserData.gamma) );
	myGraphicsData.gammaSpinner->set_float_limits( GAMMA_MIN, GAMMA_MAX );

	myGraphicsData.nvSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.topoRollout,  "num Vs:", GLUI_SPINNER_INT, &(scanUserData.numSpecVoltages) );
	myGraphicsData.nvSpinner->set_int_limits( NUM_XY_MIN, NUM_XY_MAX );


	myGraphicsData.glui->add_button_to_panel(myGraphicsData.topoRollout, "calculate G011", -1, onCalculateG011);
	myGraphicsData.glui->add_edittext_to_panel(myGraphicsData.topoRollout, "Filename", GLUI_EDITTEXT_TEXT, &scanUserData.filename);
	myGraphicsData.glui->add_button_to_panel(myGraphicsData.topoRollout, "Save", -1, onSave);

	/****************************** SPECTROSCOPY *********************************************/
	
	// myGraphicsData.vMinSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.specRollout, "vMin:", GLUI_SPINNER_FLOAT, &(scanUserData.vMin) );
	// myGraphicsData.vMinSpinner->set_int_limits( V_MIN, V_MAX );
	
	
	// myGraphicsData.glui->add_statictext_to_panel(myGraphicsData.specRollout, "         vMax == -vMin");

	// myGraphicsData.numSpecVoltagesSpinner = myGraphicsData.glui->add_spinner_to_panel(myGraphicsData.specRollout,"nv:", GLUI_SPINNER_INT, &(scanUserData.numSpecVoltages) );
	// myGraphicsData.numSpecVoltagesSpinner->set_int_limits( NUM_VOLTAGES_MIN, NUM_VOLTAGES_MAX );

	// myGraphicsData.glui->add_column(true);

	myGraphicsData.glui->sync_live();
	

	myGraphicsData.autoScaleTopo = 1;
	myGraphicsData.autoScaleSpec = 1;

	GLUI_Master.set_glutIdleFunc(idle);

	glutPostRedisplay();
	glutMainLoop();

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

	scanUserData->G11 = new complex<double>[scanUserData->numSpecVoltages];
	scanUserData->spec = new double[scanUserData->numSpecVoltages];


	initArray(scanUserData->gaps, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->bandEnergy, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->quasiEnergy, 0, scanUserData->nx*scanUserData->nx);
	initArray(scanUserData->kx, 0, scanUserData->nx);
	initArray(scanUserData->ky, 0, scanUserData->nx);
	initArray(scanUserData->G11, 0, scanUserData->numSpecVoltages);
	initArray(scanUserData->spec, 0, scanUserData->numSpecVoltages);

}

template <class T>
void displaySpecData(T *specData, int numPoints, double yMin, double yMax, double red, double green, double blue)
{
	
	int i;
	
	
	if(specData==0)
		return;
	
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//set up an orthographic projection 
	glOrtho(0, numPoints-1, yMin, yMax, -1, 1);
	
	//I draw the spectro graph in the forward and reverse directions
	//because there is a problem in how this OpenGL implementation fills
	//lines and a solution is to draw everything twice: once forward
	//and once reverse.  Still beats using windows code.
	glColor3f(red, green, blue);
	glBegin(GL_LINE_STRIP);
	//if(xStart <= xEnd){  //draw the spec without flipping the x axis
	
	for(i=0; i<numPoints; i++){
		glVertex3f((float)i, (float)specData[i], 0.0f); 
		
	}
	glEnd();
	
	glBegin(GL_LINE_STRIP);
	
	
	glEnd();
}

template <class T>
void drawImage(T* image, int cols, int rows, double min, double max){
	int i,j;
	float scale;
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//set up an orthographic projection 
	glOrtho(0, cols, 0, rows, -1, 1);

	//lode the model matrix where we can rotate or translate our image
	glMatrixMode(GL_MODELVIEW);
	//clear the model matrix stack
	glLoadIdentity(); 
	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_ALWAYS);
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	/*we don't have valid scan data*/
	if(image == 0)
		return;
	
	
	if( (max-min) != 0.0)
		scale = 1.0/(max-min);
	else
		scale = 1.0;
	
	glBegin(GL_QUADS);
	
	float color = 0;
	
	
	for(i=0; i<rows; i++){
		for(j=0; j<cols; j++){
			color = scale*(image[i*cols+j]-min);
			glColor3f(color, color, color);	
			
			glVertex3f(j, i, 0.0);
			glVertex3f(j+1, i, 0.0);
			glVertex3f(j+1, i+1, 0.0);
			glVertex3f(j, i+1, 0.0);
			
			
		}
	}
	
	glEnd();
	
	
	//glutSwapBuffers();
	return;
}

void draw(void){
	glClearColor(0,0.5, 0.6, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(scanUserData.gaps){
		glViewport(0,0, myGraphicsData.width/2.0, myGraphicsData.height/2.0);
		drawImage(scanUserData.gaps, scanUserData.nx, scanUserData.nx, -scanUserData.gap0, scanUserData.gap0);
	}
	if(scanUserData.bandEnergy){
		glViewport(myGraphicsData.width/2.0,0, myGraphicsData.width/2.0, myGraphicsData.height/2.0);
		drawImage(scanUserData.bandEnergy, scanUserData.nx, scanUserData.nx, scanUserData.bandMin, scanUserData.bandMax);
	}

	if(scanUserData.quasiEnergy){
		glViewport(myGraphicsData.width/2.0,myGraphicsData.height/2.0, myGraphicsData.width/2.0, myGraphicsData.height/2.0);
		drawImage(scanUserData.quasiEnergy, scanUserData.nx, scanUserData.nx, scanUserData.quasiMin, scanUserData.quasiMax);
	}

	if(scanUserData.spec){
		glViewport(0,myGraphicsData.height/2.0, myGraphicsData.width/2.0, myGraphicsData.height/2.0);
		displaySpecData(scanUserData.spec, scanUserData.numSpecVoltages, scanUserData.specMin, scanUserData.specMax, 1.0,0.0,0.0);
		//drawImage(scanUserData.quasiEnergy, scanUserData.nx, scanUserData.nx, scanUserData.quasiMin, scanUserData.quasiMax);
	}

	glutSwapBuffers();
}

void onResize(int width, int height){
	myGraphicsData.width = width;
	myGraphicsData.height = height;
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

	/*actually uses bare dispersion, not quasi dispersion*/
	calcG11((double)scanUserData.vMax, scanUserData.numSpecVoltages, gamma, scanUserData.bandEnergy, scanUserData.gaps, scanUserData.G11, scanUserData.nx);
	//calcG11((double)scanUserData.vMax, scanUserData.numSpecVoltages, gamma, scanUserData.quasiEnergy, scanUserData.gaps, scanUserData.G11, scanUserData.nx);
	cout << "finished G11\n";
	double v;
	double stepSize = 2*scanUserData.vMax/(scanUserData.numSpecVoltages-1);

	
	double n_squared = scanUserData.nx*scanUserData.nx;
	for(int i=0; i<scanUserData.numSpecVoltages; i++){
		v = -scanUserData.vMax + i*stepSize;
		scanUserData.spec[i] = -(1/(n_squared*PI))*scanUserData.G11[i].imag();
		//printf("%lf %lf\n", v, scanUserData.spec[i]);

	}	

	scanUserData.specMin = getMin(scanUserData.spec, scanUserData.numSpecVoltages);
	scanUserData.specMax = getMax(scanUserData.spec, scanUserData.numSpecVoltages);
	
	if(glutGetWindow() != myGraphicsData.mainWindow){
				glutSetWindow(myGraphicsData.mainWindow);
	}
	glutPostRedisplay();

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


