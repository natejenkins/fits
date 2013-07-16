#include "natestm.h"
#include <GL/glui.h>

void draw(void);
void drawTopo(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);
void drawSpectro(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);
void drawEcut(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);

void gluiCallback(int id);

graphicsData myGraphicsData;

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


int init_glui(int argc, char* argv[]){

	myGraphicsData.width = myGraphicsData.height = 700;
	glutInit(&argc, argv);
	glutInitWindowSize(myGraphicsData.width,myGraphicsData.height);
	glutInitDisplayMode(GLUT_DOUBLE);
	myGraphicsData.mainWindow = glutCreateWindow("natemodel");
	glutDisplayFunc(draw);
	glutReshapeFunc(onResize);
	
	initGL();

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
	myGraphicsData.glui->add_button_to_panel(myGraphicsData.topoRollout, "Calc Weights", -1, onCalcWeights);

	myGraphicsData.glui->sync_live();
	

	myGraphicsData.autoScaleTopo = 1;
	myGraphicsData.autoScaleSpec = 1;

	GLUI_Master.set_glutIdleFunc(idle);

	glutPostRedisplay();
	glutMainLoop();

	return(0);
}