#include "natestm.h"
#include <GL/glui.h>

void draw(void);
void drawTopo(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);
void drawSpectro(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);
void drawEcut(const graphicsData& myGraphicsData, const ScanUserData& scanUserData);

void gluiCallback(int id);
int init_glui(int argc, char* argv[]);
void post_redisplay(void);