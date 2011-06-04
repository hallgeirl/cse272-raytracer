#include <math.h>
#include <string>
#ifdef OPENMP
#include <omp.h>
#endif
#include "Miro.h"
#include "MiroWindow.h"
#include "assignment3.h"
#include <FreeImage.h>
#include "Camera.h"
#include "Image.h"
#ifdef PHOTON_MAPPING
#include "PScene.h"
#else
#include "Scene.h"
#endif

using namespace std;

void setup()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;
}


int
main(int argc, char*argv[])
{
    setup();
//    srand(time(0));
    //Initialize FreeImage
#ifndef NO_FREEIMAGE
    FreeImage_Initialise();
#endif

#ifdef __SSE4_1__
    cout << "Using SSE" << endl;
#endif
#ifdef OPENMP
    cout << "Using OpenMP with up to " << omp_get_max_threads() << " threads." << endl;
#endif
#ifdef LINUX
    srand48(time(0));
#endif
//mode = 0: Create opengl window and everything
//mode = 1: Render scenes without any GUI
//mode = 2: Other things
int mode = 0;

#if defined(NO_GFX)
    mode = 1;
#endif
#if defined(ALTERNATIVE)
    mode = 2;
#endif

cout << "Mode: " << mode << endl;

makeTask3Scene();

if (mode == 0)
{
    cout << "Executing main rendering loop" << endl;
    MiroWindow miro(&argc, argv);
    miro.mainLoop();
}
else if (mode == 1)
{
    //A1makeTeapotScene();
    cout << "Rendering without display" << endl;
    g_camera->setRenderer(Camera::RENDER_RAYTRACE);
    g_camera->click(g_scene, g_image);
    g_image->writePPM();
}
else
{
    cout << "Alternative tasks" << endl;
//    a3task1();
//    a3task2();
//    a3task3();
}

#ifndef NO_FREEIMAGE
	FreeImage_DeInitialise();
#endif
    return 0;
}

