#include <math.h>
#include <string>
#ifdef OPENMP
#include <omp.h>
#endif
#include "Miro.h"
#include "MiroWindow.h"
#include "assignment1.h"
#include <FreeImage.h>
#include "Camera.h"
#include "Image.h"
#include "Scene.h"

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
    srand(time(0));
    //Initialize FreeImage
    FreeImage_Initialise();

#ifdef __SSE4_1__
    cout << "Using SSE" << endl;
#endif
#ifdef OPENMP
    cout << "Using OpenMP with up to " << omp_get_max_threads() << " threads." << endl;
#endif


#if ! defined(NO_GFX) and ! defined(ALTERNATIVE)
    MiroWindow miro(&argc, argv);
    cout << "Executing main rendering loop" << endl;
    miro.mainLoop();
#elif ! defined(ALTERNATIVE)
    //A1makeTeapotScene();
    makeCornellScene();
    cout << "Rendering without display" << endl;
    g_camera->setRenderer(Camera::RENDER_RAYTRACE);
    g_camera->click(g_scene, g_image);
    g_image->writePPM();
#else
/*    makeCornellScene();
    cout << "Rendering without display" << endl;
    g_camera->setRenderer(Camera::RENDER_RAYTRACE);
    g_camera->click(g_scene, g_image);
    g_image->writePPM();*/

    cout << "Alternative tasks" << endl;
	makeTask1Scene();
    a1task1();
#endif

	FreeImage_DeInitialise();

    return 0;
}

