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

using namespace std;


int
main(int argc, char*argv[])
{
//    srand(time(0));
    //Initialize FreeImage
    FreeImage_Initialise();
#ifdef __SSE4_1__
    cout << "Using SSE" << endl;
#endif
#ifdef OPENMP
    cout << "Using OpenMP with up to " << omp_get_max_threads() << " threads." << endl;
#endif

	makeTask1Scene();

    MiroWindow miro(&argc, argv);
#ifndef NO_GFX
    miro.mainLoop();
#else
    g_camera->setRenderer(Camera::RENDER_RAYTRACE);
    g_camera->click(g_scene, g_image);
    g_image->writePPM();
#endif

	FreeImage_DeInitialise();

    return 0;
}

