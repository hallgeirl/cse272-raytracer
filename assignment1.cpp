#include "assignment1.h"
#include <math.h>
#include <stdio.h>
#include "Miro.h"
#include "includes.h"
#include "Emissive.h"

using namespace std;

void 
makeTask1Scene()
{
    //real squarelight
    SquareLight *l = new SquareLight;
    l->setDimensions(2,2);
    l->setPosition(Vector3(0,10,0));
    l->setNormal(Vector3(0,-1,0));
    l->setUdir(Vector3(1,0,0));
    l->setWattage(100);

    g_scene->addObject(l);

    // mirror

    TriangleMesh * mirror = new TriangleMesh;
    mirror->createSingleTriangle();
    mirror->setV1(Vector3( 5, 4, -1));
    mirror->setV2(Vector3( 5, 6, -1));
    mirror->setV3(Vector3( 5, 6, 1));
    mirror->setN1(Vector3(-1, 0, 0));
    mirror->setN2(Vector3(-1, 0, 0));
    mirror->setN3(Vector3(-1, 0, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(mirror);
    t->setMaterial(new Phong(Vector3(0), Vector3(1))); 
//    g_scene->addObject(t);

	TriangleMesh * mirror2 = new TriangleMesh;
    mirror2->createSingleTriangle();
    mirror2->setV1(Vector3( 5, 4, -1));
    mirror2->setV2(Vector3( 5, 4, 1));
    mirror2->setV3(Vector3( 5, 6, 1));
    mirror2->setN1(Vector3(-1, 0, 0));
    mirror2->setN2(Vector3(-1, 0, 0));
    mirror2->setN3(Vector3(-1, 0, 0));
    
    Triangle* t2 = new Triangle;
    t2->setIndex(0);
    t2->setMesh(mirror2);
    t2->setMaterial(new Phong(Vector3(0), Vector3(1))); 
    //g_scene->addObject(t2);

	Plane* p = new Plane();
	p->setMaterial(new Phong());
	g_scene->addObject(p);

    g_scene->preCalc();
}

void a1task1()
{
	HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));
	Vector3 shadeResult(0);

	FILE *fp = stdout;

//	fp = fopen("irrad_pathtracing.dat", "w");
    long double res = 0.;
    for (long k = 0; k < 10000000; ++k)
	{
        Ray ray = ray.diffuse(hitInfo);
		Vector3 tempShadeResult;
		if (g_scene->traceScene(ray, tempShadeResult, 0))
		{
    //        cout << tempShadeResult[0] << endl;
            res += (long double)tempShadeResult[0];
		}

		if (k % 10 == 0 )
			fprintf(fp, "%ld %2.30lf\n", k, (double)(res/((long double)k+1.)));
			//fprintf(fp, "%i %3.30lf\n", k, (double)(res/((long double)k+1.)));
			//fprintf(fp, "%i %f\n", k, shadeResult[0]/((float)k+1));
	}
	fclose(fp);
	shadeResult /= TRACE_SAMPLES; 
}

void a1task2()
{
}

double sampleLightSource(const Vector3& incidentDir)
{
/*    Vector3 
    return 100*dot(incidentDir, Vector3(0, 1, 0))/(incidentDir-);*/
    return 0;
}

void a1task3()
{
    
	HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));
	Vector3 shadeResult(0);

	FILE *fp = stdout;

	fp = fopen("irrad_pathtracing.dat", "w");
    long double res = 0.;
    for (long k = 0; k < 10000000; ++k)
	{
        Ray ray = ray.diffuse(hitInfo);
		Vector3 tempShadeResult;
		if (g_scene->traceScene(ray, tempShadeResult, 0))
		{
            res += (long double)tempShadeResult[0];
		}

		if (k % 10 == 0 )
			fprintf(fp, "%ld %2.30lf\n", k, (double)(res/((long double)k+1.)));
			//fprintf(fp, "%i %3.30lf\n", k, (double)(res/((long double)k+1.)));
			//fprintf(fp, "%i %f\n", k, shadeResult[0]/((float)k+1));
	}
	fclose(fp);
	shadeResult /= TRACE_SAMPLES; 
}

void
makeCornellScene()
{
    g_image->resize(512, 512);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(2.5, 3, 3));
    g_camera->setLookAt(Vector3(2.5, 2.5, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(90);

/*    Sphere *sp = new Sphere;
    sp->setCenter(Vector3(3, 5.4, -3));
    sp->setRadius(0.2);
    sp->setMaterial(new Phong(Vector3(1,1,1)));
    g_scene->addObject(sp);*/

    // create and place a point light source
    SquareLight* light = new SquareLight;
    light->setPosition(Vector3(2.5, 4.9, -1));
    light->setDimensions(1, 1);
    light->setWattage(40);
    light->setColor(Vector3(1, 1, 1));
    light->setUdir(Vector3(1,0,1));
    g_scene->addLight(light);
    g_scene->addObject(light);

    /*PointLight* light = new PointLight;
    light->setPosition(Vector3(2.5, 4.9, -1));
    light->setWattage(40);
    light->setColor(Vector3(1, 1, 1));
    g_scene->addLight(light);
*/
    Material* material = new Phong(Vector3(1.0f));
    TriangleMesh * mesh = new TriangleMesh;
    mesh->load("models/cornell_box_1.obj");
	addMeshTrianglesToScene(mesh, new Phong(Vector3(1,1,1)));

    mesh = new TriangleMesh;
    mesh->load("models/cornell_box_2.obj");
    addMeshTrianglesToScene(mesh, new Phong(Vector3(1,0,0)));

    mesh = new TriangleMesh;
    mesh->load("models/cornell_box_3.obj");
    addMeshTrianglesToScene(mesh, new Phong(Vector3(0,1,0)));
    
    mesh = new TriangleMesh;
    mesh->load("models/cornell_box_4.obj");
    addMeshTrianglesToScene(mesh, new Phong(Vector3(1)));

    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
makeTestScene()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));
    cout << "Test scene" << endl;
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;
    g_scene->setEnvironment(autumnHDR);
    //g_scene->setEnvironment(new CheckerBoardTexture(Vector3(1), Vector3(0), 100));

    g_image->resize(2048, 2048);

    // set up the camera
    g_camera->setBGColor(Vector3(1.0f, 1.0f, 1.0f));
    //g_camera->setEye(Vector3(-5, 2, 3));
    g_camera->setEye(Vector3(0, 6, -10));
    g_camera->setLookAt(Vector3(0, 3, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(60);


    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(0, 10, 0));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(3000);
    g_scene->addLight(light);

    Sphere *sp = new Sphere;
    sp->setCenter(Vector3(0,5,0));
    sp->setRadius(1.5);
    sp->setMaterial(new Phong(Vector3(0), Vector3(0), Vector3(1), infinity, 1.5));
    g_scene->addObject(sp);
    
//    addModel("models/teapot.obj", new Phong(Vector3(1)), g_scene, Vector3(0.f));
//    addModel("models/teapot.obj", new Phong(Vector3(1)), g_scene, Vector3(0.f));

    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, 0,  30));
    floor->setV2(Vector3( 30, 0, -30));
    floor->setV3(Vector3(-30, 0, -30));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(new TexturedPhong(new CheckerBoardTexture(Vector3(1), Vector3(0), 10))); 
    g_scene->addObject(t);

    g_scene->preCalc();
}

void
A1makeSphereScene()
{
    g_image->resize(2048, 2048);
    
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-2, 1, 5));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(500);
    g_scene->addLight(light);

    Material* mat = new Lambert(Vector3(1.0f));

    Sphere *sphere = new Sphere;
    sphere->setRadius(1.5);
    sphere->setMaterial(mat);

    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, -1.5,  10));
    floor->setV2(Vector3( 10, -1.5, -10));
    floor->setV3(Vector3(-10, -1.5, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(mat); 
    g_scene->addObject(t);
    g_scene->addObject(sphere);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}

void
A1makeTeapotScene()
{
	LoadedTexture *autumnHDR = new LoadedTexture(string("gfx/autumnforrest.hdr"));

    g_image->resize(256, 256);
    g_scene->setEnvironment(autumnHDR);
    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.2f));
    g_camera->setEye(Vector3(-2, 3, 5));
    g_camera->setLookAt(Vector3(0, 1, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(45);

    // create and place a point light source
    PointLight * light = new PointLight;
    light->setPosition(Vector3(-3, 15, 3));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(500);
    g_scene->addLight(light);

    Material* mat = new Lambert(Vector3(1.0f));

    TriangleMesh * teapot = new TriangleMesh;
    teapot->load("models/teapot.obj");
    
    // create all the triangles in the triangle mesh and add to the scene
    for (int i = 0; i < teapot->numTris(); ++i)
    {
        Triangle* t = new Triangle;
        t->setIndex(i);
        t->setMesh(teapot);
        t->setMaterial(mat); 
        g_scene->addObject(t);
    }
    
    // create the floor triangle
    TriangleMesh * floor = new TriangleMesh;
    floor->createSingleTriangle();
    floor->setV1(Vector3(  0, 0,  10));
    floor->setV2(Vector3( 10, 0, -10));
    floor->setV3(Vector3(-10, 0, -10));
    floor->setN1(Vector3(0, 1, 0));
    floor->setN2(Vector3(0, 1, 0));
    floor->setN3(Vector3(0, 1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(floor);
    t->setMaterial(mat); 
    //g_scene->addObject(t);
    
    // let objects do pre-calculations if needed
    g_scene->preCalc();
}


