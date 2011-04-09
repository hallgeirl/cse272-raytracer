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
    g_scene = new Scene;
 
// fake squarelight

	Emissive* lightMat = new Emissive;
	lightMat->setPower(25.0f);		// 100W / 4 m^2

	TriangleMesh * slight1= new TriangleMesh;
    slight1->createSingleTriangle();
    slight1->setV1(Vector3( -1, 10, -1));
    slight1->setV2(Vector3( -1, 10, 1));
    slight1->setV3(Vector3( 1, 10, 1));
    slight1->setN1(Vector3(0, -1, 0));
    slight1->setN2(Vector3(0, -1, 0));
    slight1->setN3(Vector3(0, -1, 0));
    
    Triangle* st = new Triangle;
    st->setIndex(0);
    st->setMesh(slight1);
    st->setMaterial(lightMat);
    g_scene->addObject(st);

	TriangleMesh * slight2 = new TriangleMesh;
    slight2->createSingleTriangle();
    slight2->setV1(Vector3( -1, 10, -1));
    slight2->setV2(Vector3( 1, 10, -1));
    slight2->setV3(Vector3( 1, 10, 1));
    slight2->setN1(Vector3(0, -1, 0));
    slight2->setN2(Vector3(0, -1, 0));
    slight2->setN3(Vector3(0, -1, 0));
    
    Triangle* st2 = new Triangle;
    st2->setIndex(0);
    st2->setMesh(slight2);
    st2->setMaterial(lightMat);
    g_scene->addObject(st2);

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
    t->setMaterial(new Phong(Vector3(1), Vector3(1))); 
    g_scene->addObject(t);

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
    t2->setMaterial(new Phong(Vector3(1), Vector3(1))); 
    g_scene->addObject(t2);

	Plane* p = new Plane();
	p->setMaterial(new Phong());
	g_scene->addObject(p);

    g_scene->preCalc();

	HitInfo hitInfo(0, Vector3(0), Vector3(0,1,0));
	hitInfo.material = p->getMaterial();
	Vector3 shadeResult(0);

	FILE *fp;

	fp = fopen("path_tracing_irradiance.txt", "w");

	for (int k = 0; k < 1000000; ++k)
	{
        Ray ray = ray.random(hitInfo);
		Vector3 tempShadeResult;
		if (g_scene->traceScene(ray, tempShadeResult, 0))
		{
			shadeResult += tempShadeResult;
		}
		if (k % 1000 == 0)
			fprintf(fp, "%i %f\n", k, shadeResult[0]/(k+1));
	}
	fclose(fp);
	shadeResult /= TRACE_SAMPLES; 

}

void
makeCornellScene()
{
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

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
    PointLight * light = new PointLight;
    light->setPosition(Vector3(2.5, 4.9, -1));
    light->setColor(Vector3(1, 1, 1));
    light->setWattage(40);
    g_scene->addLight(light);

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
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

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
    g_camera = new Camera;
    g_scene = new Scene;
    g_image = new Image;

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

