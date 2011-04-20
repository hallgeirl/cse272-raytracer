#include "assignment1.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <map>
#include "Miro.h"
#include "includes.h"
#include "Emissive.h"

SquareLight* g_l;

using namespace std;

struct sample
{
    double value, dist2, costheta;
    bool hit;
    bool direct;
};

typedef map<long, sample> sample_map;

struct samples
{
    sample_map X;
    double p;
    long n;
    bool isAreaSample;

    samples() { n = 0; }
};

sample samplePath()
{
    sample out;
    out.direct = true;
    HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));

    //Path trace
    Ray ray = ray.diffuse(hitInfo);
    while (true)
    {
        if (g_scene->trace(hitInfo, ray, 0, MIRO_TMAX))
        {
            //hit diffuse surface->we're done
            if (!hitInfo.material->isReflective())
            {
                double contrib = hitInfo.material->shade(ray, hitInfo, *g_scene)[0];
                out.value = contrib; //=25/PI=radiance
                out.costheta = dot(hitInfo.N, ray.d);
                out.dist2 = hitInfo.P.length2(); //x'-x = x'
                out.hit = true;

                break;
            }
            else 
            {
                ray = ray.reflect(hitInfo);
                out.direct = false;
            }
        }
        else 
        {
            out.hit = false;
            break;
        }
    }
    return out;
}

void 
makeTask1Scene()
{
    Material* m;

    //real squarelight
    SquareLight *l = new SquareLight;
    g_l = l;
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
	p->setMaterial(m=new Phong());
	g_scene->addObject(p);

    g_scene->preCalc();
}

void a1task1()
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

        //Division by 1/PI (or multiplying by PI) is neccesary because
        //E(f/p)=F=1/n*sum(f/p) and p=1/PI (distribution of rays)
		if (k % 10 == 0 )
			fprintf(fp, "%ld %2.30lf\n", k, PI*(double)(res/((long double)k+1.)));
        if (k % 1000 == 0)
            printf("%ld %2.30lf\n", k, PI*(double)(res/((long double)k+1.)));


			//fprintf(fp, "%i %3.30lf\n", k, (double)(res/((long double)k+1.)));
			//fprintf(fp, "%i %f\n", k, shadeResult[0]/((float)k+1));
	}
	fclose(fp);
}

void a1task2()
{
}

sample sampleLightSource()
{
    Vector3 lightPos = g_l->samplePhotonOrigin();

    //origin is (0,0,0)
    Vector3 l = lightPos;
    l.normalize();

    double cos_theta = dot(l, Vector3(0, 1, 0));
    
    sample out;
    out.costheta = cos_theta;
    out.dist2 = lightPos.length2();

    out.value = g_l->radiance(lightPos, l) * cos_theta*cos_theta/out.dist2;
    out.direct = true;

    return out;
}

//Evaluate the multi sample estimator using the balance heuristic
double balanceHeuristic(vector<samples> s)
{
    sample_map::iterator it;

    long N = 0;
    for (int i = 0; i < s.size(); i++)
        N += s[i].n;

    vector<double> c;

    for (int i = 0; i < s.size(); i++)
        c.push_back((double)s[i].n/(double)N);

    double res = 0;
    for (size_t i = 0; i < s.size(); i++)
    {
        /*if (i == 0) continue; //ignore contributions from direct sampling
        if (i == 1) continue; //ignore contributions from path tracing sampling*/

        for (it = s[i].X.begin(); it != s[i].X.end(); it++)
        {
            sample X = it->second;
            double p = 0;

            for (int k = 0; k < s.size(); k++)
            {
                if (s[k].isAreaSample)
                {
                    //No conversion needed
                    p += s[k].p * c[k];
                }
                else
                {
                    //convert to area probability 
                    p += (1./PI)*X.costheta * X.costheta  * c[k] / X.dist2;
                }
            }
            res += X.value / p;
        }
    }

    return res/(double)N;
}

void a1task3()
{
	Vector3 shadeResult(0);

    ofstream fp("irrad_importance.dat");

    vector<samples> allsamples;
    allsamples.push_back(samples());
    allsamples.push_back(samples());

    samples& directSamples = allsamples[0],
           & pathTraceSamples = allsamples[1];

    //we take direct samples of the light source distributed over the area
    directSamples.p = 1./4;
    directSamples.isAreaSample = true;

    //path trace samples are taken using diffuse rays, that is, with a distribution of 1/PI
    pathTraceSamples.p = 1./PI;
    pathTraceSamples.isAreaSample = false;

    //Take samples

    double indirect = 0;

    for (long k = 0; k < 10000000; ++k)
	{
		if (k % 1000 == 0 && k > 0)
        {
            double bh = balanceHeuristic(allsamples);
            fp << k << " " << bh << "\n";
            cout << k << " " << bh << endl;
        }

        sample s = samplePath();

        //we only want to importance sample the direct lighting
        if (s.hit && s.direct)
            pathTraceSamples.X[k] = s;
        else if (s.hit)
            indirect += s.value;

        pathTraceSamples.n++;

        //Sample light source
        directSamples.X[k] = sampleLightSource();
        directSamples.n++;
	}
    fp.close();
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


