#include "assignment2.h"
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
    int nrays;
    sample() { nrays = 0; } 
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
    out.nrays++;
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
                out.nrays++;
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
makeTask2Scene()
{
    Material* m;

    g_image->resize(512, 512);

    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.0f));
    g_camera->setEye(Vector3(0, -0.3, 3));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(90);
    g_scene->setBgColor(Vector3(.5,.5,.5));

    //real squarelight
    SquareLight *l = new SquareLight;
    g_l = l;
    l->setDimensions(3.5,.1);
    l->setPosition(Vector3(0,-1,0));
    l->setNormal(Vector3(0,1,0));
    l->setUdir(Vector3(1,0,0));
    l->setWattage(100);
	l->setColor(Vector3(1.f));

    g_scene->addObject(l);
	g_scene->addLight(l);

    // mirror
    TriangleMesh * mirror = new TriangleMesh;
    mirror->createSingleTriangle();
    mirror->setV1(Vector3( -2, 1, -2));
    mirror->setV2(Vector3( 2, 1, -2));
    mirror->setV3(Vector3( 2, 1, 2));
    mirror->setN1(Vector3(0, -1, 0));
    mirror->setN2(Vector3(0, -1, 0));
    mirror->setN3(Vector3(0, -1, 0));
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(mirror);
    t->setMaterial(new Phong(Vector3(0), Vector3(1))); 
    g_scene->addObject(t);

	TriangleMesh * mirror2 = new TriangleMesh;
    mirror2->createSingleTriangle();
    mirror2->setV1(Vector3( -2, 1, -2));
    mirror2->setV2(Vector3( -2, 1, 2));
    mirror2->setV3(Vector3( 2, 1, 2));
    mirror2->setN1(Vector3(0, -1, 0));
    mirror2->setN2(Vector3(0, -1, 0));
    mirror2->setN3(Vector3(0, -1, 0));

    Triangle* t2 = new Triangle;
    t2->setIndex(0);
    t2->setMesh(mirror2);
    t2->setMaterial(new Phong(Vector3(0), Vector3(1))); 
    g_scene->addObject(t2);

	// floor triangle
    TriangleMesh * square1 = new TriangleMesh;
    square1->createSingleTriangle();
    square1->setV1(Vector3( -1, 0, -1));
    square1->setV2(Vector3( 1, 0, -1));
    square1->setV3(Vector3( 1, 0, 1));
    square1->setN1(Vector3(0, 1, 0));
    square1->setN2(Vector3(0, 1, 0));
    square1->setN3(Vector3(0, 1, 0));
    
    t = new Triangle;
    t->setIndex(0);
    t->setMesh(square1);
    t->setMaterial(new Phong(Vector3(1), Vector3(0))); 
    g_scene->addObject(t);

	TriangleMesh * square2 = new TriangleMesh;
    square2->createSingleTriangle();
    square2->setV1(Vector3( -1, 0, -1));
    square2->setV2(Vector3( -1, 0, 1));
    square2->setV3(Vector3( 1, 0, 1));
    square2->setN1(Vector3(0, 1, 0));
    square2->setN2(Vector3(0, 1, 0));
    square2->setN3(Vector3(0, 1, 0));

    t2 = new Triangle;
    t2->setIndex(0);
    t2->setMesh(square2);
    t2->setMaterial(new Phong(Vector3(1), Vector3(0))); 
    g_scene->addObject(t2);

    g_scene->preCalc();
}

void a2task1()
{
    cout << "Path tracing" << endl;
	HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));
	Vector3 shadeResult(0);

	FILE *fp = stdout;

	fp = fopen("irrad_pathtracing.dat", "w");
    long double res = 0.;
    long nrays = 0;
    for (long k = 0; k < 1300000; ++k)
	{
        sample s = samplePath();
        if (s.hit)
            res += s.value;
        nrays += s.nrays;

        //Division by 1/PI (or multiplying by PI) is neccesary because
        //E(f/p)=F=1/n*sum(f/p) and p=1/PI (distribution of rays)
        if (k % 1000 == 0)
        {
			fprintf(fp, "%ld %2.30lf\n", k, PI*(double)(res/((long double)k+1.)));
            printf("%ld %2.30lf\n",  nrays, PI*(double)(res/((long double)k+1.)));
        }
	}
	fclose(fp);
}

void a2task2()
{
    cout << "Photon mapping" << endl;
	HitPoint *hp = new HitPoint;
	hp->position = Vector3(0.f);
	hp->normal = Vector3(0, 1, 0);
	hp->radius = 0.25f;

	g_scene->addHitPoint(hp);

	FILE *fp;
	fp = fopen("irrad_progphotonmapping.dat", "w");

	while (g_scene->GetPhotonsEmitted() < 100000000)
	{
		g_scene->ProgressivePhotonPass();

		printf("%ld %lf %lf %d \n", g_scene->GetPhotonsEmitted(), (double)hp->accFlux / PI / pow(hp->radius, 2) / g_scene->GetPhotonsEmitted(), hp->radius, hp->accPhotons);
		fprintf(fp, "%ld %lf %lf %d \n", g_scene->GetPhotonsEmitted(), (double)hp->accFlux / PI / pow(hp->radius, 2) / g_scene->GetPhotonsEmitted(), hp->radius, hp->accPhotons);
	}
	fclose(fp);

}

void a2task3()
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
    long nrays = 0;
    long print = -1;

    for (long k = 0; nrays < 1100000; ++k)
	{
		if (nrays >= print)
        {
            //double bh = balanceHeuristic(allsamples) + indirect/(k+1);
            double bh = 1;
            cout << nrays << " " << bh << endl;
            fp << nrays << " " << bh << "\n";
            print += 1000;
        }

        sample s = samplePath();
        nrays += s.nrays;

        //we only want to importance sample the direct lighting
        if (s.hit && s.direct)
            pathTraceSamples.X[k] = s;
        else if (s.hit)
            indirect += s.value*PI;

        pathTraceSamples.n++;

        //Sample light source
//        directSamples.X[k] = sampleLightSource();
        directSamples.n++;
        nrays += directSamples.X[k].nrays;
	}
    fp.close();
}
