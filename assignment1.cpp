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

//TODO
//Do we need to divide by cos theta when gathering samples? Hmm...

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
	l->setColor(Vector3(1.f));

    g_scene->addObject(l);
	g_scene->addLight(l);

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
    t2->setMaterial(new Phong(Vector3(0), Vector3(1))); 
    g_scene->addObject(t2);

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
    long nrays = 0;
    for (long k = 0; k < 10000000; ++k)
	{
        sample s = samplePath();
        if (s.hit)
            res += s.value;
        nrays += s.nrays;

        /*Ray ray = ray.diffuse(hitInfo);
		Vector3 tempShadeResult;
		if (g_scene->traceScene(ray, tempShadeResult, 10))
		{
            res += (long double)tempShadeResult[0];
		}*/

        //Division by 1/PI (or multiplying by PI) is neccesary because
        //E(f/p)=F=1/n*sum(f/p) and p=1/PI (distribution of rays)
		if (k % 10 == 0 )
			fprintf(fp, "%ld %2.30lf\n", k, PI*(double)(res/((long double)k+1.)));
        if (k % 1000 == 0)
        {
            printf("%ld %2.30lf\n",  nrays, PI*(double)(res/((long double)k+1.)));
        }
	}
	fclose(fp);
}

void a1task2()
{
	HitPoint *hp = new HitPoint;
	hp->position = Vector3(0.f);
	hp->normal = Vector3(0, 1, 0);
	hp->radius = 0.25f;

	g_scene->addHitPoint(hp);

	int iter = 0;

	while (g_scene->GetPhotonsEmitted() < 100000000)
	{
		g_scene->ProgressivePhotonPass();

		++iter;
		printf("Iteration: %d \n", iter);
	}
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
    out.nrays ++;

    return out;
}

//Evaluate the multi sample estimator using the balance heuristic
double balanceHeuristic(vector<samples> &s)
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
                    if (s[i].isAreaSample)
                        p += s[k].p * c[k];
                    else
                        p += (1./PI)*X.costheta * X.costheta  * c[k] / X.dist2;
                }
                else
                {
                    //convert to area probability 
                    if (!s[i].isAreaSample)
                        p += s[k].p * c[k];
                    else
                        p += X.dist2 / (X.costheta*4);
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
    long nrays = 0;


    for (long k = 0; nrays < 1300000; ++k)
	{
		if (k % 1000 == 0)
        {
            double bh = balanceHeuristic(allsamples) + indirect/(k+1);
            cout << nrays << " " << bh << endl;
            fp << nrays << " " << bh << "\n";
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
        directSamples.X[k] = sampleLightSource();
        directSamples.n++;
        nrays += directSamples.X[k].nrays;
	}
    fp.close();
}
