#include "assignment2.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include "Miro.h"
#include "includes.h"
#include "Emissive.h"
#include "Utility.h"

SquareLight* g_l;

using namespace std;


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
    t->setMaterial(new Phong(Vector3(0), Vector3(0))); 
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
    t2->setMaterial(new Phong(Vector3(0), Vector3(0))); 
    g_scene->addObject(t2);

    g_scene->preCalc();
}



sample samplePath(const Vector3& origin, const Vector3& direction)
{
    sample out;
    out.direct = true;
    HitInfo hitInfo(0, origin + Vector3(0,epsilon,0), Vector3(0,1,0));

    Ray ray(origin+Vector3(0,epsilon,0), direction); 

    //Path trace
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
            //hit reflective surface => reflect and trace again
            else 
            {
                ray = ray.reflect(hitInfo);
                out.direct = false;
                out.nrays++;
            }
        }
        //Missed the scene
        else 
        {
            out.hit = false;
            break;
        }
    }
    return out;
}

sample samplePath(const Vector3& origin)
{
    HitInfo hitInfo(0, origin + Vector3(0,epsilon,0), Vector3(0,1,0));

    //Path trace
    Ray ray = ray.diffuse(hitInfo);
    return samplePath(origin, ray.d);
}

void a2task1()
{
    cout << "Path tracing" << endl;
	HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));
	Vector3 shadeResult(0);

	ofstream fp("irrad_pathtracing.dat");
    int N = 100;

    map<double, string> outputs;

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < N; i++)
    {
        stringstream output;

        long double res = 0.;
        long nrays = 0;
        Vector3 origin(2.*(float)i/((float)N-1.)-1., 0, 0);
        #pragma omp critical
        {
            cout << "Sampling at " << origin  << endl;
        }
        long k;
        for (k = 1; k <= 10000000; ++k)
        {
            sample s = samplePath(origin);
            if (s.hit)
                res += s.value;
            nrays += s.nrays;

            //Division by 1/PI (or multiplying by PI) is neccesary because
            //E(f/p)=F=1/n*sum(f/p) and p=1/PI (distribution of rays)
            if (k % 1000 == 0)
            {
                output << (k > 0 ? "\t" : "") << PI*(double)(res/((long double)k+1.));
                if (k % 100000 == 0 && k > 0)
                    cout << k+1 << "\t" << origin.x << "\t" << PI*(double)(res/((long double)k+1.)) << "\n";
            }
        }
        outputs[origin.x] = output.str();
    }

    map<double, string>::iterator it = outputs.begin();

    while (it != outputs.end())
    {
        fp << it->second << endl;
        it++;
    }
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

        sample s = samplePath(Vector3(0,0,0));
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

void mutate(double u[3], double u_out[3])
{
    //Mutation magnitudes
    double dpos = 0.2, dtheta = 0.25, dphi = 0.1;

    u_out[0] = u[0] + (2.*dpos*frand()-dpos);
    u_out[1] = u[1] + (2.*dtheta*frand()-dtheta); 
    u_out[2] = u[2] + (2.*dphi*frand()-dphi);

    if (u_out[0] < -1) u_out[0] = -1;
    else if (u_out[0] > 1) u_out[0] = 1;

    if (u_out[1] < 0) u_out[1] += 2.*PI;
    else if (u_out[1] > 2.*PI) u_out[1] -= 2.*PI;

    if (u_out[2] < 0) u_out[2] = 0;
    else if (u_out[2] > PI/2.) u_out[2] = PI/2;
}

void a2task4()
{
    cout << "Metropolis sampling" << endl;
	HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));
	Vector3 shadeResult(0);

	ofstream fp("irrad_pathtracing.dat");
    const int N = 10;

    //Output irradiances for each point
    //There's N points distributed uniformly over the interval
    double outputs[N];

    memset(outputs, 0, sizeof(double)*N);

    double u[3];
    //Choose an initial path that gives a non-zero contribution

    double I = 0;

    Vector3 surfaceNormal(0,1,0);

    cout << "Attempting to find a starting point..." << endl;

    while (I <= 0)
    {
        //Position
        u[0] = 2.*frand() - 1;
        //Direction
        u[1] = 2.*PI*frand();       //Theta (rotation)
        u[2] = asin(sqrt(frand())); //Phi (pitch)

        sample s = samplePath(u[0], alignHemisphereToVector(surfaceNormal, u[1], u[2]));
        if (s.hit)
            I = s.value;
    }

    cout << "Found one! Contribution: " << I << endl;
    
    for (int i = 1; i < 1000000; i++)
    {
        //Mutate path. 
        double u_mutated[3];
        mutate(u, u_mutated);

        sample s = samplePath(u_mutated[0], alignHemisphereToVector(surfaceNormal, u_mutated[1], u_mutated[2]));
        if (!s.hit) s.value = 0;

        double accept = s.value / I;
        if (frand() < accept)
        {
            u[0] = u_mutated[0];
            u[1] = u_mutated[1];
            u[2] = u_mutated[2];
            I = s.value;
        }

        outputs[(int)(((u[0]+1.)/2)*N)] += I;

        for (int j = 0; j < N; j++)
            cout << outputs[j]/i << "\t";
        cout << endl;


/*        long double res = 0.;
        long nrays = 0;
        Vector3 origin(2.*(float)i/((float)N-1.)-1., 0, 0);
        #pragma omp critical
        {
            cout << "Sampling at " << origin  << endl;
        }
        long k;
        for (k = 1; k <= 10000000; ++k)
        {
            sample s = samplePath(origin);
            if (s.hit)
                res += s.value;
            nrays += s.nrays;

            //Division by 1/PI (or multiplying by PI) is neccesary because
            //E(f/p)=F=1/n*sum(f/p) and p=1/PI (distribution of rays)
            if (k % 1000 == 0)
            {
                output << (k > 0 ? "\t" : "") << PI*(double)(res/((long double)k+1.));
                if (k % 100000 == 0 && k > 0)
                    cout << k+1 << "\t" << origin.x << "\t" << PI*(double)(res/((long double)k+1.)) << "\n";
            }
        }
        outputs[origin.x] = output.str();*/
    }
}
