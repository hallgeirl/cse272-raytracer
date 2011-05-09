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

// initialize measurement points from A to B 
typedef std::vector<HitPoint*> HitPoints;
HitPoints m_hitpoints;

using namespace std;


void 
makeTask2Scene()
{
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

	// initialize measurement points from A to B 
	// might be issue with radii large than floor
	HitPoint *hp;
	float delta = 0.02;
	
	for (float i = -0.99; i < 1; i+=delta)
	{
		hp = new HitPoint;
		hp->position = Vector3(i, 0.f, 0.f);
		hp->normal = Vector3(0, 1, 0);
		hp->radius = 0.25f;

		// for task 2
		g_scene->addHitPoint(hp);
		// for task 3
		m_hitpoints.push_back(hp);
	}
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

bool UpdateMeasurementPoints(const Vector3& pos, const Vector3& power)
{
	bool hit = false;

	for (int n = 0; n < m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];
		float d = sqrt(pow(pos.x - hp->position.x, 2) +
			pow(pos.y - hp->position.y, 2) +
			pow(pos.z - hp->position.z, 2));

		if (d <= hp->radius)
		{
			float delta = (hp->accPhotons + PHOTON_ALPHA)/(hp->accPhotons + 1);
			hp->radius *= sqrt(delta);
			hp->accPhotons++;
			//need to update radius and flux
			hp->accFlux = (hp->accFlux + power.x) * delta;

			// can hit multiple measurement points
			hit = true;
		}
	}
	return hit;
}

bool SamplePhotonPath(const Ray& path, const Vector3& power)
{
    HitInfo hitInfo(0, path.o + Vector3(0,epsilon,0), Vector3(0,1,0));

	Ray ray(path.o + Vector3(0,epsilon,0), path.d);

    while (true)
    {
        if (g_scene->trace(hitInfo, ray, 0, MIRO_TMAX))
        {
            //hit diffuse surface->we're done
            if (!hitInfo.material->isReflective())
            {
                return UpdateMeasurementPoints(hitInfo.P, power);
            }
            //hit reflective surface => reflect and trace again
            else 
            {
                ray = ray.reflect(hitInfo);
            }
        }
        //Missed the scene
        else 
        {
            return false;
        }
    }
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

	ofstream fp("irrad_progphotonmapping.dat");

	while (g_scene->GetPhotonsEmitted() < 1000000)
	{
		g_scene->ProgressivePhotonPass();
		//printf("%ld %lf %lf %d \n", g_scene->GetPhotonsEmitted(), (double)hp->accFlux / PI / pow(hp->radius, 2) / g_scene->GetPhotonsEmitted(), hp->radius, hp->accPhotons);
	}

	for (int n = 0; n < g_scene->hitpoints()->size(); ++n)
	{
		HitPoint *hp = (*g_scene->hitpoints())[n];
		
		fp << (double)hp->accFlux / PI / pow(hp->radius, 2) / g_scene->GetPhotonsEmitted() << "\t" << hp->position.x << endl;
	}
	fp.close();

}

void a2task3()
{
    ofstream fp("irrad_adaptiveppm.dat");

	//find starting good path
    Vector3 power = g_l->color() * g_l->wattage(); 
	Ray goodPath;
	int m_photonsEmitted;
	float prev_di = 1;
	float prev_ai = 0;
	float mutated = 1;
	float accepted = 0;

	do
	{
        goodPath.o = g_l->samplePhotonDirection();
        goodPath.d = g_l->samplePhotonOrigin();
	} while (!SamplePhotonPath(goodPath, power));

	for (m_photonsEmitted = 0; m_photonsEmitted < 1000000; m_photonsEmitted++)
    {
        //Test new random photon
        Ray path(g_l->samplePhotonDirection(), g_l->samplePhotonOrigin());
		if (SamplePhotonPath(path, power))
		{
			goodPath = path;
			continue;
		}

		//Mutate path
		float di = prev_di + (1.f / mutated) * (prev_ai - 0.234);
		float dui = ((2 * frand() - 1) > 0 ? 1 : -1) * pow(frand(), di+1); 
		++mutated;
		path.d = goodPath.d + dui;
		prev_di = di;
		prev_ai += ((float)accepted/(float)mutated - 0.234) / (float)mutated;
		
		if (SamplePhotonPath(path, power))
		{
			goodPath = path;
			++accepted;
			continue;
		}
    }

	for (int n = 0; n < g_scene->hitpoints()->size(); ++n)
	{
		HitPoint *hp = (*g_scene->hitpoints())[n];

		float result = (double)hp->accFlux / PI / pow(hp->radius, 2) * ((float)hp->accPhotons / (float)m_photonsEmitted);
		
		fp << result << "\t" << hp->position.x << endl;
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
