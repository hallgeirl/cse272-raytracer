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
using std::min;
using std::max;


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
    mirror2->setV2(Vector3( 2, 1, 2));
    mirror2->setV3(Vector3( -2, 1, 2));
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
    square1->setV2(Vector3( 1, 0, 1));
    square1->setV3(Vector3( 1, 0, -1));
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
	float delta = 0.016;
	
	for (float i = -0.8; i < 0.8; i+=delta)
	{
		hp = new HitPoint;
		hp->position = Vector3(i, 0.f, 0.f);
		hp->normal = Vector3(0, 1, 0);
		hp->radius = 0.2f;

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
			//wait to update radius and flux
			hp->newPhotons++;
			hp->newFlux += power.x;

			// can hit multiple measurement points	
			hit = true;
		}
	}
	return hit;
}

void UpdatePhotonStats()
{
	for (int n = 0; n < m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];
		
		// only adding a ratio of the newly added photons
		float delta = (hp->accPhotons + PHOTON_ALPHA * hp->newPhotons)/(hp->accPhotons + hp->newPhotons);
		hp->radius *= sqrt(delta);
		hp->accPhotons += (int)(PHOTON_ALPHA * hp->newPhotons);
		
		// not sure about this flux acc, or about calculating the irradiance
		hp->accFlux = ( hp->accFlux + hp->newFlux) * delta;	

		// reset new values
		hp->newPhotons = 0;
		hp->newFlux = 0.f;
	}
}

bool SamplePhotonPath(const Ray& path, const Vector3& power)
{
    HitInfo hitInfo(0, path.o, Vector3(0,1,0));

	Ray ray(path.o, path.d);

    while (true)
    {
        if (g_scene->trace(hitInfo, ray, 0, MIRO_TMAX))
        {
			//return if hit triangle backface
			if (dot(ray.d, hitInfo.N) > 0)
				return false;

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

	ofstream fp("pathtracing_irrad.dat");
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
        for (k = 1; k <= 100000000; ++k)
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
		
		fp << (double)hp->accFlux / PI / pow(hp->radius, 2) / (float)g_scene->GetPhotonsEmitted() << "\t" << hp->position.x << endl;
	}
	fp.close();

}

void a2task3()
{
    ofstream fp("adaptiveppm_irrad.dat");

	//find starting good path
    Vector3 power = g_l->color() * g_l->wattage(); 
	Ray goodPath;
	int m_photonsEmitted;
	float prev_di = 1;
	float prev_ai = 0;
	float mutated = 1;
	float accepted = 0;
	float uniform = 0;

	do
	{
        goodPath.o = g_l->samplePhotonOrigin();
        goodPath.d = g_l->samplePhotonDirection();
	} while (!SamplePhotonPath(goodPath, power));

	for (m_photonsEmitted = 0; m_photonsEmitted < 1000000; m_photonsEmitted++)
    {
		if ((m_photonsEmitted + 1)% 100000 == 0)
		{
			UpdatePhotonStats();
			printf("Update Photon Stats %f percent \n", (float)m_photonsEmitted/(float)1000000);
		} 

        //Test new random photon
        Ray path(g_l->samplePhotonOrigin(), g_l->samplePhotonDirection());
		if (SamplePhotonPath(path, power))
		{
			goodPath = path;
			++uniform;
			continue;
		}

		//Mutate path
		float di = prev_di + (1.f / mutated) * (prev_ai - 0.234);
		float dui = ((2 * frand() - 1) > 0 ? 1 : -1) * pow(frand(), (1/di)+1); 
		++mutated;

		//faking mutation for now
		path.d = goodPath.d + dui;
		path.o.x = max<float>(min<float>(goodPath.o.x + dui, 1.75), -1.75);
		path.d.normalize();

		prev_di = di;
		prev_ai = accepted/mutated;

		if (SamplePhotonPath(path, power))
		{
			goodPath = path;
			++accepted;
			continue;
		}

		SamplePhotonPath(goodPath, power);		
    }

	for (int n = 0; n < g_scene->hitpoints()->size(); ++n)
	{
		HitPoint *hp = (*g_scene->hitpoints())[n];

		float result = (double)hp->accFlux / PI / pow(hp->radius, 2) / (float)m_photonsEmitted * (uniform / (float)m_photonsEmitted);
		
		fp << result << "\t" << hp->position.x << endl;
	}

    fp.close();
}

double mutate_value(double s1, double s2)
{
//    double dv = (s2-s1)*frand()+s1;
    double dv = s2*exp(-log(s2/s1)*frand());

    if (frand() < 0.5)
        return -dv;
    else
        return dv;

}

void mutate_path(const path &p0, path &p1)
{
    //Mutation magnitudes
    double dpos = 1, dtheta = .125, dphi = .125;
    double max_mutation = 1./64.;
    double min_mutation = 1./2048.;

    if (frand() < p_large)
    {
        p1.u[0] = 2.*frand() - 1.;
        p1.u[1] = 2.*frand()*PI;
        p1.u[2] = asin(sqrt(frand()));
    }
    else
    {
        bool mutated = false;
        while (!mutated)
        {
            if (frand() < p_pos)
            {
                p1.u[0] = p0.u[0] + mutate_value(min_mutation, max_mutation) * dpos;
                if (p1.u[0] > 1) p1.u[0] -= 2;
                else if (p1.u[0] < -1) p1.u[0] += 2;
                mutated = true;
            }
            else
            {
                p1.u[0] = p0.u[0];
            }

            if (frand() < p_angle)
            {
                p1.u[1] = p0.u[1] + mutate_value(min_mutation, max_mutation)*dtheta;
                if (p1.u[1] < 0) p1.u[1] += 2.*PI;
                else if (p1.u[1] > 2.*PI) p1.u[1] -= 2.*PI;

                p1.u[2] = p0.u[2] + mutate_value(min_mutation, max_mutation)*dphi;
                if (p1.u[2] > PI/2.) p1.u[2] -= PI/2.;
                else if (p1.u[2] < 0) p1.u[2] += PI/2.;
                mutated = true;
            }
            else
            {
                p1.u[1] = p0.u[1];
                p1.u[2] = p0.u[2];
            }
        }
    }
}

void a2task4()
{
    cout << "Metropolis sampling" << endl;
	HitInfo hitInfo(0, Vector3(0, epsilon, 0), Vector3(0,1,0));
	Vector3 shadeResult(0);

    long double ptracing_output[100];
    double ptracing_b = 0;

    //Get path tracing results from last iteration
    {
        cout << "Loading path tracing results..." << endl;
        char tmpbuf[100000];
        ifstream fp_tmp("pathtracing_irrad.dat");
        if (!fp_tmp.is_open())
            memset(ptracing_output,0,sizeof(long double)*100);
        else
        {
            int i = 0;

            while (!fp_tmp.getline(tmpbuf,100000).eof())
            {
                stringstream ss(tmpbuf);
                while (ss >> ptracing_output[i]) {}
                ptracing_b += (double)ptracing_output[i];
                i++;
            }
            fp_tmp.close();
        }
        ptracing_b /= 100;
    }

	ofstream fp_err("metropolis_msq.dat");
    const int N = 100; //Number of points
    const int Nseeds = 10000000;
    const int Nsamples = 100000000;

    //Output irradiances for each point
    //There's N points distributed uniformly over the interval
    long double outputs[N];
    long double error[N];

    stringstream output_strings[N];
    stringstream error_strings[N];

    memset(outputs, 0, sizeof(long double)*N);
    memset(error, 0, sizeof(long double)*N);

    Vector3 surfaceNormal(0,1,0);

    cout << "Generating path seeds..." << endl;

    path p0;
    p0.I = 0;
    long double b = 0;

    //Generate path seeds
    for (int i = 0; i < Nseeds; i++)
    {
        double u_tmp[3];
        //Position
        u_tmp[0] = 2.*frand() - 1;

        //Direction
        u_tmp[1] = 2.*PI*frand();       //Theta (rotation)
        u_tmp[2] = asin(sqrt(frand())); //Phi (pitch)
        Vector3 dir = alignHemisphereToVector(surfaceNormal, u_tmp[1], u_tmp[2]);

        sample s = samplePath(Vector3(u_tmp[0], 0, 0), dir);
    
        if (s.hit && s.value > 0)
        {
            b += s.value*PI;
            if (p0.I < s.value)
            {
                p0.I = s.value;
                memcpy(p0.u, u_tmp, 3*sizeof(double));
            }
        }
    }

    b /= Nseeds;

    cout << "b=" << b << ", b_ptracing=" << ptracing_b << ", error=" << b/ptracing_b << endl;
    cout << "p0=" << p0.I << "; [" << p0.u[0] << ", " << p0.u[1] << ", " << p0.u[2] << "]" << endl;

    for (int i = 1; i <= Nsamples; i++)
    {
        //Mutate path. 
        path p1;
        mutate_path(p0, p1);

        sample s = samplePath(Vector3(p1.u[0],0,0), alignHemisphereToVector(surfaceNormal, p1.u[1], p1.u[2]));
        
        if (!s.hit) s.value = 0;
        p1.I = s.value;
        
        double accept;
        int x1 = (int)(((p1.u[0]+1.)/2.)*N);
        int x0 = (int)(((p0.u[0]+1.)/2.)*N);
        if (x0 == N) x0--;
        if (x1 == N) x1--;

        if (p1.u[0] >= -1 && p1.u[0] <= 1 && p1.u[2] >= 0 && p1.u[2] <= PI/2.)
        {
            accept = std::min(p1.I / p0.I, 1.);
            outputs[x1] += accept*b;
        }
        else
        {
            accept = 0;
        }

        outputs[x0] += (1.-accept)*b;

        if (frand() < accept)
        {
            path_copy(p0, p1);
        }

        double msq = 0;
        if (i % 1000 == 0)
        {

            for (int j = 0; j < N; j++)
            {
                output_strings[j] << outputs[j]/(i+1)*(long double)N << " ";//(j==N-1?"":"\t");
                long double err = (outputs[j]/(i+1)*(long double)N-ptracing_output[j]);
                long double err_sq = err*err;
                error[j] = err;

                msq += err_sq;
            }
            msq /= N;

            fp_err << i << " " << msq << "\n";
        }

        if (i % 10000 == 0)
        {
            printf("Iteration %d\n", i);
            cout << "Error: " << endl;

            for (int j = 0; j < N; j++)
            {
                long double err = error[j];
                long double err_sq = err*err;
                if (i % 10000 == 0)
                {
                     printf("%6.3Lf%s", err, j==N-1?"":"\t");
                }
            }

            cout << "\nMean square error: " << msq << "\n";
            cout << "Error in b: " << b - ptracing_b << endl;
        }

        if (i == 1000000 || i == 100000000)
        {
            stringstream ss;
            ss << "metropolis_error_" << i << ".dat";
            
            ofstream f_errorgraph(ss.str().c_str());
            for (int j = 0; j < N; j++)
            {
                f_errorgraph << j+1 << " " << error[j] << "\n";
            }
        }
    }

    {
        ofstream fp("metropolis_irrad.dat");
        for (int i = 0; i < N; i++)
        {
            fp << output_strings[i].str() << "\n";
        }
        fp.close();
    }

    {
        ofstream fp("metropolis_b.dat");
        fp << "b=" << b << ", error (vs. pathtracing)=" << b-ptracing_b << "\n";
        fp.close();
    }
}
