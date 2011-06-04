#include "assignment3.h"
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
using std::min;
using std::max;

void BuildSquare(const Vector3& min, const Vector3& max, const Vector3& normal, const Material* mat)
{
    TriangleMesh * square1 = new TriangleMesh;
    square1->createSingleTriangle();
    if (min.z != max.z)
    {
        square1->setV1(Vector3( min.x, min.y, min.z));
        square1->setV2(Vector3( max.x, max.y, max.z));
        square1->setV3(Vector3( max.x, max.y, min.z));
    }
    else
    {
        square1->setV1(Vector3( min.x, min.y, min.z));
        square1->setV2(Vector3( max.x, max.y, max.z));
        square1->setV3(Vector3( max.x, min.y, min.z));
    }
    square1->setN1(normal);
    square1->setN2(normal);
    square1->setN3(normal);
    
    Triangle* t = new Triangle;
    t->setIndex(0);
    t->setMesh(square1);
    t->setMaterial(mat);
    g_scene->addObject(t);

	TriangleMesh * square2 = new TriangleMesh;
    square2->createSingleTriangle();
    if (min.z != max.z)
    {
        square2->setV1(Vector3( min.x, min.y, min.z));
        square2->setV2(Vector3( min.x, min.y, max.z));
        square2->setV3(Vector3( max.x, max.y, max.z));
    }
    else
    {
        square2->setV1(Vector3( min.x, min.y, min.z));
        square2->setV2(Vector3( min.x, max.y, max.z));
        square2->setV3(Vector3( max.x, max.y, max.z));
    }
    square2->setN1(normal);
    square2->setN2(normal);
    square2->setN3(normal);

    Triangle* t2 = new Triangle;
    t2->setIndex(0);
    t2->setMesh(square2);
    t2->setMaterial(mat); 
    g_scene->addObject(t2);
}

void 
makeTask3Scene()
{
    g_image->resize(W, H);

    // set up the camera
    g_camera->setBGColor(Vector3(0.0f, 0.0f, 0.0f));
    g_camera->setEye(Vector3(0, 0, -2));
    g_camera->setLookAt(Vector3(0, 0, 0));
    g_camera->setUp(Vector3(0, 1, 0));
    g_camera->setFOV(90);
    g_scene->setBgColor(Vector3(0));

    //real squarelight
    SquareLight *l = new SquareLight;
    g_l = l;
    l->setDimensions(20,2);
    l->setPosition(Vector3(0,2,0));
    l->setNormal(Vector3(0,-1,0));
    l->setUdir(Vector3(1,0,0));
    l->setWattage(1000);
    l->setColor(Vector3(1.f));

    g_scene->addObject(l);
	g_scene->addLight(l);
    
    double e = 0.0005;
	// Floor
	BuildSquare(Vector3(-1,-1,-1), Vector3(1,-1,1), Vector3(0,1,0), new Phong(Vector3(0.75)));
	
	// Back Wall
	BuildSquare(Vector3(-1,-1,1), Vector3(1,1,1), Vector3(0,0,-1), new Phong(Vector3(0.75)));


	// Left Wall
	BuildSquare(Vector3(-1-e,-1,-1), Vector3(-1,1,1), Vector3(1,0,0), new Phong(Vector3(0), Vector3(0.75)));
#ifdef HACKER3
    //Hackerpoint: Red wall
	BuildSquare(Vector3(1+e,-1,-1), Vector3(1,1,1), Vector3(-1,0,0), new Phong(Vector3(0.75,0,0)));
#else
	// Right Wall
	BuildSquare(Vector3(1+e,-1,-1), Vector3(1,1,1), Vector3(-1,0,0), new Phong(Vector3(0.75)));
#endif
	// Ceiling
	BuildSquare(Vector3(0.05,1,-1), Vector3(1,1,1), Vector3(0,-1,0), new Phong(0.75));
	BuildSquare(Vector3(-0.95,1,-1), Vector3(0,1,1), Vector3(0,-1,0), new Phong(0.75));

    //Hackerpoint: Sphere
#ifdef HACKER2
    Sphere* sp = new Sphere();
    sp->setCenter(Vector3(0,0,0));
    sp->setRadius(0.5);
    sp->setMaterial(new Phong(Vector3(0), Vector3(0), Vector3(1), 1, 1.5));
    g_scene->addObject(sp);
#endif
    

    g_scene->preCalc();
}

sample samplePath(const path& p, int w, int h)
{
    sample out;
    Ray ray = g_camera->eyeRay((int)(p.u[0]*(double)W), (int)(p.u[1]*(double)H), W, H, false);

    int depth = PATH_LENGTH;

    int pathpos = 2; //next random number to be used is at index 2

    HitInfo hitInfo;
    PointLight* l;
    Vector3 contribution(1./PI,1./PI,1./PI); //Current contribution (decreases if surface we hit has reflectance < 1)

    bool specular = false;
    bool diffuse = false;

    static int nrays = 0;
    nrays++;

    //Path trace
    while (depth > 0 && !out.hit)
    {
        if (g_scene->trace(hitInfo, ray, 0, MIRO_TMAX))
        {
            //Did we hit the light?
            l = dynamic_cast<PointLight*>(hitInfo.object);
            if (l != NULL)
            {
                out.value = contribution*l->radiance(hitInfo.P, ray.d);
                out.hit = true;
            }
            //hit reflective surface => reflect and trace again
            else if (hitInfo.material->isReflective())
            {
                ray = ray.reflect(hitInfo);
                contribution = contribution * hitInfo.material->getReflection();
            }
            //Hit a refractive surface?
            else if (hitInfo.material->isRefractive())
            {
                //Push the hit point inside (or outside if on the way out) the refractive object
                hitInfo.P += ray.d*epsilon*2.;
                ray = ray.refract(hitInfo);
                contribution = contribution * hitInfo.material->getRefraction();
            }
            //Did we hit a diffuse object?
            else
            {
                contribution = contribution * hitInfo.material->getDiffuse();
                ray = ray.diffuse(hitInfo, p.u[pathpos], p.u[pathpos+1]);
                pathpos += 2;
            }
            depth--;
        }
        //Missed the scene
        else 
        {
            out.hit = false;
            out.value = Vector3(0);
            break;
        }
    }
    return out;
}

float MutatePath(const float MutationSize)
{
	return ((2 * frand() - 1) > 0 ? 1 : -1) * pow(frand(), 1.f/MutationSize+1);
}

float ApplyDeltaRange(const float delta, float value, const float x1, const float x2)
{
	float range = x2 - x1;
    float delta_ = delta;
    /*if (delta < -range) delta_ = -range;
    if (delta > range) delta_ = range;*/
    float result = value + delta_;
    
	if (result < x1)
		result += range;
	else if ( result > x2)
		result -= range;

	return result;
}

bool UpdateMeasurementPoints(const Vector3& pos, const Vector3& power)
{
	bool hit = false;

	for (int n = 0; n <  g_scene->hitpoints()->size(); ++n)
	{
		HitPoint *hp = (*g_scene->hitpoints())[n];
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
    cout << "Alphas:" << endl;
	for (int n = 0; n < g_scene->hitpoints()->size(); ++n)
	{
		HitPoint *hp = (*g_scene->hitpoints())[n];
        double f_alpha = (long double)hp->accPhotons*5e-6;

        float alpha = PHOTON_ALPHA + (1.-PHOTON_ALPHA)*(1.-exp(-f_alpha));
        cout << alpha << "\t";

        // Set scaling factor for next photon pass
        float A = PI * pow(hp->radius, 2);
        float A1 = A;

        CircleSegment(Vector3(-1,0,-1), Vector3(0,0,1), hp->radius, hp->position, A1); 
        CircleSegment(Vector3(1,0,-1), Vector3(0,0,1), hp->radius, hp->position, A1);
        (*g_scene->hitpoints())[n]->scaling = A1/A;
		
		// only adding a ratio of the newly added photons
		float delta = (hp->accPhotons + alpha * hp->newPhotons)/(hp->accPhotons + hp->newPhotons);
		hp->radius *= sqrt(delta);
		hp->accPhotons += (int)(alpha * hp->newPhotons);
		
		// not sure about this flux acc, or about calculating the irradiance
		hp->accFlux = ( hp->accFlux + hp->newFlux/hp->scaling) * delta;	

		// reset new values
		hp->newPhotons = 0;
		hp->newFlux = 0.f;
	}
    cout << endl;
}

//void PrintPhotonStats(ofstream& fp, const float photonsEmitted, const float uniform)
void PrintPhotonStats(stringstream ss[100], const long double photonsEmitted, const long double uniform)
{
	for (int n = 0; n <  g_scene->hitpoints()->size(); ++n)
	{
		HitPoint *hp = (*g_scene->hitpoints())[n];

		long double A = PI * pow(hp->radius, 2);

		long double result = hp->accFlux / A / photonsEmitted * (uniform / photonsEmitted);
		
		ss[n] << result << "\t";
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

void a3task2()
{
    long double ptracing_output[100];
    ofstream msq_fp("adaptiveppm_msq.dat");

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
                i++;
            }
            fp_tmp.close();
        }
    }

	//find starting good path
    Vector3 power = g_l->color() * g_l->wattage(); 
	Ray goodPath;
	int m_photonsEmitted;
	float prev_di = 1;
	long mutated = 1;
	long accepted = 0;
	long uniform = 0;

    int Nphotons = 100000000;
    const int N = 100;
    stringstream ss_out[N];
	do
	{
        goodPath.o = g_l->samplePhotonOrigin();
        goodPath.d = g_l->samplePhotonDirection();
	} while (!SamplePhotonPath(goodPath, power));

    long double msq = 0;
	for (m_photonsEmitted = 0; m_photonsEmitted < Nphotons; m_photonsEmitted++)
    {

		if (m_photonsEmitted > 0 && m_photonsEmitted % 1000000 == 0)
		{
			UpdatePhotonStats();
            msq = 0;
            for (int i = 0; i < N; i++)
            {
                HitPoint* hp = g_scene->hitpoints()->at(i);

                long double A = PI * pow(hp->radius, 2);
                long double result = hp->accFlux / A / (long double)m_photonsEmitted * ((long double)uniform / (long double)m_photonsEmitted);

                long double err = result-ptracing_output[i];
                long double err_sq = err*err;

                msq += err_sq;
            }
            msq /= N;
            msq_fp << m_photonsEmitted << " " << msq << endl;
            cout << "Pass " << m_photonsEmitted / 100000 << " of " << Nphotons/100000 << ", msq=" << msq << endl;
		}
        if (m_photonsEmitted > 0 && m_photonsEmitted % 100000 == 0)
        {
			PrintPhotonStats(ss_out, m_photonsEmitted, (long double)uniform);
        }

        //Test random photon path
        Ray path(g_l->samplePhotonOrigin(), g_l->samplePhotonDirection());
		if (SamplePhotonPath(path, power))
		{
			goodPath = path;
			++uniform;
			continue;
		}

		//Mutatation size
		long double di = prev_di + (1. / (long double)mutated) * ((long double)accepted/(long double)mutated - 0.234);

		// Convert to spherical coords (theta phi reversed)
		float phi = acos(goodPath.d.z);
		float theta = atan2(goodPath.d.y, goodPath.d.x);

		// add mutation and convert back to cartesian coords
		path.d = alignHemisphereToVector(Vector3(0,1,0), ApplyDeltaRange(MutatePath(di),theta, 0.f, 2.*PI), ApplyDeltaRange(MutatePath(di), phi, 0, PI/2.)); 
        float mut = MutatePath(di); 
		path.o.x = ApplyDeltaRange(mut, goodPath.o.x, -1.75, 1.75);
		path.o.z = ApplyDeltaRange(MutatePath(di), goodPath.o.z, -0.05, 0.05);

		++mutated;
		prev_di = di;

		// Test mutated photon path
		if (SamplePhotonPath(path, power))
		{
			goodPath = path;
			++accepted;
			continue;
		}
		// Reuse good path
		SamplePhotonPath(goodPath, power);		
    }

	UpdatePhotonStats();
	PrintPhotonStats(ss_out, m_photonsEmitted, uniform);

    ofstream fp("adaptiveppm_irrad.dat");
    for (int i = 0; i < N; i++)
        fp << ss_out[i].str().c_str() << endl;
    fp.close();
}

double mutate_value(double s1, double s2)
{
    double dv = s2*exp(-log(s2/s1)*frand());

    if (frand() < 0.5)
        return -dv;
    else
        return dv;
}

void mutate_path(const path &p0, path &p1)
{
    //Mutation magnitudes
//    double dpos = 1, dtheta = .125, dphi = .125;
    double dpos = 1, dtheta = .5, dphi = .5;
    double max_mutation = 1./64.;
    double min_mutation = 1./1024.;

    if (frand() < p_large)
    {
        p1.init_random();
    }
    else
    {
        //Mutate position
        p1.u[0] = p0.u[0] + mutate_value(min_mutation, max_mutation) * dpos;
        if (p1.u[0] >= 1.) p1.u[0] -= 1;
        else if (p1.u[0] < 0) p1.u[0] += 1;

        p1.u[1] = p0.u[1] + mutate_value(min_mutation, max_mutation) * dpos;
        if (p1.u[1] > 1) p1.u[1] -= 1;
        else if (p1.u[1] < 0) p1.u[1] += 1;

        //Mutate directions
        for (int i = 0; i < PATH_LENGTH; i++)
        {
            int j = 2+i*2;
            p1.u[j] = p0.u[j] + mutate_value(min_mutation, max_mutation)*dtheta;
            if (p1.u[j] < 0) p1.u[j] += 1;
            else if (p1.u[j] > 1) p1.u[j] -= 1;

            p1.u[j+1] = p0.u[j+1] + mutate_value(min_mutation, max_mutation)*dtheta;
            if (p1.u[j+1] < 0) p1.u[j+1] += 1;
            else if (p1.u[j+1] > 1) p1.u[j+1] -= 1;
        }

    }
}

void a3task3()
{
    const int Nseeds = 1000000;
    const long Nsamples = 100000000;
    //Record the error after every this many samples
    const int error_interval = 100000;

    Vector3 img[W][H];
    memset(img, 0, W*H*sizeof(Vector3));

    cout << "Metropolis sampling" << endl;
    cout << Nsamples / (W*H) << " samples per pixel." << endl;

#ifdef HACKER2
    const char* version = "sphere";
#elif defined (HACKER3)
    const char* version = "red";
#else
    const char* version = "gray";
#endif

    vector<vector<Vector3> > ptracing_results;

    {
        //Load image from task 1
        ifstream ptracing;
        char filename[100];
        sprintf(filename, "pathtracing_%s.raw\0", version);
        ptracing.open(filename, ios::binary);

        int w_pt, h_pt;
        long double b_pt = 0;
        ptracing.read((char*)&w_pt, 4);
        ptracing.read((char*)&h_pt, 4);
        cout << "Loading path tracing results [width=" << w_pt << ", height " << h_pt << "]" << endl;
        for (int i = 0; i < h_pt; i++)
        {
            ptracing_results.push_back(vector<Vector3>());
            for (int j = 0; j < w_pt; j++)
            {
                Vector3 pix(0);
                ptracing.read((char*)&pix.x, sizeof(float));
                ptracing.read((char*)&pix.y, sizeof(float));
                ptracing.read((char*)&pix.z, sizeof(float));
                b_pt += pix.average();
                ptracing_results[i].push_back(pix);
            }
        }
        b_pt /= (long double)(w_pt*h_pt);
        cout << "Done reading path tracing results. b = " << b_pt << endl;
    }
    cout << "Generating path seeds..." << endl;

    long double b = 0;

    path p_init;
    p_init.I = 0;

    ofstream msq_out;

    {
        char filename[100];
        sprintf(filename, "metropolis_%s_msq.dat\0", version);
        msq_out.open(filename);
    }

    //Generate path seeds
    for (int i = 1; i <= Nseeds; i++)
    {
        path p_tmp;

        p_tmp.init_random();

        sample s = samplePath(p_tmp, W, H);

        if (s.hit)
        {
            long double I = s.value.average();
            b += I*PI;

            //is it a better path?
            if (p_init.I < I)
                path_copy(p_init, p_tmp);
        }
    }

    b /= Nseeds;

    cout << "b=" << b << endl;

    double b_result = 0;
    double msq = 0;
    //Initialize msq
    for (int y = 0; y < H; y++)
    {
        for (int x = 0; x < W; x++)
        {
            double res = ptracing_results[y][x].average();
            msq += res*res;
        }
    }
    msq /= (double)(H*W);


    path p0;
    p0.I = 0;
    path_copy(p0, p_init);

    memset(img, 0, W*H*sizeof(Vector3));

    for (long i = 1; i <= Nsamples; i++)
    {
        if (i % (Nsamples/1000) == 0)
        {
            printf("\rProgress: %f%% MSQ: %lf", (float)(100.*(double)i/(double)Nsamples), msq);
            fflush(stdout);
        }

        //Compute error vs. reference
        if (i % error_interval == 0)
        {
            msq = 0;
            bool writeImage = false;
            if (i == 100000 || i == 1000000 || i == 10000000 || i == 100000000)
                writeImage = true;

            for (int y = 0; y < H; y++)
            {
                for (int x = 0; x < W; x++)
                {
                    //Compute pixel value
                    Vector3 result = img[x][y]*b*(double)W*(double)H/(double)i;
                    msq += pow((ptracing_results[y][x] - result).average(), 2);

                    if (writeImage)
                    {
                        //Gamma correct
                        for (int i = 0; i < 3; i++)
                        {
                            result[i] = pow(result[i], 1.f/2.2f);
                        }
                        g_image->setPixel(x,y,result);
                    }
                }
            }

            msq /= (double)(W*H);
            msq_out << msq << endl;

            if (writeImage)
            {
                char filename[100];

                sprintf(filename, "metropolis_%s_%ld.ppm\0", version, i);
                cout << "\nWriting " << filename << "..." << endl;
                g_image->writePPM(filename);
            }
        }

        //Mutate path. 
        path p1;
        mutate_path(p0, p1);

        sample s = samplePath(p1, W, H);

        p1.I = s.value.average();
        p1.F = s.value;

        double accept;

        int x0 = (int)(p0.u[0]*(double)W), y0 = (int)(p0.u[1]*(double)H);
        int x1 = (int)(p1.u[0]*(double)W), y1 = (int)(p1.u[1]*(double)H);
        if (x0 == W) x0--;
        if (y0 == H) y0--;
        if (x1 == W) x1--;
        if (y1 == H) y1--;

        //Add contribution to pixels
        accept = std::min(p1.I / p0.I, 1.);
        if (p0.I > 0)
            img[x0][y0] += (1.-accept)*(p0.F / p0.I);
        if (p1.I > 0)
            img[x1][y1] += accept * (p1.F / p1.I);

        if (frand() < accept)
        {
            path_copy(p0, p1);
        }

    }

    printf("\n");
    for (int x = 0; x < W; x++)
    {
        for (int y = 0; y < H; y++)
        {
            //Compute pixel value
            Vector3 result = img[x][y]*b*(double)W*(double)H/(double)Nsamples;
            b_result += result.average();

            //Gamma correct
            for (int i = 0; i < 3; i++)
            {
                result[i] = pow(result[i], 1.f/2.2f);
            }
            g_image->setPixel(x,y,result);
        }
    }

    cout << "Resulting b: " << b_result/(double)W/(double)H << endl;
}
