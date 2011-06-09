#ifdef PHOTON_MAPPING

#include <fstream>
#include <sstream>
#include "Utility.h"
#include <cmath>
#include <iostream>
#include "Miro.h"
#include "PScene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include "Sphere.h"
#include "SquareLight.h"

#ifdef STATS
#include "Stats.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif 

//#define DEBUG_PHOTONS

#ifdef DEBUG_PHOTONS
    #ifdef OPENMP
        #define PHOTON_DEBUG(s) \
        _Pragma("omp critical") \
        cout << s << endl;
    #else
        #define PHOTON_DEBUG(s) cout << s << endl;
    #endif
#else
    #define PHOTON_DEBUG(s) 
#endif

using namespace std;

Scene * g_scene = 0;

void
Scene::openGL(Camera *cam)
{
#ifndef NO_GFX
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    cam->drawGL();

    // draw objects
    for (size_t i = 0; i < m_objects.size(); ++i)
        m_objects[i]->renderGL();

    glutSwapBuffers();
#endif
}

void
Scene::preCalc()
{
    debug("Precalcing objects\n");
    double t1 = -getTime();
    Objects::iterator it;
    for (it = m_objects.begin(); it != m_objects.end(); it++)
    {
        Object* pObject = *it;
        pObject->preCalc();
    }
    Lights::iterator lit;
    for (lit = m_lights.begin(); lit != m_lights.end(); lit++)
    {
        PointLight* pLight = *lit;
        pLight->preCalc();
    }
    t1 += getTime();
    debug("Time spent preprocessing objects and lights: %lf\n", t1);
    
    debug("Building BVH...\n");
    t1 = -getTime();
    m_bvh.build(&m_objects);
    t1 += getTime();
    debug("Done building BVH. Time spent: %lf\n", t1);

}

inline float tonemapValue(float value)
{
    return max(min(pow((double)value, 1./2.2), 1.), 0.);
}

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	int depth = TRACE_DEPTH;

    printf("Rendering Progress: %.3f%%\r", 0.0f);
    fflush(stdout);

    //For tone mapping. The Image class stores the pixels internally as 1 byte integers. We want to store the actual values first.
    int width = img->width(), height = img->height();
    Vector3 *tempImage = new Vector3[height*width];

    double t1 = -getTime();

    // loop over all pixels in the image
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Ray ray;
            Vector3 shadeResult(0.f);

            ray = cam->eyeRay(j, i, width, height, false);

			traceScene(ray, Vector3(1), depth);
			m_Points[m_Points.size()-1]->j = j;
			m_Points[m_Points.size()-1]->i = i;

			#ifdef STATS
			Stats::Primary_Rays++;
			#endif

        }
        printf("Rendering Progress: %.3f%%\r", i/float(img->height())*100.0f);
        fflush(stdout);
    }

	//if (m_Points.size() != (width*height))
	//	debug("uhohs\n");
    t1 += getTime();
    debug("Performing Adaptive passes...");
	m_pointMap.balance();

    t1 = -getTime();
	AdaptivePhotonPasses();
    t1 += getTime();

	RenderPhotonStats(tempImage, width, height);

    debug("Performing tone mapping...");

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Vector3 finalColor = tempImage[i*width+j];

            #pragma unroll(3)
            for (int k = 0; k < 3; k++)
            {
                if (finalColor[k] != finalColor[k])
                {
                    cout << "Pixel at " << j << "," << i << " is NAN!" << endl;
                }
                finalColor[k] = tonemapValue(finalColor[k]);
            }
            img->setPixel(j, i, finalColor);
        }
        #ifndef NO_GFX //If not rendering graphics to screen, don't draw scan lines (it will segfault in multithreading mode)
        img->drawScanline(i);
        #endif
    }
	m_pointMap.empty();

    printf("Rendering Progress: 100.000%%\n");
    debug("Done raytracing!\n");
    printf("Time spent raytracing image: %lf seconds.\n", t1);

#ifdef STATS
	Stats tracerStats;
	tracerStats.PrintStats();
#endif
}

bool
Scene::trace(HitInfo& minHit, const Ray& ray, float tMin, float tMax) const
{
    bool result = m_bvh.intersect(minHit, ray, tMin, tMax);
    
    //Trace the unbounded objects (like planes)
    for (int i = 0; i < m_unboundedObjects.size(); i++)
    {
        HitInfo tempMinHit;
        bool ubresult = m_unboundedObjects[i]->intersect(tempMinHit, ray, tMin, tMax);
        if (ubresult && (!result || tempMinHit.t < minHit.t))
        {
            result = true;
            minHit = tempMinHit;
            minHit.object = m_unboundedObjects[i];
        }
    }

	minHit.P -= ray.d*epsilon;

    if (result)
    {
        //Bump mapping
        float delta = 0.0001;
//        printf("object %p, material from %p\n", minHit.object, minHit.material);

        if (minHit.material->GetLookupCoordinates() == UV)
        {
            //Take a few samples to calculate the derivative
            tex_coord2d_t center = minHit.object->toUVCoordinates(minHit.P);
            float u = center.u, v = center.v;
            float u1 = minHit.material->bumpHeight2D(tex_coord2d_t(u-delta, v)), 
                  u2 = minHit.material->bumpHeight2D(tex_coord2d_t(u+delta, v)),
                  v1 = minHit.material->bumpHeight2D(tex_coord2d_t(u, v-delta)),
                  v2 = minHit.material->bumpHeight2D(tex_coord2d_t(u, v+delta));
                  
            //Approximate derivatives using central finite differences
            float dx = (u2-u1)/(2*delta),
                  dy = (v2-v1)/(2*delta);
            
            //Find two tangents
            float n[3] = { minHit.N.x, minHit.N.y, minHit.N.z };

            int m = 0;
            if (n[1] > n[0]) m = 1;
            if (n[2] > n[m]) m = 2;
            Vector3 randomVec(m == 2 ? -n[2] : 0, m == 0 ? -n[0] : 0, m == 1 ? -n[1] : 0);
           
            Vector3 t1 = cross(minHit.N, randomVec);
            minHit.N += dx*(cross(minHit.N, t1))-dy*(cross(minHit.N, cross(minHit.N, t1)));
            minHit.N.normalize();
        }   
        //Todo: implement for 3D
        //bumpHeight = minHit.material->bumpHeight3D(tex_coord3d_t(minHit.P.x, minHit.P.y, minHit.P.z));
    }
    return result;
}

bool Scene::traceScene(const Ray& ray, Vector3 contribution, int depth)
{
    HitInfo hitInfo;
    bool hit = false;
    
    if (depth >= 0)
    {
		if (trace(hitInfo, ray))
		{
            hit = true;

			--depth;

			//if diffuse material, send trace with RandomRay generate by Monte Carlo
			if (hitInfo.material->isDiffuse())
			{
				g_scene->addPoint(hitInfo.P, hitInfo.N, ray.d, contribution.average(), INITIAL_RADIUS, false);
			}
			
			//if reflective material, send trace with ReflectRay
			if (hitInfo.material->isReflective())
			{
				Ray reflectRay = ray.reflect(hitInfo);
				traceScene(reflectRay, contribution*hitInfo.material->getReflection(), depth);
			}

			//if refractive material, send trace with RefractRay
			if (hitInfo.material->isRefractive())
			{
			    float Rs = ray.getReflectionCoefficient(hitInfo); //Coefficient from fresnel

		        if (Rs > 0.01)
		        {
					//Send a reflective ray (Fresnel reflection)
					Ray reflectRay = ray.reflect(hitInfo);
		            traceScene(reflectRay, contribution * hitInfo.material->getRefraction() * Rs, depth);
		        }
			    
				Ray	refractRay = ray.refract(hitInfo);
				traceScene(refractRay, contribution * hitInfo.material->getRefraction() * (1.f-Rs), depth);
			}

            PointLight *l = dynamic_cast<PointLight*>(hitInfo.object);
			if (l != NULL)
			{
				//this means we hit an emissive material (light), so create a default measurement point
				g_scene->addPoint(Vector3(0.f), Vector3(0.f), Vector3(0.f), 0.f, 0.f, true);
                g_scene->m_Points.back()->accFlux = l->radiance(hitInfo.P, ray.d)*contribution.average();
			}
		}
	}
    
    return hit;
}

float Mutate(const float MutationSize)
{
	return ((2 * frand() - 1) > 0 ? 1 : -1) * pow(frand(), 1.f/MutationSize+1);
}

float ApplyDeltaRange(const float delta, float value, const float x1, const float x2)
{
	float range = x2 - x1;
    float result = value + delta;
    
	if (result < x1)
		result += range;
	else if ( result > x2)
		result -= range;

	return result;
}

Path MutatePath(const Path& goodPath, const float MutationSize)
{
	Path mutatedPath;

	float theta = atan2(goodPath.Direction.y, goodPath.Direction.x);
	float phi = acos(goodPath.Direction.z);

	//Align to squarelight direction
	mutatedPath.Direction = alignHemisphereToVector(Vector3(0,-1,0), ApplyDeltaRange(Mutate(MutationSize),theta, 0.f, 2.*PI),
																		ApplyDeltaRange(Mutate(MutationSize), phi, 0, PI/2.));

	// mutate phi/theta or random numbers?
	for (int i = 0; i < (TRACE_DEPTH_PHOTONS); ++i)
	{
		//theta = u1 (0 - 2PI), phi = u2 (0 - PI/2)
		mutatedPath.u[i*2] = ApplyDeltaRange(Mutate(MutationSize), goodPath.u[i*2], 0.f, 2.*PI);
		mutatedPath.u[i*2+1] = ApplyDeltaRange(Mutate(MutationSize), goodPath.u[i*2+1], 0, PI/2.);
	} 

	mutatedPath.Origin.x = ApplyDeltaRange(Mutate(MutationSize), goodPath.Origin.x, -10, 10);
	mutatedPath.Origin.y = goodPath.Origin.y;
	mutatedPath.Origin.z = ApplyDeltaRange(Mutate(MutationSize), goodPath.Origin.z, -1, 1);

	return mutatedPath;
}

bool Scene::UpdateMeasurementPoints(const Vector3& pos, const Vector3& normal, const Vector3& power)
{
	bool hit = false;
    const int npoints = W*H;

	NearestPoints np;

    np.dist2 = new float[npoints+1];
    np.index = new Point*[npoints+1];
	m_pointMap.find_points(&np, pos, normal, max_radius+epsilon, npoints);
	//m_pointMap.find_points(&np, pos, normal, INITIAL_RADIUS, npoints);

	for (int i=1; i<=np.found; i++) {
		Point *hp = np.index[i];

//	for (int n = 0; n <  m_Points.size(); ++n)
//	{
//		Point *hp = m_Points[n];

		// skip the measurement points that did not hit a surface
		if (hp->bLight)
        {
			continue;
        }

		if(dot(hp->normal, normal) < epsilon)
    		continue;

		float d = (pow(pos.x - hp->position.x, 2) +
			pow(pos.y - hp->position.y, 2) +
			pow(pos.z - hp->position.z, 2));

		if (d <= hp->radius*hp->radius)
		{
			//wait to update radius and flux * BRDF
			hp->newPhotons++;
//			hp->newFlux += power.x * hp->brdf;

            //The BRDF are taken into account with russian roulette.
			hp->newFlux += power.x;

			// can hit multiple measurement points	
			hit = true;
		}
	}
    delete[] np.dist2;
    delete[] np.index;

	return hit;
}

void Scene::UpdatePhotonStats()
{
	max_radius = 0.f;
	for (int n = 0; n < m_Points.size(); ++n)
	{
		Point *hp = m_Points[n];
		if (hp->bLight)
			continue;

		if (hp->radius > max_radius)
			max_radius = hp->radius;

		// continue if no new photons have been added
		if(hp->newPhotons == 0)
			continue;

        double f_alpha = (long double)hp->accPhotons*5e-6;

        float alpha = PHOTON_ALPHA;
//        float alpha = PHOTON_ALPHA + (1.-PHOTON_ALPHA)*(1.-exp(-f_alpha));

        // Set scaling factor for next photon pass
		hp->scaling = AdjustCorners(hp->radius, hp->position, hp->normal);;

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
}

void Scene::RenderPhotonStats(Vector3 *tempImage, const int width, const int height)
{
	// initialize for now
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
		{
			tempImage[i*width+j] = 0.f;
		}
	}

	int n;
	for (n = 0; n <  m_Points.size(); ++n)
	{
		Point *hp = m_Points[n];

		if (hp->bLight)
		{
			tempImage[hp->i*width+hp->j] = hp->accFlux;
			continue;
		}

		if (hp->accFlux < epsilon)
		{
			tempImage[hp->i*width+hp->j] = 0.f;
			continue;
		}

		long double A = PI * pow(hp->radius, 2);

		long double result = hp->accFlux / A / (long double)m_photonsEmitted * ((long double)m_photonsUniform / (long double)m_photonsEmitted)*hp->brdf;

		tempImage[hp->i*width+hp->j] = Vector3(result)/PI;
	}
	//if (n != (width*height))
	//	debug("Measurement points do not equal image dimensions");

	double sum = 0;
	for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
		{
			sum += tempImage[i*width+j].average();
		}
	}

	cout << "Average Radiance: " << sum/(double)(width*height) << endl;
}

void Scene::AdaptivePhotonPasses()
{
    Vector3 ptracing_results[W][H];
    Vector3 tempImage[W*H];
    stringstream msq_out;

    //Record the error after every this many samples
    const int error_interval = 100000;

#ifdef HACKER2
    const char* version = "sphere";
#elif defined (HACKER3)
    const char* version = "red";
#else
    const char* version = "gray";
#endif

    {
        //Load image from task 1
        ifstream ptracing;
        char filename[100];
        sprintf(filename, "pathtracing_%s.raw\0", version);
        ptracing.open(filename, ios::binary);

        int w_pt, h_pt;
        long double b_pt = 0;
        if (!ptracing.is_open())
        {
            w_pt = W; h_pt = H;
        }
        ptracing.read((char*)&w_pt, 4);
        ptracing.read((char*)&h_pt, 4);
        cout << "Loading path tracing results [width=" << w_pt << ", height " << h_pt << "]" << endl;

        for (int i = 0; i < h_pt; i++)
        {
            for (int j = 0; j < w_pt; j++)
            {
                Vector3 pix(0);
                if (ptracing.is_open())
                {
                    ptracing.read((char*)&pix.x, sizeof(float));
                    ptracing.read((char*)&pix.y, sizeof(float));
                    ptracing.read((char*)&pix.z, sizeof(float));
                }
                b_pt += pix.average();
                ptracing_results[j][i] = pix;
            }
        }
        b_pt /= (long double)(w_pt*h_pt);
        cout << "Done reading path tracing results. b = " << b_pt << endl;
    }

	PointLight *light = m_lights[0];
	
    Vector3 power = light->color() * light->wattage();

	Path goodPath;
	float prev_di = 1;
	long mutated = 1;
	long accepted = 0;

    //int Nphotons = 100000001;
    int Nphotons = 1000001;

	//find starting good path
	do
	{
        goodPath.Origin = light->samplePhotonOrigin();
        goodPath.Direction = light->samplePhotonDirection();
	} while (tracePhoton(goodPath, goodPath.Origin, goodPath.Direction, power, 0) == 0);

    long double msq = 0;
	for (m_photonsEmitted = 0; m_photonsEmitted < Nphotons; m_photonsEmitted++)
    {
		if (m_photonsEmitted > 0 && m_photonsEmitted % 10000 == 0)
			UpdatePhotonStats();

        if (m_photonsEmitted > 0 && m_photonsEmitted % 10000 == 0)
        {

			debug("Photons emitted %d of %d [%f%%] MSQ: %Lf      \r", m_photonsEmitted, Nphotons, 100.f*(float)m_photonsEmitted/(float)Nphotons, msq);
			//PrintPhotonStats();
        }

        //Compute error vs. reference
        if ((m_photonsEmitted+1) % error_interval == 0)
        {
            long i = m_photonsEmitted+1;
            printf("\n");
            RenderPhotonStats(tempImage, W, H);

            msq = 0;
            bool writeImage = false;
            if (i == 100000 || i == 1000000 || i % 10000000 == 0 || i == 100000000)
                writeImage = true;

            for (int y = 0; y < H; y++)
            {
                for (int x = 0; x < W; x++)
                {
                    Vector3 result = tempImage[x+y*W];
                    msq += pow((ptracing_results[x][y] - result).average(), 2);

                    if (writeImage)
                    {
                        //Gamma correct
                        for (int i = 0; i < 3; i++)
                        {
                            result[i] = pow(abs(result[i]), 1.f/2.2f);
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

                sprintf(filename, "adaptiveppm_%s_%ld.ppm\0", version, i);
                cout << "Writing " << filename << "..." << endl;
                g_image->writePPM(filename);
            }
        }
    



        //Test random photon path
        Path uniformPath(light->samplePhotonOrigin(), light->samplePhotonDirection());
		if (tracePhoton(uniformPath, uniformPath.Origin, uniformPath.Direction, power, 0) > 0)
		{
			goodPath = uniformPath;
			++m_photonsUniform;
			continue;
		}

		//Mutatation size
		long double di = prev_di + (1. / (long double)mutated) * ((long double)accepted/(long double)mutated - 0.234);

		// add mutation and convert back to cartesian coords
		Path mutatedPath = MutatePath(goodPath, di);

		++mutated;
		prev_di = di;

		// Test mutated photon path
		if (tracePhoton(mutatedPath, mutatedPath.Origin, mutatedPath.Direction, power, 0) > 0)
		{
			goodPath = mutatedPath;
			++accepted;
			continue;
		}
		// Reuse good path
		tracePhoton(goodPath, goodPath.Origin, goodPath.Direction, power, 0);
    }

	UpdatePhotonStats();

    //write msq to file
    {
        ofstream msq_outfile;
        char filename[100];
        sprintf(filename, "bidirectional_%s_msq.dat\0", version);
        msq_outfile.open(filename);
        msq_outfile << msq_out.str().c_str();
    }
}

//Trace a single photon through the scene
int Scene::tracePhoton(const Path& path, const Vector3& position, const Vector3& direction, const Vector3& power, int depth)
{
    if (depth >= TRACE_DEPTH_PHOTONS) return 0;
    PHOTON_DEBUG(endl << "tracePhoton(): pos " << position << ", dir " << direction << ", pwr " << power << ", depth " << depth);

    //Create a ray to trace the scene with
    Ray ray(position+epsilon*direction, direction);
    HitInfo hit;

	++depth;
    if (trace(hit, ray, 0.0f, MIRO_TMAX))
    {
		//SPECIAL: return if hit triangle backface
		if (dot(ray.d, hit.N) > 0)
			return 0;

        //Do russian roulette
        //Choose a random kind of ray - transmission, diffuse or reflective. Or absorb.
        //[ --diffuse-- | --specular (refl.)-- | --transmission-- | --absorb-- ]
        float prob[3], rnd = frand();
        Vector3 diffuseColor;
        if (hit.material->GetLookupCoordinates() == UV)
            diffuseColor = hit.material->diffuse2D(hit.object->toUVCoordinates(hit.P));
        else
            diffuseColor = hit.material->diffuse3D(tex_coord3d_t(hit.P.x, hit.P.y, hit.P.z));

        prob[0] = diffuseColor.average();
        prob[1] = prob[0] + hit.material->getReflection().average();
        prob[2] = prob[1] + hit.material->getRefraction().average();

        PHOTON_DEBUG("rnd = " << rnd << " F[reflect] = " << prob[0] << " F[refract] = " << prob[1] << "F[absorb] = " << prob[2]);

        if (rnd > prob[2])
        {
            //Absorb. Do nothing.
            PHOTON_DEBUG("Absorbed.");
            return 0;
        }

        if (rnd < prob[0])
        {
            PHOTON_DEBUG("Diffuse contribution.");
            int nPhotons = 0;

            //Diffuse, but store direct lighting for progressive mapping
			// only increment photons stored if hit a measurement point
			if (UpdateMeasurementPoints(hit.P, hit.N, power))
            {
				nPhotons++;
            }

#ifdef STATS
			Stats::Photon_Bounces++;
#endif
            //Shoot out a new diffuse photon 
			//depth-1 since already incremented
			Ray r = Ray::alignToVector(hit.N, hit.P, path.u[(depth-1)*2], path.u[(depth-1)*2+1]);
            r.isDiffuse = true;
            HitInfo diffHit;
            PHOTON_DEBUG("Tracing diffuse photon");
            return nPhotons + tracePhoton(path, r.o, r.d, diffuseColor*power/prob[0], depth);
        }
        else if (rnd < prob[1])
        {

#ifdef STATS
			Stats::Photon_Bounces++;
#endif
            //Reflect.
            Ray refl = ray.reflect(hit);
            PHOTON_DEBUG("Tracing reflected photon");
            return tracePhoton(path, hit.P, refl.d, power, depth);
        }
        else if (rnd < prob[2])
        {

#ifdef STATS
			Stats::Photon_Bounces++;
#endif
            //Transmit (refract)
            float Rs = ray.getReflectionCoefficient(hit); //Coefficient from fresnel
			
			//Fresnel reflection
			if (frand() < Rs)
			{
                Ray refl = ray.reflect(hit);
                PHOTON_DEBUG("Tracing reflected photon (Fresnel reflection)");
                return tracePhoton(path, hit.P, refl.d, power, depth);
			}
			else
			{
                Ray refr = ray.refract(hit);
                PHOTON_DEBUG("Tracing refracted photon");
                return tracePhoton(path, hit.P, refr.d, power, depth);
            }
        }
    }
#   ifdef DEBUG_PHOTONS
    else { PHOTON_DEBUG("Missed scene."); }
#   endif
	return 0;
}

Vector3
Scene::getEnvironmentMap(const Ray & ray)
{
	Vector3 envResult;
	//Environment mapping here
	if (m_environment != 0)
	{
		tex_coord2d_t coords;
		float phi = atan2(ray.d.x, ray.d.z) + m_environmentRotation.x + PI; //Phi is in [0, 2PI]
		float theta = asin(ray.d.y) + m_environmentRotation.y;    //
        if (theta > PI/2.0f) 
        {
            phi += PI;
            theta -= 2.0f*(theta-PI/2.0f);
        }
		if (phi > 2.0f*PI) phi -= (2.0f*PI); // Force phi to be in [0, 2PI]
	
		//Calculate texture coordinates for where the ray hits the "sphere"
		coords.u = phi / (2.0f * PI);
		coords.v = theta / PI + 0.5;
		//And just look up the shading value in the texture.
        if (!ray.isDiffuse)
    		envResult = m_environment->lookup2D(coords);
        else
            envResult = m_environment->lowresLookup2D(coords);
	}
	else
	{
		envResult = m_bgColor; 
	}
	return envResult;
}
#endif
