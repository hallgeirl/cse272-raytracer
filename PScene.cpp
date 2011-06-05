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

inline float tonemapValue(float value, float maxIntensity)
{
    return value;
//    return sigmoid(20*value-2.5);
    //return std::min(pow(value / maxIntensity, 0.35f)*1.1f, 1.0f);
}

void
Scene::raytraceImage(Camera *cam, Image *img)
{
	int depth = TRACE_DEPTH;
    float minIntensity = infinity, maxIntensity = -infinity;

    printf("Rendering Progress: %.3f%%\r", 0.0f);
    fflush(stdout);

    //For tone mapping. The Image class stores the pixels internally as 1 byte integers. We want to store the actual values first.
    int width = img->width(), height = img->height();
    Vector3 *tempImage = new Vector3[height*width];

    double t1 = -getTime();

    // loop over all pixels in the image
    #ifdef OPENMP
    #pragma omp parallel for schedule(dynamic, 2)
    #endif
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Ray ray;
            Vector3 shadeResult(0.f);

            ray = cam->eyeRay(j, i, width, height, false);

			traceScene(ray, shadeResult, depth);

			#ifdef STATS
			Stats::Primary_Rays++;
			#endif

        }
        #ifdef OPENMP
        if (omp_get_thread_num() == 0)
        #endif
        {
            printf("Rendering Progress: %.3f%%\r", i/float(img->height())*100.0f);
            fflush(stdout);
        }
    }

	if (m_hitpoints.size() != (width*height))
		debug("uhohs\n");
    t1 += getTime();
    debug("Performing Adaptive passes...");

    t1 = -getTime();
	AdaptivePhotonPasses();
    t1 += getTime();

	RenderPhotonStats(tempImage, width, height, minIntensity, maxIntensity);

    debug("Performing tone mapping...");

    #ifdef OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Vector3 finalColor = tempImage[i*width+j];

            #pragma unroll(3)
            for (int k = 0; k < 3; k++)
            {
                if (finalColor[k] != finalColor[k])
                    finalColor[k] = maxIntensity;
                
                finalColor[k] = tonemapValue(finalColor[k], maxIntensity);
            }
            img->setPixel(j, i, finalColor);
        }
        #ifndef NO_GFX //If not rendering graphics to screen, don't draw scan lines (it will segfault in multithreading mode)
        img->drawScanline(i);
        #endif
    }

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

bool Scene::traceScene(const Ray& ray, Vector3& shadeResult, int depth)
{
    HitInfo hitInfo;
	shadeResult = Vector3(0.f);
    bool hit = false;
    
    if (depth >= 0)
    {
		if (trace(hitInfo, ray))
		{
            hit = true;
			bool bModified = false;

			--depth;

			//if diffuse material, send trace with RandomRay generate by Monte Carlo
			if (hitInfo.material->isDiffuse())
			{
				g_scene->addHitPoint(new HitPoint(hitInfo.P, hitInfo.N, ray.d, hitInfo.material->getDiffuse()[0]/PI, 0.25f));

				bModified = true;
/*#ifdef PHOTON_MAPPING

				float pos[3] = {hitInfo.P.x, hitInfo.P.y, hitInfo.P.z};
				float normal[3] = {hitInfo.N.x, hitInfo.N.y, hitInfo.N.z};
				float irradiance[3] = {0,0,0};
				float caustic[3] = {0,0,0};
            
				m_photonMap.irradiance_estimate(irradiance, pos, normal, PHOTON_MAX_DIST, PHOTON_SAMPLES);
				m_causticMap.irradiance_estimate(caustic, pos, normal, PHOTON_MAX_DIST, PHOTON_SAMPLES);

                //irradiance_estimate does the dividing by PI and all that
				shadeResult += Vector3(irradiance[0]+caustic[0], irradiance[1]+caustic[1], irradiance[2]+caustic[2]);
#endif*/
			}
			
			//if reflective material, send trace with ReflectRay
			if (hitInfo.material->isReflective())
			{
				Vector3 reflectResult;
				Ray reflectRay = ray.reflect(hitInfo);
				if (traceScene(reflectRay, reflectResult, depth))
				{
					shadeResult += hitInfo.material->getReflection() * reflectResult;
				}
				bModified = true;
			}

			//if refractive material, send trace with RefractRay
			if (hitInfo.material->isRefractive())
			{
			    float Rs = ray.getReflectionCoefficient(hitInfo); //Coefficient from fresnel

		        if (Rs > 0.01)
		        {
					//Send a reflective ray (Fresnel reflection)
					Vector3 reflectResult;
					Ray reflectRay = ray.reflect(hitInfo);
		            if (traceScene(reflectRay, reflectResult, depth))
			        {
				        shadeResult += hitInfo.material->getRefraction() * reflectResult * Rs;
			        }
		        }
			    
				Vector3 refractResult;
				Ray	refractRay = ray.refract(hitInfo);
				if (traceScene(refractRay, refractResult, depth))
				{
					shadeResult += hitInfo.material->getRefraction() * refractResult * (1.f-Rs);
				}

				bModified = true;
			}
			if (!bModified)
			{
				//this means we hit an emissive material, so create a default measurement point
				g_scene->addHitPoint(new HitPoint());
			}
		}
		else
		{
            if (m_environment != 0)
            {
                shadeResult = getEnvironmentMap(ray);
                hit = true;
            }
            else hit = false;
			g_scene->addHitPoint(new HitPoint());
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
    float delta_ = delta;
    /*if (delta < -range) delta_ = -range;
    if (delta > range) delta_ = range;*/
    float result = value + delta_;
    
	if (result < x1)
		result += range;
	else if ( result > x2)
		result -= range;
/*    if (delta > 1000 || delta < -1000)
        cout << delta << endl; */

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
	for (int i = 0; i < (TRACE_DEPTH_PHOTONS-1); ++i)
	{
		//theta = u1 (0 - 2PI), pi = u2 (0 - PI/2)
		mutatedPath.u[i*2] = ApplyDeltaRange(Mutate(MutationSize), goodPath.u[i*2], 0.f, 2.*PI);
		mutatedPath.u[i*2+1] = ApplyDeltaRange(Mutate(MutationSize), goodPath.u[i*2+1], 0, PI/2.);
	} 

	mutatedPath.Origin.x = ApplyDeltaRange(Mutate(MutationSize), goodPath.Origin.x, -10, 10);
	mutatedPath.Origin.z = ApplyDeltaRange(Mutate(MutationSize), goodPath.Origin.z, -1, 1);

	return mutatedPath;
}

bool Scene::UpdateMeasurementPoints(const Vector3& pos, const Vector3& power)
{
	bool hit = false;

	for (int n = 0; n <  m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];

		// skip the measurement points that did not hit a surface
		if (!hp->bHit)
			continue;

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

void Scene::UpdatePhotonStats()
{
    cout << "Photon Stats:" << endl;
	for (int n = 0; n < m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];
		if (!hp->bHit)
			continue;

        double f_alpha = (long double)hp->accPhotons*5e-6;

        float alpha = PHOTON_ALPHA + (1.-PHOTON_ALPHA)*(1.-exp(-f_alpha));

        // Set scaling factor for next photon pass
        float A = PI * pow(hp->radius, 2);

        //CircleSegment(Vector3(-1,0,-1), Vector3(0,0,1), hp->radius, hp->position, A1); 
        //CircleSegment(Vector3(1,0,-1), Vector3(0,0,1), hp->radius, hp->position, A1);
        //m_hitpoints[n]->scaling = A1/A;
		
		// only adding a ratio of the newly added photons
		float delta = (hp->accPhotons + alpha * hp->newPhotons)/(hp->accPhotons + hp->newPhotons);
		hp->radius *= sqrt(delta);
		hp->accPhotons += (int)(alpha * hp->newPhotons);
		
		// not sure about this flux acc, or about calculating the irradiance
		hp->accFlux = ( hp->accFlux + hp->newFlux/hp->scaling) * delta;	

        cout << hp->accFlux << "\t" << hp->accPhotons << "\t" << hp->accPhotons << "\n";

		// reset new values
		hp->newPhotons = 0;
		hp->newFlux = 0.f;
	}
    cout << endl;
}

//void PrintPhotonStats(ofstream& fp, const float photonsEmitted, const float m_photonsUniform)
void Scene::PrintPhotonStats()
{
	for (int n = 0; n <  m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];

		long double A = PI * pow(hp->radius, 2);

		long double result = hp->accFlux / A / (long double)m_photonsEmitted * ((long double)m_photonsUniform / (long double)m_photonsEmitted);
	}
}

void Scene::RenderPhotonStats(Vector3 *tempImage, const int width, const int height, float minIntensity, float maxIntensity)
{
	float localMaxIntensity = -infinity, localMinIntensity = infinity;

	int n;
	for (n = 0; n <  m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];

		if (!hp->bHit)
		{
			tempImage[n] = m_bgColor;
			continue;
		}

		//tempImage[n] = Vector3(1.f);

		long double A = PI * pow(hp->radius, 2);

		long double result = hp->accFlux / A / (long double)m_photonsEmitted * ((long double)m_photonsUniform / (long double)m_photonsEmitted);

		tempImage[n] = Vector3(result);

		if (tempImage[n].x < minIntensity) minIntensity = tempImage[n].x;
        if (tempImage[n].x > maxIntensity) maxIntensity = tempImage[n].x;
	}
	if (n != (width*height))
		debug("Measurement points do not equal image dimensions");

	/*for (int k = 0; k < 3; k++)
    {
        if (shadeResult[k] > localMaxIntensity)
            localMaxIntensity = shadeResult[k];
        if (shadeResult[k] < localMinIntensity)
            localMinIntensity = shadeResult[k];
    }*/
}


bool Scene::SamplePhotonPath(const Path& path, const Vector3& power)
{
	return (tracePhoton(path, path.Origin, path.Direction, power, 0, false) > 0);

/*    HitInfo hitInfo(0, path.o, Vector3(0,1,0));

	Ray ray(path.o, path.d);

    while (true)
    {
        if (g_scene->trace(hitInfo, ray, 0, MIRO_TMAX))
        {
			//return if hit triangle backface
			if (dot(ray.d, hitInfo.N) > 0)
				return false;

            //hit diffuse surface->we're done
            if (hitInfo.material->isDiffuse())
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
    }*/
}

void Scene::AdaptivePhotonPasses()
{
	PointLight *light = m_lights[0];
	
    Vector3 power = light->color() * light->wattage();

	Path goodPath;
	float prev_di = 1;
	long mutated = 1;
	long accepted = 0;

    int Nphotons = 100000;

	//find starting good path
	do
	{
        goodPath.Origin = light->samplePhotonOrigin();
        goodPath.Direction = light->samplePhotonDirection();
	} while (!SamplePhotonPath(goodPath, power));

    long double msq = 0;
	for (m_photonsEmitted = 0; m_photonsEmitted < Nphotons; m_photonsEmitted++)
    {
		printf("photons emitted: %d\n", m_photonsEmitted);
		if (m_photonsEmitted > 0 && m_photonsEmitted % 1000000 == 0)
		{
			UpdatePhotonStats();
		}
        /*if (m_photonsEmitted > 0 && m_photonsEmitted % 100000 == 0)
        {
			PrintPhotonStats();
        }*/

        //Test random photon path
        Path uniformPath(light->samplePhotonOrigin(), light->samplePhotonDirection());
		if (SamplePhotonPath(uniformPath, power))
		{
			goodPath = uniformPath;
			++m_photonsUniform;
			continue;
		}

		//Mutatation size
		long double di = prev_di + (1. / (long double)mutated) * ((long double)accepted/(long double)mutated - 0.234);

		// Convert to spherical coords (theta phi reversed)
		/*float phi = acos(goodPath.d.z);
		float theta = atan2(goodPath.d.y, goodPath.d.x);*/

		// add mutation and convert back to cartesian coords
		Path mutatedPath = MutatePath(goodPath, di);

		++mutated;
		prev_di = di;

		// Test mutated photon path
		if (SamplePhotonPath(mutatedPath, power))
		{
			goodPath = mutatedPath;
			++accepted;
			continue;
		}
		// Reuse good path
		SamplePhotonPath(goodPath, power);		
    }

	UpdatePhotonStats();
	PrintPhotonStats();
}

void Scene::ProgressivePhotonPass()
{
	traceProgressivePhotons();
	
	//iterate through all of the scene hitpoints
	for (int n = 0; n < m_hitpoints.size(); ++n)
	{
		HitPoint *hp = m_hitpoints[n];

		float pos[3] = {hp->position.x, hp->position.y, hp->position.z};
		float normal[3] = {hp->normal.x, hp->normal.y, hp->normal.z};
		float irradiance[3] = {0,0,0};
    
		int M = m_photonMap.irradiance_estimate(irradiance, pos, normal, hp->radius, PHOTON_SAMPLES, false);
		
		//only adding a ratio of the newly added photons
		float delta = (hp->accPhotons + PHOTON_ALPHA * M)/(hp->accPhotons + M);
		hp->radius *= sqrt(delta);
		hp->accPhotons += (int)(PHOTON_ALPHA * M);
		
		//not sure about this flux acc, or about calculating the irradiance
		hp->accFlux = ( hp->accFlux + irradiance[0]/hp->scaling ) * delta;	
	}

	m_photonMap.empty();
} 

//Shoot out all photons and trace them
void Scene::traceProgressivePhotons()
{
    if (PhotonsPerLightSource == 0) 
    {
        m_photonMap.balance();
        return;
    }
    
    int photonsAdded = 0; //Photons added to the scene
    
    for (int l = 0; l < m_lights.size(); l++)
    {
        PointLight *light = m_lights[l];
    
        #ifdef OPENMP
        #pragma omp parallel for
        #endif
        for (int i = 0; i < PhotonsPerLightSource; i++)
        {
            //Create a new photon
            Vector3 power = light->color() * light->wattage(); 
            Vector3 dir = light->samplePhotonDirection();
            Vector3 pos = light->samplePhotonOrigin();
			//printf("squarelight photon position: %f %f %f and direction %f %f %f \n", pos.x, pos.y, pos.z, dir.x, dir.y, dir.z);
            //tracePhoton(pos, dir, power, 0);
        }

		m_photonsEmitted += PhotonsPerLightSource;
    }
	// do not scale photons in progressive photon mapping
    // m_photonMap.scale_photon_power(1.0f/(float)totalPhotons);
    m_photonMap.balance();
    #ifdef VISUALIZE_PHOTON_MAP
    debug("Rebuilding BVH for visualization. Number of objects: %d\n", m_objects.size());
    m_bvh.build(&m_objects);

    #endif
}

//Trace a single photon through the scene
int Scene::tracePhoton(const Path& path, const Vector3& position, const Vector3& direction, const Vector3& power, int depth, bool bCausticRay)
{
    PHOTON_DEBUG(endl << "tracePhoton(): pos " << position << ", dir " << direction << ", pwr " << power << ", depth " << depth);
    if (depth > TRACE_DEPTH_PHOTONS) return 0;

    //Create a ray to trace the scene with
    Ray ray(position+epsilon*direction, direction);
    HitInfo hit;

	++depth;
    if (trace(hit, ray, 0.0f, MIRO_TMAX))
    {
		//SPECIAL: return if hit triangle backface
		if (dot(ray.d, hit.N) > 0)
			return 0;

        //Do "russian roulette but not really"
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
            //Diffuse.
            //only store indirect lighting -- but store direct lighting for progressive mapping
            //if (depth > 1)
            {
				UpdateMeasurementPoints(hit.P, power);
				nPhotons++;
                //float pos[3] = {hit.P.x, hit.P.y, hit.P.z}, dir[3] = {direction.x, direction.y, direction.z}, pwr[3] = {power.x, power.y, power.z};
#               //ifdef OPENMP
#               //pragma omp critical
#              // endif
                {

                    /*PHOTON_DEBUG("Storing photon at " << hit.P << ". Surface normal " << hit.N << " Is caustic photon: " << bCausticRay);
                    if (bCausticRay)
						m_causticMap.store(pwr, pos, dir);
					else
						m_photonMap.store(pwr, pos, dir);

					nPhotons++;*/

#                   ifdef VISUALIZE_PHOTON_MAP
                    Sphere* sp = new Sphere;
                    sp->setCenter(hit.P); sp->setRadius(0.02f);
                    Vector3 ref = power; //Use the normalized power as the reflectance for visualization.
                    ref.normalize(); sp->setMaterial(new Phong(ref)); addObject(sp);
#                   endif
                }
            }
            /*else
            {
			    //Caustic Rays only send rays from specular surfaces
			    if (bCausticRay)
				    return 0;
		    }*/

#ifdef STATS
			Stats::Photon_Bounces++;
#endif
            //Shoot out a new diffuse photon 
			//depth-1 since already incremented
			Ray r = Ray::alignToVector(hit.N, hit.P, path.u[(depth-1)*2], path.u[(depth-1)*2+1]);
            r.isDiffuse = true;
            HitInfo diffHit;
            PHOTON_DEBUG("Tracing diffuse photon");
            return nPhotons + tracePhoton(path, r.o, r.d, diffuseColor*power/prob[0], depth, bCausticRay);
        }
        else if (rnd < prob[1])
        {
			//only caustics should count this first bounce
			//if (!bCausticRay && depth == 1)
				//return 0;

#ifdef STATS
			Stats::Photon_Bounces++;
#endif
            //Reflect.
            Ray refl = ray.reflect(hit);
            PHOTON_DEBUG("Tracing reflected photon");
            return tracePhoton(path, hit.P, refl.d, power, depth, bCausticRay);
        }
        else if (rnd < prob[2])
        {
			//only caustics should count this first bounce
			//if (!bCausticRay && depth == 1)
			//	return 0;

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
                return tracePhoton(path, hit.P, refl.d, power, depth, bCausticRay);
			}
			else
			{
                Ray refr = ray.refract(hit);
                PHOTON_DEBUG("Tracing refracted photon");
                return tracePhoton(path, hit.P, refr.d, power, depth, bCausticRay);
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
