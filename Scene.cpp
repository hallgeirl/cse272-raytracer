
#ifndef PHOTON_MAPPING
#include "Utility.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "Miro.h"
#include "Scene.h"
#include "Camera.h"
#include "Image.h"
#include "Console.h"
#include "Sphere.h"
#include "DirectionalAreaLight.h"

#ifdef STATS
#include "Stats.h"
#endif

#ifdef OPENMP
#include <omp.h>
#endif 

#define PATH_TRACING

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
//    value = min(max(value,0.f), 1.f);
    //gamma correction
    return min(max(pow(value, 1./2.2), 0.), 1.);;
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
        float localMaxIntensity = -infinity,
              localMinIntensity = infinity;

        for (int j = 0; j < width; ++j)
        {
            Ray ray;
            Vector3 tempShadeResult;
            Vector3 shadeResult(0.f);

			#if defined (PATH_TRACING) || defined(DOF)
			for (int k = 0; k < TRACE_SAMPLES; ++k)
			{
                ray = cam->eyeRay(j, i, width, height, false);
				if (traceScene(ray, tempShadeResult, depth))
				{
					shadeResult += tempShadeResult;
				}
			}
			shadeResult /= TRACE_SAMPLES; 
            tempImage[i*width+j] = shadeResult;
			#else
            ray = cam->eyeRay(j, i, width, height, false);
			if (traceScene(ray, shadeResult, depth))
                tempImage[i*width+j] = shadeResult;
            else
                tempImage[i*width+j] = m_bgColor;

			#endif // PATH_TRACING
            for (int k = 0; k < 3; k++)
            {
                if (shadeResult[k] > localMaxIntensity)
                    localMaxIntensity = shadeResult[k];
                if (shadeResult[k] < localMinIntensity)
                    localMinIntensity = shadeResult[k];
            }
            #ifdef OPENMP
            #pragma omp critical
			#endif
            {
                if (localMinIntensity < minIntensity) minIntensity = localMinIntensity;
                if (localMaxIntensity > maxIntensity) maxIntensity = localMaxIntensity;
            }

        }
        #ifdef OPENMP
        if (omp_get_thread_num() == 0)
        #endif
        {
            printf("Rendering Progress: %.3f%%\r", i/float(img->height())*100.0f);
            fflush(stdout);
        }
    }
    debug("Performing tone mapping...");
    t1 += getTime();

    ofstream img_raw("pathtracing.raw", ios::binary);

    img_raw.write((char*)&width, 4);
    img_raw.write((char*)&height, 4);

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            Vector3 finalColor = tempImage[i*width+j];

            //store raw data in rgb tuples
            img_raw.write((char*)&finalColor.x, sizeof(float));
            img_raw.write((char*)&finalColor.y, sizeof(float));
            img_raw.write((char*)&finalColor.z, sizeof(float));

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
    img_raw.close();

    printf("Rendering Progress: 100.000%%\n");
    debug("Done raytracing!\n");
    printf("Time spent raytracing image: %lf seconds.\n", t1);
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
		// AL: shouldn't decrementing depth be independent if there was a trace hit?
		if (trace(hitInfo, ray))
		{
            hit = true;

			--depth;
            const PointLight *l = dynamic_cast<PointLight*>(hitInfo.object);
            if (l)
            {
                shadeResult = Vector3(l->radiance(l->samplePhotonOrigin(), l->position() - hitInfo.P));
                return true;
            }

            double prob[3];
            prob[0] = hitInfo.material->getDiffuse().average();
            prob[1] = prob[0] + hitInfo.material->getReflection().average();
            prob[2] = prob[1] + hitInfo.material->getRefraction().average();

            double rnd = frand();
            if (rnd > prob[2])
            {
                return false;
            }

            //Diffuse reflection.
            if (rnd < prob[0])
			{
				Vector3 diffuseResult;
				Ray diffuseRay = ray.diffuse(hitInfo);

				if (traceScene(diffuseRay, diffuseResult, depth))
					shadeResult = diffuseResult * hitInfo.material->getDiffuse() / hitInfo.material->getDiffuse().average();
            }
            else if (rnd < prob[1])
            {
				Vector3 reflectResult;
				Ray reflectRay = ray.reflect(hitInfo);
				if (traceScene(reflectRay, reflectResult, depth))
					shadeResult = reflectResult * hitInfo.material->getReflection() / hitInfo.material->getReflection().average();
            }
            else if (rnd < prob[2])
            {
                //Push the hit point inside the refractive object (or outside if on the way out)
                hitInfo.P += ray.d*epsilon*2.;
			    float Rs = ray.getReflectionCoefficient(hitInfo); //Coefficient from fresnel

                if (frand() < Rs)
                {
					//Send a reflective ray (Fresnel reflection)
					Vector3 reflectResult;
					Ray reflectRay = ray.reflect(hitInfo);
		            if (traceScene(reflectRay, reflectResult, depth))
				        shadeResult = reflectResult * hitInfo.material->getRefraction() / hitInfo.material->getRefraction().average();
		        }
                else
                {
                    Vector3 refractResult;
                    Ray	refractRay = ray.refract(hitInfo);
                    if (traceScene(refractRay, refractResult, depth))
                    {
                        shadeResult = refractResult * hitInfo.material->getRefraction() / hitInfo.material->getRefraction().average();
                    }
                }
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
		}
	}
    
    return hit;
}

Vector3 Scene::getEnvironmentMap(const Ray & ray)
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
