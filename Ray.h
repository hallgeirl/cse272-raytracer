#ifndef CSE168_RAY_H_INCLUDED
#define CSE168_RAY_H_INCLUDED

#include <stdlib.h>
#include <iostream>
#include "Vector3.h"
#include "Material.h"
#include "Utility.h"

#ifdef STATS
#include "Stats.h"
#endif

#include "SSE.h"
    
//! Contains information about a ray hit with a surface.
/*!
  HitInfos are used by object intersection routines. They are useful in
  order to return more than just the hit distance.
 */
class HitInfo
{
    public:
        float t;                            //!< The hit distance
        Vector3 P;                          //!< The hit point
        Vector3 N;                          //!< Shading normal vector
        const Material* material;           //!< Material of the intersected object
        Object  * object;             //!< Material of the intersected object

        //! Default constructor.
        explicit HitInfo(float t = 0.0f,
                const Vector3& P = Vector3(),
                const Vector3& N = Vector3(0.0f, 1.0f, 0.0f)) :
            t(t), P(P), N(N), material (0)
    {
        // empty
    }
};

class Ray
{
    public:
        bool isDiffuse;
        Vector3 o,      //!< Origin of ray
                d;      //!< Direction of ray


#ifdef __SSE4_1__
        __m128 d_SSE, o_SSE, d_SSE_rcp;

        void setupSSE()
        {
            static const __m128 _zero = _mm_setzero_ps();
            float _o[4] = {o.x, o.y, o.z, o.x};
            float _d[4] = {d.x, d.y, d.z, d.x};

            d_SSE = _mm_sub_ps(_zero, _mm_loadu_ps(_d));
            d_SSE_rcp = _mm_rcp_ps(_mm_loadu_ps(_d));
            o_SSE = _mm_loadu_ps(_o);
        }
#endif

    Ray() : o(), d(Vector3(0.0f,0.0f,1.0f))
    {
        isDiffuse = false;

#ifdef __SSE4_1__
        setupSSE();
#endif
    }

    Ray(const Vector3& o, const Vector3& d) : o(o), d(d)
    {
        isDiffuse = false;
#ifdef __SSE4_1__
        setupSSE();
#endif
    }

        //Returns a ray that is aligned to a vector v, given spherical coordinates theta and phi.
        static Ray alignToVector(const Vector3& v, const Vector3& origin, float theta, float phi) 
        {
            Vector3 dir = alignHemisphereToVector(v, theta, phi);

            return Ray(origin, dir);
        }

        //Diffuse ray with predetermined random numbers
        Ray diffuse(const HitInfo & hitInfo, float u1, float u2) const
        {
            //bias to the surface normal
            float phi = asin(sqrt(u1));
            float theta = 2.0f * PI * u2;

            Ray random = alignToVector(hitInfo.N, hitInfo.P, theta, phi);
            random.isDiffuse = true;
            random.o += random.d*epsilon;

            return random;
        }

        Ray diffuse(const HitInfo & hitInfo) const
        {
            //bias to the surface normal
            float phi = asin(sqrt(frand()));
            float theta = 2.0f * PI * frand();

            Ray random = alignToVector(hitInfo.N, hitInfo.P, theta, phi);
            random.isDiffuse = true;

            return random;
        }

        //Shoots a reflection ray. If path tracing is enabled, shoot a random ray according to the glossyness of the material.
        Ray reflect(const HitInfo & hitInfo) const
        {

//#ifdef PATH_TRACING
//            //Generate randomized reflection ray based on glossyness
//            //bias to the perfectly reflected ray
//            float phi = acos(pow(frand(), 1.0f/(1.0f+hitInfo.material->getShininess())));
//            float theta = 2.0f * PI * frand();
//
//            //Direction of perfect reflection
//            Vector3 d_reflect = d - 2 * dot(hitInfo.N, d) * hitInfo.N;
//
//            return alignToVector(d_reflect, hitInfo.P, theta, phi);
//#else
            Vector3 d_r = d - 2. * dot(hitInfo.N, d) * hitInfo.N;
			d_r.normalize();
            Ray reflect(hitInfo.P, d_r);
            return reflect;
//#endif
        }

        //Reflection coefficient from the Fresnel equations
        float getReflectionCoefficient(const HitInfo &hitInfo) const
        {
            float n1, n2;
            Vector3 n;
                       
            // if ray enters object, else ray exits object
            if ( dot(d, hitInfo.N) < 0)
            {
                n1 = 1.0f;
                n2 = hitInfo.material->getRefractionIndex();
                n = hitInfo.N;
            }
            else 
            {
                n1 = hitInfo.material->getRefractionIndex();
                n2 = 1.0f;
                n = -hitInfo.N;
            }
            
            float cosTheta = dot(-d, n);
            float sinTheta = sin(acos(cosTheta));
            float powsomething = std::pow((n1/n2) * sinTheta, 2.f);
            
            //If above critical angle, return 1.
            if (powsomething > 1.f) return 1;
            
            float sqrtSinTheta = sqrt(1.f - (powsomething));

            //The equation is (n1*cos(th) - n2 * sqrt(1-((n1/n2)*sin(th))^2)) / (n1*cos(th) + n2 * sqrt(1-((n1/n2)*sin(th))^2))
            float res = std::pow((n1*cosTheta - sqrtSinTheta)/(n1*cosTheta + sqrtSinTheta), 2.f);

            return res;
        }

        Ray refract(const HitInfo & hitInfo) const
        {
            float n1, n2;
            Vector3 n;

            // if ray enters object, else ray exits object
            if ( dot(d, hitInfo.N) < 0)
            {
                n1 = 1.0f;
                n2 = hitInfo.material->getRefractionIndex();
                n = hitInfo.N;
            }
            else 
            {
                n1 = hitInfo.material->getRefractionIndex();
                n2 = 1.0f;
                n = -hitInfo.N;
            }

            float energy = 1 - (pow(n1, 2) * (1 - pow(dot(d, n), 2)) / pow(n2, 2));

            // Total internal reflection: all of the energy is reflected
            if (energy < 0)
            {
//                std::cout << "int refl" << std::endl;
                return reflect(hitInfo);
            }

            Vector3 d_r = n1 * (d - n * dot(d, n)) / n2 - n * sqrt(energy);

//            #ifdef PATH_TRACING
//            float phi = acos(pow(frand(), 1.0f/(1.0f+hitInfo.material->getShininess())));
//            float theta = 2.0f * PI * frand();
//
//            return alignToVector(d_r, hitInfo.P, theta, phi);
//            #else
            return Ray(hitInfo.P, d_r);
//            #endif
        }

};


#endif // CSE168_RAY_H_INCLUDED
