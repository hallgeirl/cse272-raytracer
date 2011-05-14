#ifndef __UTILITY__H_
#define __UTILITY__H_

#include <cstdlib>
#include <cmath>
#include "Material.h"
#include "Vector3.h"
#include "Matrix4x4.h"

double getTime();
void getEigenVector(const float (&A)[3][3], float (&outV)[3], float lambda);
void addMeshTrianglesToScene(TriangleMesh * mesh, Material * material);
void addModel(const char* filename, Material *mat, Scene* scene, Vector3 position, float rotY=0, Vector3 scale=Vector3(1,1,1));

//Returns random number between 0 and 1
inline float frand()
{
#ifdef LINUX
    return drand48();
#else
    return (float)rand() / (float)RAND_MAX;
#endif
}

inline float sigmoid(float x)
{
    return 1.f/(1.f+exp(-x));
}

//Get two tangents to the normal vector that is supplied. Store the results in t1 and t2.
inline void getTangents(const Vector3& normal, Vector3& t1, Vector3& t2)
{
    //determine one surface tangent from cross of v and one of the unit vectors 
    t1 = cross(Vector3(0,0,1), normal);
    if (t1.length2() < 1e-6) t1 = cross(Vector3(0, 1, 0), normal);
    t2 = cross(t1, normal);
}

//Aligns a hemisphere to the vector v, and returns the vector formed by the spherical coordinates theta and phi.
inline Vector3 alignHemisphereToVector(const Vector3& v, float theta, float phi) 
{
    //convert spherical coords to cartesian coords
    float u1 = sin(phi) * cos(theta);
    float u2 = sin(phi) * sin(theta);
    float u3 = cos(phi);

    //determine one surface tangent from cross of v and one of the unit vectors 
    Vector3 t1 = cross(Vector3(0,0,1), v);
    if (t1.length2() < 1e-6) t1 = cross(Vector3(0, 1, 0), v);

    //Aligned direction determined by combination of cartesian coordiantes along tangents
    Vector3 aligned_d(u1*t1 + u2*cross(t1, v) + u3*v);

    aligned_d.normalize();
    return aligned_d;
}

//Returns a random direction in the hemisphere that is oriented in the direction specified
inline Vector3 sampleHemisphereDirection(const Vector3& hemisphereOrientation)
{
    //bias to the surface normal
    float x, y, z;
    do
    { 
        x = 2*frand() - 1;
        y = 2*frand() - 1;
        z = 2*frand() - 1;
    } while (x*x + y*y + z*z > 1.0f && dot(Vector3(x, y, z), hemisphereOrientation) < 0);

    return Vector3(x,y,z).normalize();
}

//Returns a random direction in the hemisphere that is oriented in the direction specified
inline Vector3 sampleSphericalDirection()
{
    //bias to the surface normal
    float x, y, z;
    do
    { 
        x = 2*frand() - 1;
        y = 2*frand() - 1;
        z = 2*frand() - 1;
    } while (x*x + y*y + z*z > 1.0f);

    return Vector3(x,y,z).normalize();
}

inline VectorR2 sampleDisc(float radius)
{
	float x_rand, y_rand;
	Vector3 new_eye;
	do {
		x_rand = (2*frand() - 1) * radius;
		y_rand = (2*frand() - 1) * radius;
	} while (x_rand*x_rand + y_rand*y_rand > radius*radius);

    VectorR2 v;
    v.x = x_rand; v.y = y_rand;

    return v;
}

inline Matrix4x4
translate(float x, float y, float z)
{
    Matrix4x4 m;
    m.setColumn4(Vector4(x, y, z, 1));
    return m;
}


inline Matrix4x4
scale(float x, float y, float z)
{
    Matrix4x4 m;
    m.m11 = x;
    m.m22 = y;
    m.m33 = z;
    return m;
}

// angle is in degrees
inline Matrix4x4
rotate(float angle, float x, float y, float z)
{
    float rad = angle*(PI/180.);
    
    float x2 = x*x;
    float y2 = y*y;
    float z2 = z*z;
    float c = cos(rad);
    float cinv = 1-c;
    float s = sin(rad);
    float xy = x*y;
    float xz = x*z;
    float yz = y*z;
    float xs = x*s;
    float ys = y*s;
    float zs = z*s;
    float xzcinv = xz*cinv;
    float xycinv = xy*cinv;
    float yzcinv = yz*cinv;
    
    Matrix4x4 m;
    m.set(x2 + c*(1-x2), xy*cinv+zs, xzcinv - ys, 0,
          xycinv - zs, y2 + c*(1-y2), yzcinv + xs, 0,
          xzcinv + ys, yzcinv - xs, z2 + c*(1-z2), 0,
          0, 0, 0, 1);
    return m;
}

inline float
CircleSegment(const Vector3& rayOrigin, const Vector3& rayDir, const float radius, const Vector3 center)
{
    const Vector3 toO = rayOrigin - center;

    const float a = rayDir.length2();
    const float b = dot(2*rayDir, toO);
    const float c = toO.length2() - radius*radius;

    const float discrim = b*b-4.0f*a*c;

    if (discrim < 0)
        return 1.f;   // quadratic equation would yield imaginary numbers

    const float sqrt_discrim = sqrt(discrim);

    // solve the quadratic equation
    const float t[2] = {(-b-sqrt_discrim)/(2.0f*a), (-b+sqrt_discrim)/(2.0f*a)};

    // since we know that discrim >= 0, t[0] < t{1]
    // return the t closest to us that is within range

	Vector3 intersectP1 = rayOrigin + rayDir * t[0];
	Vector3 intersectP2 = rayOrigin + rayDir * t[1];

	float theta = acos(dot((intersectP1 - center).normalize(), (intersectP2 - center).normalize()));
	float segArea = 0.5f * (theta - sin(theta)) * pow(radius, 2);

	return 1.f - segArea / (PI * pow(radius, 2));
}

#endif
