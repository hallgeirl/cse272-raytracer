#ifndef _SQUARELIGHT_H_
#define _SQUARELIGHT_H_
#include "PointLight.h"
#include "Emissive.h"
#include <iostream>
#include "TriangleMesh.h"
#include "Triangle.h"
#include <cmath>
#include <stdio.h>

// Work around for windows error 
using std::min;
using std::max;

class SquareLight : public PointLight
{
public:
    SquareLight() { hasTangent1 = false; m_normal = Vector3(0,1,0); m_material = new Emissive; }

    virtual void setWattage(float f) 
    {
        m_wattage = f; //Divide by PI for diffuse light sources
        ((Emissive*)m_material)->setPower(f/(PI*m_dimensions[0]*m_dimensions[1])); 
    }

    virtual Vector3 coordsMin() const { return m_minCorner; }

    virtual Vector3 coordsMax() const { return m_maxCorner; }

    virtual Vector3 center() const    { return m_position;  }

    Vector3 getTangentU() const { return m_tangent1; } 
    Vector3 getTangentV() const { return m_tangent1; } 

    virtual float radiance(const Vector3& x, const Vector3& incdir) const 
    { 
        return m_wattage / (PI*m_dimensions[0]*m_dimensions[1]); 
    } 

    void setNormal(Vector3 n) 
    { 
        m_normal = n;
    }
    
    Vector3 getNormal()
    {
        return m_normal;
    }

    virtual bool intersect(HitInfo& hit, const Ray& ray, float tmin, float tmax)
    {
        if (dot(ray.d, m_normal) > 0) return false;
        HitInfo tmphit;
        bool hit1 = m_triangles[0].intersect(hit, ray, tmin, tmax);
     
        if (m_triangles[1].intersect(tmphit, ray, tmin, tmax))
        {
            if ((hit1 && tmphit.t < hit.t) || !hit1)
                hit = tmphit;

            hit1 = true;
        }

        if (hit1) hit.material = m_material;

        return hit1;
    }

    void setDimensions(float width, float height) { m_dimensions[0] = width; m_dimensions[1] = height; }
    float area() { return m_dimensions[0] * m_dimensions[1]; }
    void setUdir(const Vector3& udir) { m_tangent1 = udir; hasTangent1 = true; }

    //Get the position on the square light given two coordinates in [0,1]
    virtual Vector3 getPhotonOrigin(double u1, double u2) const  
    {
        float du = m_dimensions[0];
        float dv = m_dimensions[1];

		// Add epsilon in normal direction to prevent photons from hitting light mesh
        return m_position + m_tangent1*((u1-0.5)*du) + m_tangent2*((u2-0.5)*dv) + m_normal*epsilon;
    }
    virtual Vector3 samplePhotonOrigin(int sampleNumber = 0, int totalSamples = 1) const  
    {
        //Take samples within a subdivided rectangle. For simplicity we assume that the light is square so we have nxn cells.
        //First find the cell dimensions
        float sideLength = sqrt((float)totalSamples);
        float du = m_dimensions[0] / sideLength;
        float dv = m_dimensions[1] / sideLength;

        //Sample index in grid
        int sx = sampleNumber % int(sideLength);
        int sy = sampleNumber / int(sideLength);

        float u = (du*frand()) + sx * du - m_dimensions[0]/2.0f;
        float v = (dv*frand()) + sy * dv - m_dimensions[1]/2.0f;

		// Add epsilon in normal direction to prevent photons from hitting light mesh
        return m_position + u*m_tangent1 + v*m_tangent2 + m_normal * epsilon;
    }

    virtual Vector3 samplePhotonDirection() const
    {
        //bias to the light normal
        float phi = asin(sqrt(frand()));
        float theta = 2.0f * PI * (frand());

        return alignHemisphereToVector(m_normal, theta, phi);
    }

    virtual Vector3 samplePhotonDirection(double u1, double u2) const
    {
        //bias to the light normal
        float phi = asin(sqrt(u1));
        float theta = 2.0f * PI * u2;

        return alignHemisphereToVector(m_normal, theta, phi);
    }

    virtual void preCalc()
    {
        if (!hasTangent1)
        {
            std::cerr << "Warning (SquareLight): No u-tangent set; choosing arbitrary tangents." << std::endl;
            getTangents(m_normal, m_tangent1, m_tangent2);
        }
        else
        {
            m_tangent1.normalize();
            m_tangent2 = cross(m_tangent1, m_normal);
            if (std::abs(dot(m_tangent1, m_normal)) > 1e-5)
            {
                std::cerr << "Warning (SquareLight): Chosen u-tangent not perpendicular to normal. Correcting." << std::endl;
                m_tangent1 = cross(m_tangent1, m_normal);
            }
        }

        float du = m_dimensions[0] / 2.;
        float dv = m_dimensions[1] / 2.;

        Vector3 v1 = m_position + du*m_tangent1 + dv*m_tangent2;
        Vector3 v2 = m_position - du*m_tangent1 - dv*m_tangent2;

		m_minCorner.set(min(v1.x, v2.x), min(v1.y, v2.y), min(v1.z, v2.z));
        m_maxCorner.set(max(v1.x, v2.x), max(v1.y, v2.y), max(v1.z, v2.z));

        m_mesh[0].createSingleTriangle();
        m_mesh[0].setV1(m_position + du*m_tangent1 + dv*m_tangent2);
        m_mesh[0].setV2(m_position + du*m_tangent1 - dv*m_tangent2);
        m_mesh[0].setV3(m_position - du*m_tangent1 + dv*m_tangent2);
        m_mesh[0].setN(m_normal);

        m_mesh[1].createSingleTriangle();
        m_mesh[1].setV1(m_position - du*m_tangent1 - dv*m_tangent2);
        m_mesh[1].setV2(m_position - du*m_tangent1 + dv*m_tangent2);
        m_mesh[1].setV3(m_position + du*m_tangent1 - dv*m_tangent2);
        m_mesh[1].setN(m_normal);

        m_triangles[0].setIndex(0);
        m_triangles[0].setMesh(&m_mesh[0]);
        m_triangles[0].setMaterial(m_material);
        m_triangles[1].setIndex(0);
        m_triangles[1].setMesh(&m_mesh[1]);
        m_triangles[1].setMaterial(m_material);

    }
    
protected:
    TriangleMesh m_mesh[2];
    Triangle m_triangles[2];
    Vector3 m_normal; Vector3 m_tangent1, m_tangent2;
    float m_dimensions[2];
    bool hasTangent1;
    Vector3 m_minCorner, m_maxCorner;
};


#endif
