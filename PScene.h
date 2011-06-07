#ifdef PHOTON_MAPPING

#ifndef CSE168_PSCENE_H_INCLUDED
#define CSE168_PSCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "BVH.h"
#include "Texture.h"
#include "PhotonMap.h"
#include "PointMap.h"

class Camera;
class Image;

struct Path
{
	Vector3 Origin;
	Vector3 Direction;
	double u[TRACE_DEPTH_PHOTONS*2]; 

	Path(const Vector3& inOrigin, const Vector3& inDirection)
		:Origin(inOrigin), Direction(inDirection)
	{
		init_random();
	}

	Path()
	{
		init_random();
	}

	void init_random()
	{
		for (int i = 0; i < TRACE_DEPTH_PHOTONS*2; ++i)
		{
			u[i] = frand();
		}
	}
};

/*struct Point
{
	Vector3 position;
	Vector3 normal;
	Vector3 dir;
	float brdf;
	//float pixel_wgt;
	double radius;
	int accPhotons;
	int newPhotons;
	long double accFlux;
	long double newFlux;
    float scaling;
	bool bHit;

	Point()
		:position(0.f), normal(0.f), dir(0.f), brdf(1.f), radius(0.f), accPhotons(0), newPhotons(0), accFlux(0.f), newFlux(0.f), scaling(1), bHit(false)
	{}

		Point(const Vector3& inPosition, const Vector3& inNormal, const Vector3& inDir, const float inBRDF, const float inRadius)
		:position(inPosition), normal(inNormal), dir(inDir), brdf(inBRDF), radius(inRadius), accPhotons(0), newPhotons(0), accFlux(0.f), newFlux(0.f), scaling(1), bHit(true)
	{}
};*/

typedef std::vector<Point*> Points;

class Scene
{
public:
	Scene() 
		: m_photonMap(PhotonsPerLightSource*TRACE_DEPTH*MaxLights+MaxLights*10000), m_causticMap(CausticPhotonsPerLightSource*TRACE_DEPTH*MaxLights+MaxLights*10000), m_pointMap(W*H), m_environment(0), m_bgColor(Vector3(0.0f)), m_photonsEmitted(0), m_photonsUniform(0)
	{}
	// TODO: need right image dimensions
    void addObject(Object* pObj)        
    { 
        if (pObj->isBounded()) m_objects.push_back(pObj);
        else m_unboundedObjects.push_back(pObj);
		if (pObj->getMaterial()->isRefractive() || pObj->getMaterial()->isReflective()) m_specObjects.push_back(pObj);
    }
    const Objects* objects() const      {return &m_objects;}
    const Objects* specObjects() const      {return &m_specObjects;}

    void addLight(PointLight* pObj)     {m_lights.push_back(pObj);}
    const Lights* lights() const        {return &m_lights;}

	void addPoint(const Vector3& inPosition, const Vector3& inNormal, const Vector3& inDir, const float inBRDF, const float inRadius, const bool inbHit)	
	{
		if (inbHit)
		{
			Point* hp = m_pointMap.store(inPosition, inNormal, inDir, inRadius, inBRDF, inbHit);
			if (hp != NULL)
				m_Points.push_back(hp);
			else 
				m_Points.push_back(new Point());

//			m_Points.push_back(new Point(inPosition, inNormal, inDir, inBRDF, inRadius));
		}
		else
			m_Points.push_back(new Point());

	}
	//const Points* Points() const	{return &m_Points;}

    void generatePhotonMap();
	void AdaptivePhotonPasses();
    void ProgressivePhotonPass();

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    bool trace(HitInfo& minHit, const Ray& ray,
               float tMin = 0.0f, float tMax = MIRO_TMAX) const;
	bool traceScene(const Ray& ray, Vector3& shadeResult, int depth);

	void UpdatePhotonStats();
	void PrintPhotonStats();
	void RenderPhotonStats(Vector3 *tempImage, const int width, const int height, float minIntensity, float maxIntensity);
	bool SamplePhotonPath(const Path& path, const Vector3& power);
	bool UpdateMeasurementPoints(const Vector3& pos, const Vector3& normal, const Vector3& power);
    void traceProgressivePhotons();
    int tracePhoton(const Path& path, const Vector3& position, const Vector3& direction, const Vector3& power, int depth, bool bCausticRay=false);
	long int GetPhotonsEmitted() { return m_photonsEmitted; }

	void setEnvironment(Texture* environment) { m_environment = environment; }
	Vector3 getEnvironmentMap(const Ray & ray);

    void setBgColor(Vector3 color) { m_bgColor = color; }

    void setEnvironmentRotation(float phi, float theta) { m_environmentRotation.x = phi; m_environmentRotation.y = theta; }

protected:
    VectorR2 m_environmentRotation;
    Objects m_objects;
    Objects m_specObjects;
    Objects m_unboundedObjects;
    Photon_map m_photonMap;
    Photon_map m_causticMap;
	Point_map m_pointMap;
    BVH m_bvh;
    Lights m_lights;
	Points m_Points;
    Texture * m_environment; //Environment map
    Vector3 m_bgColor;       //Background color (for when environment map is not available)

    static const int MaxLights = 10;

    static const int PhotonsPerLightSource = 1;
    static const int CausticPhotonsPerLightSource = 1;

	long int m_photonsEmitted;
	long int m_photonsUniform;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
#endif
