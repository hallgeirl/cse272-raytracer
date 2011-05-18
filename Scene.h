#ifndef CSE168_SCENE_H_INCLUDED
#define CSE168_SCENE_H_INCLUDED

#include "Miro.h"
#include "Object.h"
#include "PointLight.h"
#include "BVH.h"
#include "Texture.h"
#include "PhotonMap.h"

class Camera;
class Image;

struct HitPoint
{
	Vector3 position;
	Vector3 normal;
	Vector3 dir;
	float brdf;
	//int pixel_x, pixel_y;
	//float pixel_wgt;
	double radius;
	int accPhotons;
	int newPhotons;
	long double accFlux;
	long double newFlux;
    float scaling;

	HitPoint()
		:position(0.f), normal(0.f), dir(0.f), brdf(1.f), radius(0.f), accPhotons(0), newPhotons(0), accFlux(0.f), newFlux(0.f), scaling(1)
	{}
};

typedef std::vector<HitPoint*> HitPoints;

class Scene
{
public:
	Scene() 
		: m_photonMap(PhotonsPerLightSource*TRACE_DEPTH*MaxLights+MaxLights*10000), m_causticMap(CausticPhotonsPerLightSource*TRACE_DEPTH*MaxLights+MaxLights*10000), m_environment(0), m_bgColor(Vector3(0.0f)), m_photonsEmitted(0)
	{}
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

	void addHitPoint (HitPoint* hp)		{m_hitpoints.push_back(hp);}
	const HitPoints* hitpoints() const	{return &m_hitpoints;}

    void generatePhotonMap();
    void ProgressivePhotonPass();

    void preCalc();
    void openGL(Camera *cam);

    void raytraceImage(Camera *cam, Image *img);
    bool trace(HitInfo& minHit, const Ray& ray,
               float tMin = 0.0f, float tMax = MIRO_TMAX) const;
	bool traceScene(const Ray& ray, Vector3& shadeResult, int depth);

    void tracePhotons();
    void traceCausticPhotons();
    void traceProgressivePhotons();
    int tracePhoton(const Vector3& position, const Vector3& direction, const Vector3& power, int depth, bool bCausticRay=false);
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
    BVH m_bvh;
    Lights m_lights;
	HitPoints m_hitpoints;
    Texture * m_environment; //Environment map
    Vector3 m_bgColor;       //Background color (for when environment map is not available)

    static const int MaxLights = 10;

    static const int PhotonsPerLightSource = 100000;
    static const int CausticPhotonsPerLightSource = 100000;

	long int m_photonsEmitted;
};

extern Scene * g_scene;

#endif // CSE168_SCENE_H_INCLUDED
