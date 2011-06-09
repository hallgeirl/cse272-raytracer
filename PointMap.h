//----------------------------------------------------------------------------
// pointmap.cc
// A copy implementation of the photon map data structure
//
// Henrik Wann Jensen - February 2001
//----------------------------------------------------------------------------

#ifndef __POINT_MAP_H__
#define __POINT_MAP_H__

#include "Vector3.h"

struct Point
{

	Vector3 position;
	Vector3 normal;
	Vector3 dir;
	float brdf;
	double radius;
	int accPhotons;
	int newPhotons;
	long double accFlux;
	long double newFlux;
    float scaling;
	bool bLight;
    short plane;                  // splitting plane for kd-tree
	unsigned char theta, phi;
	int i, j;

	Point()
		:position(0.f), normal(0.f), dir(0.f), brdf(1.f), radius(0.f), accPhotons(0), newPhotons(0), accFlux(0.f), newFlux(0.f), 
			scaling(1), bLight(true), theta(0), phi(0), i(0), j(0)
	{}

	Point(const Vector3& inPosition, const Vector3& inNormal, const Vector3& inDir, const float inBRDF, const float inRadius)
		:position(inPosition), normal(inNormal), dir(inDir), brdf(inBRDF), radius(inRadius), accPhotons(0), newPhotons(0), accFlux(0.f), newFlux(0.f), 
			scaling(1), bLight(false), theta(0), phi(0), i(0), j(0)
	{}
};

/* This structure is used only to locate the
 * nearest photons
*/
//******************************
typedef struct NearestPoints {
//******************************
    int max;
    int found;
    int got_heap;
    Vector3 pos;
    float *dist2;
    Point **index;
} NearestPoints;


/* This is the Point_map class
 */
//*****************
class Point_map {
//*****************
public:
  Point_map( int max_point );
  ~Point_map();

  Point* store(
    const Vector3& pos,            // measurement point direction
    const Vector3& normal,         // surface normal
    const Vector3& dir,			   // point direction
	const float radius,			   // initial radius
	const float brdf,			   // surface brdf
	const bool bLight);			   // if valid point

  void empty();

  void balance(void);              // balance the kd-tree (before use!)

  void find_points(
	NearestPoints* np,
    const Vector3& pos,            // surface position
    const Vector3& normal,         // surface normal at pos
    const float max_dist,          // max distance to look for photons
    const int npoints) const; // if the flux should be normalized by density    

  void locate_points(
    NearestPoints *const np,      // np is used to locate the photons
    const int index, const Vector3& normal ) const;       // call with index = 1

private:

  void balance_segment(
    Point **pbal,
    Point **porg,
    const int index,
    const int start,
    const int end );

  void median_split(
    Point **p,
    const int start,
    const int end,
    const int median,
    const int axis );
  
  Point *points;

  int stored_points;
  int half_stored_points;
  int max_points;
  int prev_scale;

  float costheta[256];
  float sintheta[256];
  float cosphi[256];
  float sinphi[256];
  
  Vector3 bbox_min;		// use bbox_min;
  Vector3 bbox_max;		// use bbox_max;
};

#endif
