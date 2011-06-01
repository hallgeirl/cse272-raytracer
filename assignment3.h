#ifndef __ASSIGNMENT3_H_
#define __ASSIGNMENT3_H_
#include <map>
#include <cstring>
#include "Material.h"

void makeTask3Scene();
void a3task1();
void a3task2();
void a3task3();
void BuildSquare(const Vector3& min, const Vector3& max, const Vector3& normal, const Material* mat);

struct sample
{
    double value, dist2, costheta;
    bool hit;
    bool direct;
    int nrays;
    sample() 
		:nrays(0), hit(false)
	{} 
};

typedef std::map<long, sample> sample_map;

struct samples
{
    sample_map X;
    double p;
    long n;
    bool isAreaSample;

    samples() { n = 0; }
};

struct path
{
    double u[3]; //position, theta, phi
    double I;
}; 

inline void path_copy(path& to, path& from)
{
    memcpy(&to, &from, sizeof(path));
}

static const double p_large = 0.6;
static const double p_pos = 1;
static const double p_angle = 1;
#endif
