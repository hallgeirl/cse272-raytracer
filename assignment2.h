#ifndef __ASSIGNMENT2_H_
#define __ASSIGNMENT2_H_
#include <map>
#include <cstring>

void makeTask2Scene();
void a2task1();
void a2task2();
void a2task3();
void a2task4();

struct sample
{
    double value, dist2, costheta;
    bool hit;
    bool direct;
    int nrays;
    sample() { nrays = 0; } 
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
