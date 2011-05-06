#ifndef __ASSIGNMENT2_H_
#define __ASSIGNMENT2_H_
#include <map>

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

#endif
