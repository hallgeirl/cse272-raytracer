#ifndef __ASSIGNMENT3_H_
#define __ASSIGNMENT3_H_
#include <map>
#include <cstring>
#include "Material.h"
#include "Utility.h"

void makeTask3Scene();
void a3task1();
void a3task2();
void a3task3();
void a3hacker1();
void BuildSquare(const Vector3& min, const Vector3& max, const Vector3& normal, const Material* mat);

#define W 512
#define H 512

//for bidirectional path tracing
struct hit_point
{
    Vector3 x;
    Vector3 contrib;
    Vector3 N;

    hit_point() : contrib(0), x(0), N(0)
    {}
};

struct sample
{
    Vector3 value;
    bool hit;
    bool direct;
    sample() 
		:hit(false), value(0), direct(false)
	{} 
};

typedef std::map<long, sample> sample_map;

struct samples
{
    sample_map X;
    double p;
    long n;

    samples() { n = 0; }
};

//Number of random directions
//Regular metropolis
#define PATH_LENGTH 9

//For bidirectional path tracing
#define LIGHT_PATH_LENGTH 4
#define EYE_PATH_LENGTH 5


class path
{
public:
    path() : F(0) { I = 0; }
    double u[(PATH_LENGTH+1)*2+2]; //position, [theta, phi]*9 (one direction for starting point, 8 directions for subsequent points). Starting point direction may or may not be used.
    double I;
    Vector3 F;

    void init_random()
    {
        for (int i = 0; i < (PATH_LENGTH+1)*2+2; i++)
        {
            u[i] = frand();
        }
    }
    void print()
    {
        printf("Path: I=%lf, u=[", I);
        for (int i = 0; i < 2+(PATH_LENGTH+1)*2; i++)
        {
            printf("%s%lf", (i == 0 ?"" : ", "), u[i]);
        }
        printf("]\n");
    }
}; 

inline void path_copy(path& to, path& from)
{
    memcpy(&to, &from, sizeof(path));
}

//static const double p_large = 1;
static const double p_large = 0.005;
#endif
