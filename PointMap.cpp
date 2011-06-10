//----------------------------------------------------------------------------
// photonmap.cc
// An example implementation of the photon map data structure
//
// Henrik Wann Jensen - February 2001
//----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Vector3.h"
#include "PointMap.h"
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

/* This is the constructor for the photon map.
 * To create the photon map it is necessary to specify the
 * maximum number of photons that will be stored
 */
//************************************************
Point_map :: Point_map( const int max_point )
    //************************************************
{
    stored_points = 0;
    half_stored_points = 0;
    prev_scale = 1;
    max_points = max_point;

    points = (Point*)malloc( sizeof( Point ) * ( max_points+1 ) );

    if (points == NULL) {
        fprintf(stderr,"Out of memory initializing point map\n");
        exit(-1);
    }

    bbox_min.x = bbox_min.y = bbox_min.z = 1e8f;
    bbox_max.x = bbox_max.y = bbox_max.z = -1e8f;

    //----------------------------------------
    // initialize direction conversion tables
    //----------------------------------------

    for (int i=0; i<256; i++) {
        double angle = double(i)*(1.0/256.0)*M_PI;
        costheta[i] = cos( angle );
        sintheta[i] = sin( angle );
        cosphi[i]   = cos( 2.0*angle );
        sinphi[i]   = sin( 2.0*angle );
    }
}


//*************************
Point_map :: ~Point_map()
    //*************************
{
    free( points );
}


/* irradiance_estimate computes an irradiance estimate
 * at a given surface position
 */
//**********************************************
void Point_map :: find_points(
        NearestPoints* np,
        const Vector3& pos,            // surface position
        const float max_dist,          // max distance to look for photons
        const int npoints) const     // number of photons to use
//**********************************************
{
    np->pos = pos; 
    np->max = npoints;
    np->found = 0;
    np->got_heap = 0;
    np->dist2[0] = max_dist*max_dist;

    // locate the nearest points
    locate_points( np, 1);
}


/* locate_points finds the nearest photons in the
 * photon map given the parameters in np
 */
//******************************************
void Point_map :: locate_points(
        NearestPoints *const np,
        const int index) const
//******************************************
{
    Point *p = &points[index];
    float dist1;

//    for (int i = 0; i < stored_points; i++)
//    {
//        Point* pp = &points[i];
//        if (pp->j == 176 && pp->i == 303)
//            cout << endl << "IT'S HERE!!!" << endl;
//    }
//    cout << stored_points << endl;
//    for (int i = 0; i < stored_points; i++)
//    {
//        Point* pp = &points[i];
//        if (pp->j == 176 && pp->i == 303)
//            cout << endl << "IT'S HERE!!!" << endl;
//    }

    if (index<half_stored_points) {

        switch(p->plane)
        {
            case 0:dist1 = np->pos.x - p->position.x; break;
            case 1:dist1 = np->pos.y - p->position.y; break;
            case 2:dist1 = np->pos.z - p->position.z; break;
        }
        //dist1 = np->pos[ p->plane ] - p->pos[ p->plane ];

        if (dist1>0.0) { // if dist1 is positive search right plane
            locate_points( np, 2*index+1);
            if ( dist1*dist1 < np->dist2[0] )
                locate_points( np, 2*index);
        } else {         // dist1 is negative search left first
            locate_points( np, 2*index);
            if ( dist1*dist1 < np->dist2[0] )
                locate_points( np, 2*index+1);
        }
    }

    // compute squared distance between current photon and np->pos

    dist1 = p->position.x - np->pos.x;
    float dist2 = dist1*dist1;
    dist1 = p->position.y - np->pos.y;
    dist2 += dist1*dist1;
    dist1 = p->position.z - np->pos.z;
    dist2 += dist1*dist1;


    if ( dist2 < np->dist2[0]) {
        // we found a photon :) Insert it in the candidate list

        if ( np->found < np->max ) {
            // heap is not full; use array
            np->found++;
            np->dist2[np->found] = dist2;
            np->index[np->found] = p;
        } else {
            int j,parent;

            if (np->got_heap==0) { // Do we need to build the heap?
                // Build heap
                float dst2;
                Point *point;
                int half_found = np->found>>1;
                for ( int k=half_found; k>=1; k--) {
                    parent=k;
                    point = np->index[k];
                    dst2 = np->dist2[k];
                    while ( parent <= half_found ) {
                        j = parent+parent;
                        if (j<np->found && np->dist2[j]<np->dist2[j+1])
                            j++;
                        if (dst2>=np->dist2[j])
                            break;
                        np->dist2[parent] = np->dist2[j];
                        np->index[parent] = np->index[j];
                        parent=j;
                    }
                    np->dist2[parent] = dst2;
                    np->index[parent] = point;
                }
                np->got_heap = 1;
            }

            // insert new point into max heap
            // delete largest element, insert new and reorder the heap

            parent=1;
            j = 2;
            while ( j <= np->found ) {
                if ( j < np->found && np->dist2[j] < np->dist2[j+1] )
                    j++;
                if ( dist2 > np->dist2[j] )
                    break;
                np->dist2[parent] = np->dist2[j];
                np->index[parent] = np->index[j];
                parent = j;
                j += j;
            }
            np->index[parent] = p;
            np->dist2[parent] = dist2;

            np->dist2[0] = np->dist2[1];
        }
    }
}


/* store puts a point into the flat array that will form
 * the final kd-tree.
 *
 * Call this function to store a point.
 */
//***************************
Point* Point_map :: store(
        const Vector3& pos,            
        const Vector3& normal,         
        const Vector3& dir,			   
        const float radius,			   
        const float brdf,
        const bool bLight,
        int x, int y) //x,y is pixel position
//***************************
{
    if (stored_points>=max_points)
        return NULL;

//    if (bLight)
//        printf("Invalid point created in point map");

    stored_points++;
    //Point *const node = &points[stored_points];
    Point *node = &points[stored_points];

    node->position = pos;
    node->normal = normal;
    node->dir = dir;
    node->radius = radius;
    node->brdf = brdf;
    node->bLight = bLight;
    node->accFlux = 0.f;
    node->accPhotons = 0;
    node->newFlux = 0.f;
    node->newPhotons = 0;
    node->scaling = 1.f;
    node->i = y;
    node->j = x;



    if (node->position.x < bbox_min.x)
        bbox_min.x = node->position.x;
    if (node->position.x > bbox_max.x)
        bbox_max.x = node->position.x;

    if (node->position.y < bbox_min.y)
        bbox_min.y = node->position.y;
    if (node->position.y > bbox_max.y)
        bbox_max.y = node->position.y;

    if (node->position.z < bbox_min.z)
        bbox_min.z = node->position.z;
    if (node->position.z > bbox_max.z)
        bbox_max.z = node->position.z;

    int theta = int( acos(dir.z)*(256.0/M_PI) );
    if (theta>255)
        node->theta = 255;
    else
        node->theta = (unsigned char)theta;

    int phi = int( atan2(dir.y,dir.x)*(256.0/(2.0*M_PI)) );
    if (phi>255)
        node->phi = 255;
    else if (phi<0)
        node->phi = (unsigned char)(phi+256);
    else
        node->phi = (unsigned char)phi;

    return node;
}

/* empty the  flat array 
 * and reset kd-tree bounds
 *
 * Call this function to empty the photon map.
 */
//***************************
void Point_map :: empty()
{
    stored_points = 0;
    half_stored_points = 0;
    prev_scale = 1;

    bbox_min.x = bbox_min.y = bbox_min.z = 1e8f;
    bbox_max.x = bbox_max.y = bbox_max.z = -1e8f;
}

/* balance creates a left balanced kd-tree from the flat photon array.
 * This function should be called before the photon map
 * is used for rendering.
 */
//******************************
void Point_map :: balance(void)
    //******************************
{
    //  cout << "Balancing tree. Stored photons " << stored_photons << endl;
    if (stored_points>1) {
        // allocate two temporary arrays for the balancing procedure
        Point **pa1 = (Point**)malloc(sizeof(Point*)*(stored_points+1));
        Point **pa2 = (Point**)malloc(sizeof(Point*)*(stored_points+1));


        for (int i=0; i<=stored_points; i++)
        {
            pa2[i] = &points[i];
        }

        balance_segment( pa1, pa2, 1, 1, stored_points );
        free(pa2);

        // reorganize balanced kd-tree (make a heap)
        int d, j=1, foo=1;
        Point foo_point = points[j];

        for (int i=1; i<=stored_points; i++) {
            d=pa1[j]-points;
            pa1[j] = NULL;
            if (d != foo)
                points[j] = points[d];
            else {
                points[j] = foo_point;

                if (i<stored_points) {
                    for (;foo<=stored_points; foo++)
                        if (pa1[foo] != NULL)
                            break;
                    foo_point = points[foo];
                    j = foo;
                }
                continue;
            }
            j = d;
        }
        free(pa1);
    }

    half_stored_points = stored_points/2-1;
}


#define swap(ph,a,b) { Point *ph2=ph[a]; ph[a]=ph[b]; ph[b]=ph2; }
float DetermineAxis(Point *p, const int axis)
{
    switch(axis)
    {
        case 0: return p->position.x; break;
        case 1: return p->position.y; break;
        case 2: default: return p->position.z; break;
    }
}

float DetermineAxis(const Vector3& v, const int axis)
{
    switch(axis)
    {
        case 0: return v.x; break;
        case 1: return v.y; break;
        case 2: default: return v.z; break;
    }
}

// median_split splits the photon array into two separate
// pieces around the median with all photons below the
// the median in the lower half and all photons above
// than the median in the upper half. The comparison
// criteria is the axis (indicated by the axis parameter)
// (inspired by routine in "Algorithms in C++" by Sedgewick)
//*****************************************************************
void Point_map :: median_split(
        Point **p,
        const int start,               // start of photon block in array
        const int end,                 // end of photon block in array
        const int median,              // desired median number
        const int axis )               // axis to split along
//*****************************************************************
{
    int left = start;
    int right = end;

    while ( right > left ) {
        const float v = DetermineAxis(p[right], axis);
        int i=left-1;
        int j=right;
        for (;;) {
            while ( DetermineAxis(p[++i], axis) < v )
                ;
            while ( DetermineAxis(p[--j], axis) > v && j>left )
                ;
            if ( i >= j )
                break;
            swap(p,i,j);
        }

        swap(p,i,right);
        if ( i >= median )
            right=i-1;
        if ( i <= median )
            left=i+1;
    }
}


// See "Realistic image synthesis using Photon Mapping" chapter 6
// for an explanation of this function
//****************************
void Point_map :: balance_segment(
        Point **pbal,
        Point **porg,
        const int index,
        const int start,
        const int end )
//****************************
{
    //--------------------
    // compute new median
    //--------------------

    int median=1;
    while ((4*median) <= (end-start+1))
        median += median;

    if ((3*median) <= (end-start+1)) {
        median += median;
        median += start-1;
    } else	
        median = end-median+1;

    //--------------------------
    // find axis to split along
    //--------------------------

    int axis=2;
    if ((bbox_max.x-bbox_min.x)>(bbox_max.y-bbox_min.y) &&
            (bbox_max.x-bbox_min.x)>(bbox_max.z-bbox_min.z))
        axis=0;
    else if ((bbox_max.y-bbox_min.y)>(bbox_max.z-bbox_min.z))
        axis=1;

    //------------------------------------------
    // partition photon block around the median
    //------------------------------------------

    median_split( porg, start, end, median, axis );

    pbal[ index ] = porg[ median ];
    pbal[ index ]->plane = axis;

    //----------------------------------------------
    // recursively balance the left and right block
    //----------------------------------------------

    if ( median > start ) {
        // balance left segment
        if ( start < median-1 ) {
            const float tmp = DetermineAxis(bbox_max, axis);
            switch(axis)
            {
                case 0: bbox_max.x = DetermineAxis(pbal[index], axis); break;
                case 1: bbox_max.y = DetermineAxis(pbal[index], axis); break;
                case 2: bbox_max.z = DetermineAxis(pbal[index], axis); break;
            }
            balance_segment( pbal, porg, 2*index, start, median-1 );
            switch(axis)
            {
                case 0: bbox_max.x = tmp; break;
                case 1: bbox_max.y = tmp; break;
                case 2: bbox_max.z = tmp; break;
            }
        } else {
            pbal[ 2*index ] = porg[start];
        }
    }

    if ( median < end ) {
        // balance right segment
        if ( median+1 < end ) {
            const float tmp = DetermineAxis(bbox_min, axis);		
            switch(axis)
            {
                case 0: bbox_min.x = DetermineAxis(pbal[index], axis); break;
                case 1: bbox_min.y = DetermineAxis(pbal[index], axis); break;
                case 2: bbox_min.z = DetermineAxis(pbal[index], axis); break;
            }
            balance_segment( pbal, porg, 2*index+1, median+1, end );
            switch(axis)
            {
                case 0: bbox_min.x = tmp; break;
                case 1: bbox_min.y = tmp; break;
                case 2: bbox_min.z = tmp; break;
            }
        } else {
            pbal[ 2*index+1 ] = porg[end];
        }
    }	
}
