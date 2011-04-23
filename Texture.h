#ifndef _TEXTURE_H_
#define _TEXTURE_H_

#define NO_FREEIMAGE
#include "Phong.h"
#include <Perlin.h>
#include <string>
#include <iostream>
#include <vector>
#include <list>
#include <stdlib.h>
#ifndef NO_FREEIMAGE
#include <FreeImage.h>
#endif
#include <Worley.h>
#include <cstdlib>
#include <cmath>
#include "Utility.h"


//Generate turbulent perlin noise at (x,y,z).
//The returned number is between -1 and 1.
//frequencyIncrease indicates how quickly the frequency increases over multiple turbulence iterations.
//amplitudeFalloff indicates how quickly the amplitude of the higher frequency noise values decrease.
//iterations determines the number of iterations to use for turbulence.
inline float generateNoise(float x, float y, float z, float initialFrequency, float frequencyIncrease, float amplitudeFalloff, int iterations)
{
    float amplitude = 1;
    float frequency = initialFrequency;
    float value = 0;
    float max_val = 0;
    
    for (int i = 0; i < iterations; i++)
    {
        value += amplitude * PerlinNoise::noise(x*frequency, y*frequency, z*frequency);
        max_val += amplitude;
        
        frequency *= frequencyIncrease;
        amplitude *= amplitudeFalloff;
    }
  
    return value/max_val;
}

//Used for the grid in cellular textures
typedef struct gridcell_s
{
    std::vector<tex_coord2d_t> points;
    void addPoint(tex_coord2d_t point) { points.push_back(point); }
} gridcell_t;

class Grid
{
public:
    Grid(int gridWidth, int gridHeight);
    void addPoint(tex_coord2d_t point);

    const std::vector<tex_coord2d_t> &getPoints(int cellI, int cellJ) const { return m_grid[cellI][cellJ].points; }
    int getWidth() const { return m_gridWidth; }
    int getHeight() const { return m_gridHeight; }
private:
    gridcell_t ** m_grid;
    int m_gridWidth, m_gridHeight;
};

class Texture
{
public:
    virtual float bumpHeight2D(const tex_coord2d_t & coords) const { return 0; }
    virtual float bumpHeight3D(const tex_coord3d_t & coords) const { return 0; }
	virtual LookupCoordinates GetLookupCoordinates() const = 0;
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const { return Vector3(0,0,0); } // Look up the color value for a specified position
    virtual Vector3 lookup3D(const tex_coord3d_t & coords) const { return Vector3(0,0,0); } // For 3D textures
    virtual Vector3 lowresLookup2D(const tex_coord2d_t & coords) const { return lookup2D(coords); }
    virtual Vector3 lowresLookup3D(const tex_coord2d_t & coords) const { return lookup2D(coords); }
};

class Texture2D : public Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const {return UV; }
};

class Texture3D : public Texture
{
public:
	virtual LookupCoordinates GetLookupCoordinates() const {return UVW; }
};

class CellularTexture2D : public Texture2D
{
public:
    CellularTexture2D(int points, int gridWidth, int gridHeight);
    float * getClosestDistances(const tex_coord2d_t &point, int n) const; //Get the n closest distances from a given point
    
    //Populate the grid with grid points. Override to control distribution of points.
    virtual void populateGrid(int points); 
    //Typically, this function combines the values gotten from getClosestDistances() and maps it to a color value.
    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const; 
    
protected:
    std::vector<tex_coord2d_t> m_points;
    
    Grid m_grid;
};

class CheckerBoardTexture : public Texture2D
{
protected:
    float m_scale;
    Vector3 m_color1, m_color2;
public:
    CheckerBoardTexture(Vector3 color1 = Vector3(1), Vector3 color2 = Vector3(0), float scale=1) 
    { 
        m_color1 = color1; 
        m_color2 = color2; 
        m_scale = scale; 
    }

    virtual Vector3 lookup2D(const tex_coord2d_t & coords) const
    {
        float u = std::abs(m_scale*(coords.u));
        float v = std::abs(m_scale*(coords.v));
        if (coords.u < 0) u += m_scale;
        if (coords.v < 0) v += m_scale;
        return ((int)((int)u+(int)v) % 2 == 0 ? m_color1 : m_color2); 
    }
};

#ifndef NO_FREEIMAGE

//Texture loaded from a file
class LoadedTexture : public Texture2D
{
public:
    LoadedTexture(std::string filename);
	~LoadedTexture();

    Vector3 lookup(const tex_coord2d_t & coords, bool lowres) const;
    Vector3 lookup2D(const tex_coord2d_t & coords) const { return lookup(coords, false); }
    Vector3 lowresLookup2D(const tex_coord2d_t & coords) const { return lookup(coords, true); }
    
protected:
    static Vector3 getPixel(FIBITMAP* bm, int x, int y);     //Get the tonemapped pixel
    static void setPixel(FIBITMAP* bm, Vector3& value, int x, int y);
    float tonemapValue(float val) const;

	FIBITMAP* m_bitmap; 
    FIBITMAP* m_lowres;
    static const int LOWRES_WIDTH = 24;	
	float m_maxIntensity;
};

#endif

//Shading model that also does textures
class TexturedPhong : public Phong
{
public:
    TexturedPhong(Texture * texture,
            			const Vector3 & specularColor = Vector3(0),    //Reflectivity for each color. 1 is fully reflective, 0 is fully non-reflective.
			            const Vector3 & transparentColor = Vector3(0), //Transparency for each color. 1 is fully transparent (refracting according to refractIndex), 0 is fully opaque.
			            const float shinyness = 1.0f,
			            const float refractIndex = 1);
	virtual LookupCoordinates GetLookupCoordinates() const { return m_texture->GetLookupCoordinates(); }

    virtual Vector3 diffuse2D(const tex_coord2d_t & texture_coords) const;
    virtual Vector3 diffuse3D(const tex_coord3d_t & texture_coords) const;
    virtual float   bumpHeight2D(const tex_coord2d_t & texture_coords) const;
    virtual float   bumpHeight3D(const tex_coord3d_t & texture_coords) const;

protected:
    Texture * m_texture;
};

#endif
