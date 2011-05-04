#include "TriangleMesh.h"
#include "Triangle.h"
#include "Scene.h"

void 
createTriangleMesh(const char* filename, Material *mat, Scene* scene, Vector3 position, float rotY, Vector3 scale)
{
	TriangleMesh * mesh = new TriangleMesh();
    Matrix4x4 m_rot(Vector4(cos(rotY),0,-sin(rotY),0), Vector4(0,1,0,0), Vector4(sin(rotY),0,cos(rotY),0), Vector4(0, 0, 0, 1));
    Matrix4x4 m_trans(Vector4(1,0,0,0), Vector4(0,1,0,0), Vector4(0,0,1,0), Vector4(position.x, position.y, position.z, 1));
    Matrix4x4 m_scale(Vector4(scale.x,0,0,0), Vector4(0,scale.y,0,0), Vector4(0,0,scale.z,0), Vector4(0, 0, 0, 1));
	mesh->load(filename, m_trans*m_rot*m_scale);
    for (int i = 0; i < mesh->numTris(); i++)
    {
		Triangle *tri = new Triangle();
		tri->setMesh(mesh);
		tri->setIndex(i);
		tri->setMaterial(mat);
		g_scene->addObject(tri);
    }
}


TriangleMesh::TriangleMesh() :
    m_normals(0),
    m_vertices(0),
    m_texCoords(0),
    m_normalIndices(0),
    m_vertexIndices(0),
    m_texCoordIndices(0)
{

}

TriangleMesh::~TriangleMesh()
{
    delete [] m_normals;
    delete [] m_vertices;
    delete [] m_texCoords;
    delete [] m_normalIndices;
    delete [] m_vertexIndices;
    delete [] m_texCoordIndices;
}

void TriangleMesh::setVertex(int index, const Vector3 &v)
{
    m_vertices[index] = v;

    #ifdef __SSE4_1__
    m_SSEvertices[index] = _mm_set_ps(v.x, v.y, v.z, 0.0f);
    #endif
}

void TriangleMesh::setNormal(int index, const Vector3 &n)
{
    m_normals[index] = n;

    #ifdef __SSE4_1__
    m_SSEnormals[index] = _mm_set_ps(n.x, n.y, n.z, 0.0f);
    #endif
}
