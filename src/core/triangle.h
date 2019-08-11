#pragma once

#include "intersection.h"
#include "ray.h"
#include "vectors.h"

// triangle object

class AABB;

class TriangleObject
{
public:
    int  materialID() const;
    void setMaterialID(int id);
    void setVertexAlpha(const vec3& alpha);  // for consistent normal

    const int& primitiveID() const;
    void setPrimitiveID(int id);

    TriangleObject();
    TriangleObject(vec3 a, vec3 b, vec3 c);
    TriangleObject(vec3 a, vec3 b, vec3 c, vec3 na, vec3 nb, vec3 nc);

    float area() const;
    vec3  baryCenter() const;
    vec3  normal() const;
    vec3  consistentNormal(const vec3& wi,
                           float       u,
                           float       v,
                           vec3*       refl = nullptr,
                           float*      cos_wi = nullptr) const;

    // interpolation using barycentric coordinates
    vec3 normal(float u, float v) const;
    vec3 point(float u, float v) const;

    // uniform sampling
    vec3 samplePoint() const;
    vec3 samplePoint(float& u, float& v) const;
    vec3 samplePoint(float& u, float& v, float r0, float r1) const;

    AABB boundingBox() const;

    bool intersect(Ray& ray, HitInfo& hit_info) const;

    std::string toString() const;

    static void setIntersectionEpsilon(float eps);

public:
    vec3 a, b, c;         // vertices
    vec3 n, na, nb, nc;   // face normal and vertex normals
    vec3 m_vertex_alpha;  // for consistent normal

private:
    int m_matID;
    int m_objID;

    static float s_intersection_epsilon;
};

typedef TriangleObject Triangle;