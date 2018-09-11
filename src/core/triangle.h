#ifndef triangle_h__
#define triangle_h__

#include "vectors.h"

// triangle object

class Ray;
class AABB;
class HitInfo;

class TriangleObject
{
    int mat_id;

public:
    vec3 a, b, c;
    vec3 n, na, nb, nc;
    vec3 m_vertex_alpha;

    int materialID() const;
    void setMaterialID(int id);
    void setVertexAlpha(const vec3& alpha);

    TriangleObject();
    TriangleObject(vec3 a_, vec3 b_, vec3 c_);
    TriangleObject(vec3 a_, vec3 b_, vec3 c_, vec3 na_, vec3 nb_, vec3 nc_);

    HitInfo intersect(const HitInfo &hit, const Ray &ray) const;
    vec3 getCenter() const;
    float getArea() const;
    vec3 getFaceNormal() const;
    vec3 samplePoint() const;
    vec3 samplePoint(float& u, float& v) const;
    vec3 normal(float u, float v) const;
    vec3 point(float u, float v) const;
    vec3 consistentNormal(const vec3& wi, float u, float v,
        vec3 *refl = nullptr, float *cos_wi = nullptr) const;
    AABB boundingBox() const;
};

#endif // triangle_h__
