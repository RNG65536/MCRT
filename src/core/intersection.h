#ifndef intersection_h__
#define intersection_h__

#include "vectors.h"

// intersection record

class TriangleObject;

class HitInfo
{
private:
    // intersection information
    float t;//parametric distance of intersection
    vec3 nl;//normal vector at hit point
    float u, v;
    const TriangleObject *objectTriangle;
    int prim_id = -1;
    int mat_id = -1;

public:
    HitInfo();
//     HitInfo(float t_, vec3 nl_, const TriangleObject* object_);
    HitInfo(float t_, float u_, float v_, const vec3& nl_, const TriangleObject* object_);
    operator bool() const;

//     vec3 getColor() const;
//     float getEmission() const;
//     MaterialType getMaterialType() const;
//     AbstractBSDF *getMaterialBSDF() const;
//     const Material& getMaterial() const;
    void setDistance(float dist);
    float getDistance() const;
    void setShadingNormal(const vec3& n);
    vec3 getShadingNormal() const;
    vec3 getFaceNormal() const;
    vec3 getConsistentNormal(const vec3& wi,
        vec3 *refl = nullptr, float *cos_wi = nullptr) const;
    const TriangleObject *getObject() const;
    int materialID() const;
    void setMaterialID(int id);
    int primitiveID() const;
    void setPrimitiveID(int id);
};

#endif // intersection_h__
