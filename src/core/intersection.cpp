#include "constants.h"
#include "intersection.h"
#include "triangle.h"

HitInfo::operator bool() const
{
//     return objectTriangle != nullptr;
    return prim_id != -1;
}

// HitInfo::HitInfo(float t_, vec3 nl_, const TriangleObject* object_) : t(t_), nl(nl_), objectTriangle(object_)
HitInfo::HitInfo(float t_, float u_, float v_, const vec3& nl_, const TriangleObject* object_) : t(t_), u(u_), v(v_), nl(nl_), objectTriangle(object_)
{

}

// HitInfo::HitInfo() : t(NUM_INFINITY), nl(vec3(0, 0, 0)),
// objectTriangle(nullptr)
HitInfo::HitInfo() : t(NUM_INFINITY), u(0), v(0),
objectTriangle(nullptr)
{

}

// vec3 HitInfo::getColor() const
// {
//     return  objectTriangle->getMaterial().getColor();
// }
// float HitInfo::getEmission() const{
//     return  objectTriangle->getMaterial().getEmission();
// }
// MaterialType HitInfo::getMaterialType() const{
//     return  objectTriangle->mat.getType();
// }

// const Material& HitInfo::getMaterial() const
// {
//     return objectTriangle->getMaterial();
// }

vec3 HitInfo::getShadingNormal() const
{
    return nl;
//     return objectTriangle->normal(u, v);
}

float HitInfo::getDistance() const
{
    return t;
}

const TriangleObject * HitInfo::getObject() const
{
    return objectTriangle;
}

vec3 HitInfo::getConsistentNormal(const vec3& wi, vec3 *refl, float *cos_wi) const
{
    return objectTriangle->consistentNormal(wi, u, v, refl, cos_wi);
}

vec3 HitInfo::getFaceNormal() const
{
    return objectTriangle->getFaceNormal();
}

void HitInfo::setMaterialID(int id)
{
    mat_id = id;
}

int HitInfo::materialID() const
{
    return mat_id;
}

void HitInfo::setPrimitiveID(int id)
{
    prim_id = id;
}

int HitInfo::primitiveID() const
{
    return prim_id;
}

void HitInfo::setShadingNormal(const vec3& n)
{
    nl = n;
}

void HitInfo::setDistance(float dist)
{
    t = dist;
}

// AbstractBSDF * Hit_CUDA::getMaterialBSDF() const
// {
//     return objectTriangle->mat.getBSDF();
// }
