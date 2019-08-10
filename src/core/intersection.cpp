#include <cassert>
#include "constants.h"
#include "intersection.h"
#include "triangle.h"

HitInfo::HitInfo() : m_t(NUM_INFINITY), m_u(0), m_v(0), m_hit_object(nullptr)
{
}

HitInfo::HitInfo(
    float t, float u, float v, const vec3& nl, const TriangleObject* hit_object)
    : m_t(t), m_u(u), m_v(v), m_nl(nl), m_hit_object(hit_object)
{
    assert(hit_object);
    m_fn = hit_object->normal();
}

HitInfo::operator bool() const
{
    return m_objID != -1;
}

void HitInfo::reset()
{
    m_t = NUM_INFINITY;
    m_u = 0.0f;
    m_v = 0.0f;
    m_objID = -1;
    m_matID = -1;
    m_hit_object = nullptr;
}

float HitInfo::distance() const
{
    return m_t;
}

const vec3& HitInfo::shadingNormal() const
{
    return m_nl;
}

const vec3& HitInfo::faceNormal() const
{
    return m_fn;
}

vec3 HitInfo::consistentNormal(const vec3& wi, vec3* refl, float* cos_wi) const
{
    return m_hit_object->consistentNormal(wi, m_u, m_v, refl, cos_wi);
}

const TriangleObject* HitInfo::triangleObject() const
{
    return m_hit_object;
}

void HitInfo::setMaterialID(int id)
{
    m_matID = id;
}

int HitInfo::materialID() const
{
    return m_matID;
}

void HitInfo::setPrimitiveID(int id)
{
    m_objID = id;
}

bool HitInfo::operator<(const HitInfo& b) const
{
    return m_t < b.m_t;
}

int HitInfo::primitiveID() const
{
    return m_objID;
}
