#pragma once

#include "vectors.h"

// intersection record

class TriangleObject;

class HitInfo
{
public:
    HitInfo();
    HitInfo(float                 t,
            float                 u,
            float                 v,
            const vec3&           nl,
            const TriangleObject* hit_object);
    operator bool() const;

    void reset();

    float distance() const;

    const vec3& shadingNormal() const;
    const vec3& faceNormal() const;

    vec3 consistentNormal(const vec3& wi,
                          vec3*       refl = nullptr,
                          float*      cos_wi = nullptr) const;

    const TriangleObject* triangleObject() const;
    int                   materialID() const;
    void                  setMaterialID(int id);
    int                   primitiveID() const;
    void                  setPrimitiveID(int id);

    bool operator<(const HitInfo& b) const;

public:
    float m_t;       // ray distance
    vec3  m_nl;      // shading normal
    vec3  m_fn;      // face normal
    float m_u, m_v;  // barycentric coordinates

    int m_objID = -1;
    int m_matID = -1;

    const TriangleObject* m_hit_object = nullptr;
};
