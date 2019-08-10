#pragma once

#include "intersection.h"
#include "ray.h"
#include "vectors.h"

class TriangleObject;
class BoundedTriangle;

// axis aligned bounding box

class AABB
{
public:
    AABB();
    AABB(const float3& min, const float3& max);
    AABB(const AABB& aabb);
    AABB& operator=(const AABB& aabb);

    const float3& min() const;
    const float3& max() const;

    void reset();

    void   expandBy(const AABB& box);
    void   expandBy(const float3& point);
    float3 center() const;
    float3 diagonal() const;
    float  surfaceArea() const;
    float  surfaceAreaCost() const;

    bool isValid() const;
    int  largestDimension() const;

    bool intersect(Ray& ray, float& tnear) const;

    // MUST use this for AABB culling because
    // the ray's origin can be inside the AABB
    bool overlap(const Ray& ray) const;

    void intersectBy(const AABB& box);

    AABB intersect(const AABB& box) const;

    AABB clipAxis(float min, float max, int axis) const;
    AABB clipX(float min, float max) const;
    AABB clipY(float min, float max) const;
    AABB clipZ(float min, float max) const;

    std::string toString() const;

private:
    float3 m_min;
    float3 m_max;
};

AABB expandBy(const AABB& a, const AABB& b);
AABB expandBy(const AABB& a, const float3& point);
AABB getAABB(const TriangleObject& object);
