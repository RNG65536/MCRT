#pragma once

#include "vectors.h"

class TriangleObject;

// axis aligned bounding box

class AABB
{
public:
    AABB();
    AABB(const vec3& min, const vec3& max);
    AABB(const AABB& aabb);
    AABB& operator=(const AABB& aabb);

    const vec3& min() const;
    const vec3& max() const;

    void enclose(const AABB& b);
    void enclose(const vec3& p);
    vec3 center() const;
    vec3 diagonal() const;

private:
    vec3 m_min, m_max;
};

AABB enclose(const AABB& a, const AABB& b);
AABB enclose(const AABB& a, const vec3& point);
AABB getAABB(const TriangleObject& object);
