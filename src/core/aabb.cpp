#include <cmath>
#include "aabb.h"
#include "constants.h"
#include "triangle.h"

AABB::AABB()
{
    m_min = vec3(NUM_INFINITY, NUM_INFINITY, NUM_INFINITY);
    m_max = vec3(-NUM_INFINITY, -NUM_INFINITY, -NUM_INFINITY);
}

AABB::AABB(const AABB& aabb) : m_min(aabb.m_min), m_max(aabb.m_max)
{
}

AABB::AABB(const vec3& min, const vec3& max) : m_min(min), m_max(max)
{
}

AABB& AABB::operator=(const AABB& aabb)
{
    m_min = aabb.m_min;
    m_max = aabb.m_max;
    return *this;
}

void AABB::enclose(const vec3& p)
{
    m_min = f_min(m_min, p);
    m_max = f_max(m_max, p);
}

void AABB::enclose(const AABB& b)
{
    m_min = f_min(m_min, b.m_min);
    m_max = f_max(m_max, b.m_max);
}

vec3 AABB::diagonal() const
{
    return m_max - m_min;
}

vec3 AABB::center() const
{
    return (m_max + m_min) * 0.5f;
}

const vec3& AABB::min() const
{
    return m_min;
}

const vec3& AABB::max() const
{
    return m_max;
}

AABB enclose(const AABB& a, const AABB& b)
{
    return AABB(f_min(a.min(), b.min()), f_max(a.max(), b.max()));
}

AABB enclose(const AABB& a, const vec3& point)
{
    return AABB(f_min(a.min(), point), f_max(a.max(), point));
}

AABB getAABB(const TriangleObject& object)
{
    return AABB(f_min(f_min(object.a, object.b), object.c),
                f_max(f_max(object.a, object.b), object.c));
}
