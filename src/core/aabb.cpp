#include <cassert>
#include <cmath>
#include <sstream>
#include "aabb.h"
#include "constants.h"
#include "numeric.h"
#include "triangle.h"

// test if AABB and ray overlap in the range of [tmin, tmax]
// the actual hit distance is not needed
static bool rayAABBOverlap(const Ray&    r,
                           const float3& bounds_min,
                           const float3& bounds_max)
{
#if 0
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    const float3 bounds[2] = {bounds_min, bounds_max};

    tmin = (bounds[r.dir_sign.x].x - r.orig.x) * r.inv_dir.x;
    tmax = (bounds[1 - r.dir_sign.x].x - r.orig.x) * r.inv_dir.x;
    tymin = (bounds[r.dir_sign.y].y - r.orig.y) * r.inv_dir.y;
    tymax = (bounds[1 - r.dir_sign.y].y - r.orig.y) * r.inv_dir.y;

    if ((tmin > tymax) || (tymin > tmax))
    {
        return false;
    }
    if (tymin > tmin)
    {
        tmin = tymin;
    }
    if (tymax < tmax)
    {
        tmax = tymax;
    }

    tzmin = (bounds[r.dir_sign.z].z - r.orig.z) * r.inv_dir.z;
    tzmax = (bounds[1 - r.dir_sign.z].z - r.orig.z) * r.inv_dir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
    {
        return false;
    }
    if (tzmin > tmin)
    {
        tmin = tzmin;
    }
    if (tzmax < tmax)
    {
        tmax = tzmax;
    }

    // modification in order to handle cases where
    // ray origin is inside of the box
    if (tmin > r.tmax)
    {
        return false;
    }
    if (tmax < r.tmin)
    {
        return false;
    }
    return true;
#else
    // modified in order to handle cases where
    // ray origin is inside of the box.
    // prioritizes far culling.
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    const float3 bounds[2] = {bounds_min, bounds_max};

    tmin = (bounds[r.dir_sign.x].x - r.orig.x) * r.inv_dir.x;
    tymin = (bounds[r.dir_sign.y].y - r.orig.y) * r.inv_dir.y;
    tmax = (bounds[1 - r.dir_sign.x].x - r.orig.x) * r.inv_dir.x;
    tymax = (bounds[1 - r.dir_sign.y].y - r.orig.y) * r.inv_dir.y;

    if ((tmin > tymax) || (tymin > tmax))
    {
        return false;
    }

    if (tymin > tmin)
    {
        tmin = tymin;
    }
    // far culling
    if (tmin > r.tmax)
    {
        return false;
    }

    if (tymax < tmax)
    {
        tmax = tymax;
    }
    // near culling
    if (tmax < r.tmin)
    {
        return false;
    }

    tzmin = (bounds[r.dir_sign.z].z - r.orig.z) * r.inv_dir.z;
    tzmax = (bounds[1 - r.dir_sign.z].z - r.orig.z) * r.inv_dir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
    {
        return false;
    }

    if (tzmin > tmin)
    {
        tmin = tzmin;
    }
    // far culling
    if (tmin > r.tmax)
    {
        return false;
    }

    if (tzmax < tmax)
    {
        tmax = tzmax;
    }
    // near culling
    if (tmax < r.tmin)
    {
        return false;
    }

    return true;
#endif
}

static bool rayAABBIntersect(const Ray&    r,
                             const float3& bounds_min,
                             const float3& bounds_max,
                             float&        t_near)
{
    float tmin, tmax, tymin, tymax, tzmin, tzmax;

    const float3 bounds[2] = {bounds_min, bounds_max};

    tmin = (bounds[r.dir_sign.x].x - r.orig.x) * r.inv_dir.x;
    tmax = (bounds[1 - r.dir_sign.x].x - r.orig.x) * r.inv_dir.x;
    tymin = (bounds[r.dir_sign.y].y - r.orig.y) * r.inv_dir.y;
    tymax = (bounds[1 - r.dir_sign.y].y - r.orig.y) * r.inv_dir.y;

    if ((tmin > tymax) || (tymin > tmax))
    {
        t_near = NUM_INFINITY;
        return false;
    }
    if (tymin > tmin)
    {
        tmin = tymin;
    }
    if (tymax < tmax)
    {
        tmax = tymax;
    }

    tzmin = (bounds[r.dir_sign.z].z - r.orig.z) * r.inv_dir.z;
    tzmax = (bounds[1 - r.dir_sign.z].z - r.orig.z) * r.inv_dir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
    {
        t_near = NUM_INFINITY;
        return false;
    }
    if (tzmin > tmin)
    {
        tmin = tzmin;
    }
    if (tzmax < tmax)
    {
        tmax = tzmax;
    }

    // in order to handle cases where
    // ray origin is inside of the box
    float rayTMin = r.tmin;
    float rayTMax = r.tmax;
    if (tmin > rayTMin && tmin <= rayTMax)
    {
        t_near = tmin;
        return true;
    }
    else if (tmax > rayTMin && tmax <= rayTMax)
    {
        t_near = tmax;
        return true;
    }
    else
    {
        t_near = NUM_INFINITY;
        return false;
    }
}

AABB::AABB() : m_min(NUM_INFINITY), m_max(-NUM_INFINITY)
{
}

AABB::AABB(const AABB& aabb) : m_min(aabb.m_min), m_max(aabb.m_max)
{
}

AABB::AABB(const float3& min, const float3& max) : m_min(min), m_max(max)
{
}

AABB& AABB::operator=(const AABB& aabb)
{
    m_min = aabb.m_min;
    m_max = aabb.m_max;
    return *this;
}

void AABB::reset()
{
    m_min = float3(NUM_INFINITY);
    m_max = float3(-NUM_INFINITY);
}

void AABB::expandBy(const AABB& box)
{
#if 0
    m_min = f_min(m_min, box.m_min);
    m_max = f_max(m_max, box.m_max);
#else
    // this AABB could be invalid, but the floating point
    // arithmetic would handle it
    assert(b.isValid());

    m_min.x = f_min(m_min.x, box.m_min.x);
    m_min.y = f_min(m_min.y, box.m_min.y);
    m_min.z = f_min(m_min.z, box.m_min.z);

    m_max.x = f_max(m_max.x, box.m_max.x);
    m_max.y = f_max(m_max.y, box.m_max.y);
    m_max.z = f_max(m_max.z, box.m_max.z);
#endif
}

void AABB::expandBy(const float3& point)
{
#if 0
    m_min = f_min(m_min, p);
    m_max = f_max(m_max, p);
#else
    // this AABB could be invalid, but the floating point
    // arithmetic would handle it

    m_min.x = f_min(m_min.x, point.x);
    m_min.y = f_min(m_min.y, point.y);
    m_min.z = f_min(m_min.z, point.z);

    m_max.x = f_max(m_max.x, point.x);
    m_max.y = f_max(m_max.y, point.y);
    m_max.z = f_max(m_max.z, point.z);
#endif
}

float3 AABB::diagonal() const
{
    return m_max - m_min;
}

float AABB::surfaceArea() const
{
    float3 edge = m_max - m_min;
    return fabs(edge.x * edge.y + edge.y * edge.z + edge.x * edge.z) * 2.0f;
}

float AABB::surfaceAreaCost() const
{
    float3 edge = m_max - m_min;
    float  cost = edge.x * edge.y + edge.y * edge.z + edge.x * edge.z;
    assert(cost >= 0.0f);
    return cost;
}

bool AABB::isValid() const
{
    return m_max.x >= m_min.x && m_max.y >= m_min.y && m_max.z >= m_min.z;
}

int AABB::largestDimension() const
{
    assert(isValid());
    float e0 = m_max.x - m_min.x;
    float e1 = m_max.y - m_min.y;
    float e2 = m_max.z - m_min.z;

    int   max_dim = 0;
    float max_e = e0;

    if (e1 > max_e)
    {
        max_e = e1;
        max_dim = 1;
    }

    if (e2 > max_e)
    {
        max_e = e2;
        max_dim = 2;
    }

    return max_dim;
}

bool AABB::intersect(Ray& ray, float& tnear) const
{
    assert(this->isValid());

#if 0
	Logger::debug() << "    ray = " << ray.toString() << endl;
	Logger::debug() << "    AABB = " << this->toString() << endl;
#endif

    bool hit = rayAABBIntersect(ray, this->min(), this->max(), tnear);

    if (hit && tnear > ray.tmin && tnear < ray.tmax)
    {
        return true;
    }
    else
    {
        tnear = NUM_INFINITY;
        return false;
    }
}

AABB AABB::intersect(const AABB& box) const
{
    assert(box.isValid());

    AABB ret;

    ret.m_min.x = f_max(m_min.x, box.m_min.x);
    ret.m_max.x = f_min(m_max.x, box.m_max.x);
    if (ret.m_min.x > ret.m_max.x)
    {
        ret.reset();
        return ret;
    }

    ret.m_min.y = f_max(m_min.y, box.m_min.y);
    ret.m_max.y = f_min(m_max.y, box.m_max.y);
    if (ret.m_min.y > ret.m_max.y)
    {
        ret.reset();
        return ret;
    }

    ret.m_min.z = f_max(m_min.z, box.m_min.z);
    ret.m_max.z = f_min(m_max.z, box.m_max.z);
    if (ret.m_min.z > ret.m_max.z)
    {
        ret.reset();
        return ret;
    }

    return ret;
}

bool AABB::overlap(const Ray& ray) const
{
    assert(this->isValid());

#if 0
    Logger::debug() << "    ray = " << ray.toString() << endl;
    Logger::debug() << "    AABB = " << this->toString() << endl;
#endif

    bool isOverlapping = rayAABBOverlap(ray, m_min, m_max);
#if 0
    if (isOverlapping)
    {
        Logger::debug() << "    ray and AABB do overlap" << endl;
    }
    else
    {
        Logger::debug() << "    ray and AABB do NOT overlap" << endl;
    }
#endif
    return isOverlapping;
}

void AABB::intersectBy(const AABB& box)
{
    assert(box.isValid());

    m_min.x = f_max(m_min.x, box.m_min.x);
    m_max.x = f_min(m_max.x, box.m_max.x);
    if (m_min.x > m_max.x)
    {
        this->reset();
        return;
    }

    m_min.y = f_max(m_min.y, box.m_min.y);
    m_max.y = f_min(m_max.y, box.m_max.y);
    if (m_min.y > m_max.y)
    {
        this->reset();
        return;
    }

    m_min.z = f_max(m_min.z, box.m_min.z);
    m_max.z = f_min(m_max.z, box.m_max.z);
    if (m_min.z > m_max.z)
    {
        this->reset();
        return;
    }
}

float3 AABB::center() const
{
    return (m_min + m_max) * 0.5f;
}

const float3& AABB::min() const
{
    return m_min;
}

const float3& AABB::max() const
{
    return m_max;
}

AABB AABB::clipAxis(float min, float max, int axis) const
{
#if 0
    assert(min < max);
    float3 min_v = m_min;
    float3 max_v = m_max;
    switch (axis)
    {
    case 0:
        min_v.x = clampf(min_v.x, min, max);
        max_v.x = clampf(max_v.x, min, max);
        break;
    case 1:
        min_v.y = clampf(min_v.y, min, max);
        max_v.y = clampf(max_v.y, min, max);
        break;
    case 2:
        min_v.z = clampf(min_v.z, min, max);
        max_v.z = clampf(max_v.z, min, max);
        break;
    default:
        Logger::info() << "AABB invalid axis being clipped against"
            << std::endl;
        exit(-1);
        break;
    }
    return AABB(min_v, max_v);
#else
    assert(min < max);
    assert(axis >= 0 && axis <= 2);
    float3 min_v = m_min;
    float3 max_v = m_max;
    min_v[axis] = clampf(min_v[axis], min, max);
    max_v[axis] = clampf(max_v[axis], min, max);
    return AABB(min_v, max_v);
#endif
}

AABB AABB::clipX(float min, float max) const
{
    assert(min < max);
    float3 min_v = m_min;
    float3 max_v = m_max;
    min_v.x = clampf(min_v.x, min, max);
    max_v.x = clampf(max_v.x, min, max);
    return AABB(min_v, max_v);
}

AABB AABB::clipY(float min, float max) const
{
    assert(min < max);
    float3 min_v = m_min;
    float3 max_v = m_max;
    min_v.y = clampf(min_v.y, min, max);
    max_v.y = clampf(max_v.y, min, max);
    return AABB(min_v, max_v);
}

AABB AABB::clipZ(float min, float max) const
{
    assert(min < max);
    float3 min_v = m_min;
    float3 max_v = m_max;
    min_v.z = clampf(min_v.z, min, max);
    max_v.z = clampf(max_v.z, min, max);
    return AABB(min_v, max_v);
}

std::string AABB::toString() const
{
    std::stringstream ss;
    ss << "min(" << m_min.x << ", " << m_min.y << ", " << m_min.z << "), ";
    ss << "max(" << m_max.x << ", " << m_max.y << ", " << m_max.z << ")";
    return ss.str();
}

AABB expandBy(const AABB& a, const AABB& b)
{
    return AABB(f_min(a.min(), b.min()), f_max(a.max(), b.max()));
}

AABB expandBy(const AABB& a, const float3& point)
{
    return AABB(f_min(a.min(), point), f_max(a.max(), point));
}

AABB getAABB(const TriangleObject& object)
{
    return AABB(f_min(f_min(object.a, object.b), object.c),
                f_max(f_max(object.a, object.b), object.c));
}
