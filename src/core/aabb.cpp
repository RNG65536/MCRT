#include <cmath>
#include "aabb.h"
#include "constants.h"
#include "numeric.h"
#include "ray.h"
#include "triangle.h"

AABB enclose(const AABB& a, const AABB& b) { return AABB(f_min(a._min, b._min), f_max(a._max, b._max)); }
AABB enclose(const AABB& a, const vec3& point) { return AABB(f_min(a._min, point), f_max(a._max, point)); }

AABB getAABB(const TriangleObject& object)
{
    return AABB(
        f_min(f_min(object.a, object.b), object.c),
        f_max(f_max(object.a, object.b), object.c)
        );
}


float AABB::surface_area() const
{
    vec3 d = _max - _min;
    return 2.0f * (d.x*d.y + d.y*d.z + d.x*d.z);
}

void AABB::enclose(const vec3& p)
{
    _min = f_min(_min, p); _max = f_max(_max, p);
}

void AABB::enclose(const AABB& b)
{
    _min = f_min(_min, b._min); _max = f_max(_max, b._max);
}

bool AABB::onFrame(const vec3& p) const
{
    float u = (p.x - _min.x) / (_max.x - _min.x);
    float v = (p.y - _min.y) / (_max.y - _min.y);
    float w = (p.z - _min.z) / (_max.z - _min.z);
    return
        u<0.01 && v<0.01 ||
        u<0.01 && v>0.99 ||
        u>0.99 && v<0.01 ||
        u>0.99 && v>0.99 ||
        u<0.01 && w<0.01 ||
        u<0.01 && w>0.99 ||
        u>0.99 && w<0.01 ||
        u>0.99 && w>0.99 ||
        v < 0.01 && w<0.01 ||
        v<0.01 && w>0.99 ||
        v>0.99 && w<0.01 ||
        v>0.99 && w > 0.99;
}

bool AABB::nearEdge(vec3& a, vec3& b, const vec3& p)
{
    vec3 e = b - a;
    float u = dot(e, p - a) / e.lengthSquared();
    vec3 c = a + e*u;
    return (p - u).lengthSquared() < 0.05*0.05;
}

size_t AABB::smallestDimension() const
{
    float dx = fabs(_max.x - _min.x);
    float dy = fabs(_max.y - _min.y);
    float dz = fabs(_max.z - _min.z);
    if (dx < dy && dx < dz) return 0;
    if (dy < dz)  return 1; return 2;
}

size_t AABB::largestDimension() const
{
    float dx = fabs(_max.x - _min.x);
    float dy = fabs(_max.y - _min.y);
    float dz = fabs(_max.z - _min.z);
    if (dx > dy && dx > dz) return 0;
    if (dy > dz)  return 1; return 2;
}

vec3 AABB::getNormal(const vec3& p0) const
{
    vec3 c(getCenter());
    if (fabs(p0.x - _min.x) < NUM_EPS || fabs(p0.x - _max.x) < NUM_EPS)
        return p0.x < c.x ? vec3(-1, 0, 0) : vec3(1, 0, 0);
    else if (fabs(p0.y - _min.y) < NUM_EPS || fabs(p0.y - _max.y) < NUM_EPS)
        return p0.y < c.y ? vec3(0, -1, 0) : vec3(0, 1, 0);
    else
        return p0.z < c.z ? vec3(0, 0, -1) : vec3(0, 0, 1);
}

bool AABB::intersect(const Ray& ray, const vec3& inverseDirection, float closestKnownT) const
{
    float nearT, farT;
    return intersect(ray, inverseDirection, closestKnownT, nearT, farT);
}

bool AABB::intersect(const Ray& ray, const vec3& inverseDirection, float closestKnownT, float& nearT, float& farT) const
{
    //this implementation is slightly faster
    bool xDirNegative = ray.d.x < 0;
    bool yDirNegative = ray.d.y < 0;
    bool zDirNegative = ray.d.z < 0;
    float tmin, tmax, tymin, tymax, tzmin, tzmax;
    if (xDirNegative) {
        tmin = (_max.x - ray.o.x) * inverseDirection.x;
        tmax = (_min.x - ray.o.x) * inverseDirection.x;
    }
    else {
        tmin = (_min.x - ray.o.x) * inverseDirection.x;
        tmax = (_max.x - ray.o.x) * inverseDirection.x;
    }
    if (yDirNegative) {
        tymin = (_max.y - ray.o.y) * inverseDirection.y;
        tymax = (_min.y - ray.o.y) * inverseDirection.y;
    }
    else {
        tymin = (_min.y - ray.o.y) * inverseDirection.y;
        tymax = (_max.y - ray.o.y) * inverseDirection.y;
    }
    //         float tmin = ((xDirNegative ? _max.x : _min.x) - ray.o.x) * inverseDirection.x;
    //         float tmax = ((xDirNegative ? _min.x : _max.x) - ray.o.x) * inverseDirection.x;
    //         float tymin = ((yDirNegative ? _max.y : _min.y) - ray.o.y) * inverseDirection.y;
    //         float tymax = ((yDirNegative ? _min.y : _max.y) - ray.o.y) * inverseDirection.y;
    if (tmin > tymax || tymin > tmax) { return false; }
    if (tymin > tmin) { tmin = tymin; }
    if (tymax < tmax) { tmax = tymax; }
    if (zDirNegative) {
        tzmin = (_max.z - ray.o.z) * inverseDirection.z;
        tzmax = (_min.z - ray.o.z) * inverseDirection.z;
    }
    else {
        tzmin = (_min.z - ray.o.z) * inverseDirection.z;
        tzmax = (_max.z - ray.o.z) * inverseDirection.z;
    }
    //         float tzmin = ((zDirNegative ? _max.z : _min.z) - ray.o.z) * inverseDirection.z;
    //         float tzmax = ((zDirNegative ? _min.z : _max.z) - ray.o.z) * inverseDirection.z;
    if (tmin > tzmax || tzmin > tmax) { return false; }
    if (tzmin > tmin) { tmin = tzmin; }
    if (tzmax < tmax) { tmax = tzmax; }
    nearT = tmin;
    farT = tmax;
    return (tmin < closestKnownT) //this accounts for early termination
        && (tmax > NUM_EPS);
}

bool AABB::unbounded() const
{
    return _min.x == -NUM_INFINITY || _min.y == -NUM_INFINITY || _min.z == -NUM_INFINITY || _max.x == NUM_INFINITY || _max.y == NUM_INFINITY || _max.z == NUM_INFINITY;
}

bool AABB::isEnclosing(const vec3& p) const
{
    return p.x > _min.x&&p.x<_max.x&&p.y>_min.y&&p.y<_max.y&&p.z>_min.z&&p.z < _max.z;
}

bool AABB::isValid() const
{
    return _max.x >= _min.x && _max.y >= _min.y && _max.z >= _min.z;
}

vec3 AABB::getCenter() const
{
    return vec3(_max + _min) * 0.5f;
}

void AABB::invalidify()
{
    _min = vec3(NUM_INFINITY, NUM_INFINITY, NUM_INFINITY); _max = vec3(-NUM_INFINITY, -NUM_INFINITY, -NUM_INFINITY);
}

AABB & AABB::operator=(const AABB &aabb)
{
    _min = aabb._min; _max = aabb._max; return *this;
}

AABB::AABB(const AABB &aabb) : _min(aabb._min), _max(aabb._max)
{

}

AABB::AABB(const vec3 &min_, const vec3 &max_) : _min(min_), _max(max_)
{

}

AABB::AABB()
{
    invalidify();
}

float AABB::diagonalLength() const
{
    return length(_max - _min);
}
