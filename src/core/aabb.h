#ifndef aabb_h__
#define aabb_h__

#include "vectors.h"

// axis aligned bounding box

class TriangleObject;
class Ray;

class AABB
{
public:
    vec3 _min, _max;
    AABB();
    AABB(const vec3 &min_, const vec3 &max_);
    AABB(const AABB &aabb);
    AABB &operator=(const AABB &aabb);
    void invalidify();
    vec3 getCenter() const;
    bool isValid() const;
    bool isEnclosing(const vec3& p) const;
    bool unbounded() const;
    bool intersect(const Ray& ray, const vec3& inverseDirection, float closestKnownT, float& nearT, float& farT) const;
    bool intersect(const Ray& ray, const vec3& inverseDirection, float closestKnownT) const;
    vec3 getNormal(const vec3& p0) const;
    size_t largestDimension() const;
    size_t smallestDimension() const;
    bool nearEdge(vec3& a, vec3& b, const vec3& p);
    bool onFrame(const vec3& p) const;
    void enclose(const AABB& b);
    void enclose(const vec3& p);
    float surface_area() const;
    float diagonalLength() const;
};
AABB enclose(const AABB& a, const AABB& b);
AABB enclose(const AABB& a, const vec3& point);
AABB getAABB(const TriangleObject& object);

#endif // aabb_h__
