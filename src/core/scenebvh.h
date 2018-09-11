#ifndef scenebvh_h__
#define scenebvh_h__

#include <vector>
#include "triangle.h"
#include "intersection.h"

class Ray;

// scene management with bvh

class SceneBVH
{
private:
    void *device = nullptr;
    void *RtcScene = nullptr;

    const std::vector<TriangleObject>& triangleObjectPool;

public:
    SceneBVH(const std::vector<TriangleObject>& triangles);
    ~SceneBVH();

    HitInfo intersect(const Ray& ray) const;
    bool occluded(const Ray& ray, float tfar) const;
};

#endif // scenebvh_h__
