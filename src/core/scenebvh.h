#pragma once

#include <vector>
#include "intersection.h"
#include "triangle.h"

class Ray;

// scene traversal with acceleration structure

class SceneBVH
{
public:
    SceneBVH(const std::vector<TriangleObject>& triangles);
    ~SceneBVH();

    HitInfo intersect(const Ray& ray) const;
    bool    occluded(const Ray& ray, float tfar) const;

private:
    const std::vector<TriangleObject>& m_triangles_ref;

    // geometry and acceleration info
    void* m_ctx = nullptr;
};
