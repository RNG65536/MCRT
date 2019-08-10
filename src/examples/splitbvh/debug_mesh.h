#pragma once

#include <vector>
#include "triangle.h"
#include "aabb.h"

class DebugMesh
{
public:
    DebugMesh();

    std::vector<Triangle>::const_iterator begin() const
    {
        return m_triangles.cbegin();
    }

    std::vector<Triangle>::const_iterator end() const
    {
        return m_triangles.cend();
    }

    const AABB& boundingBox() const
    {
        return m_boundingBox;
    }

    const std::vector<Triangle>& triangles() const
    {
        return m_triangles;
    }

    bool intersect(Ray& ray, HitInfo& hit_info) const;

private:
    void initMesh();
    void initRandomMesh();

    std::vector<Triangle> m_triangles;
    AABB                  m_boundingBox;
};
