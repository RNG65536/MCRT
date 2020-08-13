#pragma once

#include <memory>
#include "aabb.h"
#include "triangle.h"

// reference to a triangle that is bounded by a clipping aabb
class BoundedTriangle
{
public:
    BoundedTriangle(const std::shared_ptr<Triangle>& triangle, const AABB& aabb);

    BoundedTriangle(const std::shared_ptr<BoundedTriangle>& triangle,
                    const AABB&                             aabb);

    const AABB& clippedBoundingBox() const;

    const std::shared_ptr<Triangle>& originalTriangle() const;

    bool intersect(Ray& ray, HitInfo& hit_info) const;

    bool intersectBound(Ray& ray) const;

    // clip the triangle in the min-max slab and return an even tighter bound
    AABB clipBoundTight(float min, float max, int axis) const;

private:
    const std::shared_ptr<Triangle> m_triangle;
    AABB                            m_clipped_bound;  // clipped bound
};
