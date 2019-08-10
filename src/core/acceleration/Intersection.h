#pragma once

#include "Ray.h"

namespace debug
{
// Moller-Trumbore for ray-triangle intersection
// and BVH acceleration for ray-mesh intersection
class Intersection
{
public:
    virtual float intersectWith(const Ray& ray) const = 0;
};
}
