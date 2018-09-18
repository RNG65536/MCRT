#pragma once

#include "vectors.h"

// ray

class Ray
{
public:
    Ray();
    Ray(const vec3& o_, const vec3& d_);

    // normalize direction
    Ray normalize() const;

    // point t units away from the ray origin
    vec3 at(float t) const;

public:
    // origin
    vec3 orig;

    // direction
    vec3 dir;
};
