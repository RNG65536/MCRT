#pragma once

#include <string>
#include "vectors.h"

// ray

class Ray
{
public:
    Ray();
    Ray(const float3& o_, const float3& d_);

    // normalize direction
    Ray normalize() const;

    // point t units away from the ray origin
    vec3 at(float t) const;

    static void set_ray_epsilon(float eps);

    std::string toString() const;

public:
    float3 orig;      // origin
    float  tmin;      // minimum valid t
    float3 dir;       // direction
    float  tmax;      // maximum valid t
    float3 inv_dir;   // direction reciprocal
    byte3  dir_sign;  // direction sign

    static float s_ray_epsilon;
};
