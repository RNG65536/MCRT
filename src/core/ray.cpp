#include "ray.h"

#include <cassert>
#include <sstream>

#include "constants.h"

float Ray::s_ray_epsilon = NUM_EPS_RAY;

void Ray::set_ray_epsilon(float eps)
{
    s_ray_epsilon = eps;
}

Ray::Ray() : tmin(s_ray_epsilon), tmax(NUM_INFINITY)
{
}

Ray::Ray(const Ray& ray) : Ray(ray.orig, ray.dir, ray.tmin, ray.tmax)
{
}

Ray::Ray(const float3& o_, const float3& d_)
    : Ray(o_, d_, s_ray_epsilon, NUM_INFINITY)
{
}

Ray::Ray(const float3& o_, const float3& d_, float tmin_, float tmax_)
    : orig(o_), dir(d_), tmin(tmin_), tmax(tmax_)

{
#ifndef NDEBUG
    float length = d_.length();
    assert(dir.isNormalized());
#else
    // direction must be normalized
#endif

    float inv_x = dir.x == 0.0f ? 0.0f : 1.0f / dir.x;
    float inv_y = dir.y == 0.0f ? 0.0f : 1.0f / dir.y;
    float inv_z = dir.z == 0.0f ? 0.0f : 1.0f / dir.z;
    inv_dir = float3(inv_x, inv_y, inv_z);

    dir_sign.x = dir.x < 0.0f;
    dir_sign.y = dir.y < 0.0f;
    dir_sign.z = dir.z < 0.0f;
}

float3 Ray::at(float t) const
{
    return orig + dir * t;
}

Ray Ray::normalize() const
{
    return Ray(orig, dir.normalize());
}

std::string Ray::toString() const
{
    std::stringstream ss;
    ss << "o(" << orig.x << ", " << orig.y << ", " << orig.z << "), ";
    ss << "d(" << dir.x << ", " << dir.y << ", " << dir.z << ")";
    ss << "t[" << tmin << ", " << tmax << "]";
    return ss.str();
}
