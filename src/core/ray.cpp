#include "ray.h"

vec3 Ray::at(float t) const
{
    return orig + dir * t;
}

Ray Ray::normalize() const
{
    return Ray(orig, dir.normalize());
}

Ray::Ray(const vec3& o_, const vec3& d_) : orig(o_), dir(d_.normalize())
{
}

Ray::Ray()
{
}
