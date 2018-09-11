#include "ray.h"

vec3 Ray::at(float t) const
{
    return o + d * t;
}

Ray Ray::norm() const
{
    return Ray(o, d.norm());
}

Ray::Ray(const vec3 &o_, const vec3 &d_) :o(o_), d(d_.norm())
{

}

Ray::Ray()
{

}
