#ifndef ray_h__
#define ray_h__

#include "vectors.h"

// ray

class Ray
{
public:
    vec3 o, d;
    Ray();
    Ray(const vec3 &o_, const vec3 &d_);
    Ray norm() const;
    vec3 at(float t) const;
};

#endif // ray_h__
