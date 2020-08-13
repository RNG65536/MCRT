#include <cfloat>
#include <cmath>
#include "constants.h"
#include "geometry.h"
#include "numeric.h"
#include "ray.h"
#include "sample.h"
#include "volumetric.h"

const int       SSS_DEPTH = 100;
VolumetricProps default_medium;

float sampleSegment(float epsilon, float sigma)
{
    return -logf(f_max(epsilon, FLT_MIN)) / sigma;
}

vec3 sampleHG(float g /* mean cosine */, float e1, float e2)
{
    // float s=2.0*e1-1.0, f = (1.0-g*g)/(1.0+g*s), cost =
    // 0.5*(1.0/g)*(1.0+g*g-f*f), sint = std::sqrtf(1.0-cost*cost);
    float s = 1.0f - 2.0f * e1, denom = 1.0f + g * s;
    float cost = (s + 2.0f * g * g * g * (-1.0f + e1) * e1 + g * g * s +
                  2.0f * g * (1.0f - e1 + e1 * e1)) /
                 (denom * denom),
          sint = std::sqrt(1.0f - cost * cost);
    return vec3(
        cosf(2.0f * NUM_PI * e2) * sint, sinf(2.0f * NUM_PI * e2) * sint, cost);
}

float scatter(const Ray& r,
              Ray&       sRay,
              float      sigma_s,
              float&     s,
              float      e0,
              float      e1,
              float      e2)
{
    s = sampleSegment(e0, sigma_s);
    vec3 x = r.orig + r.dir * s;
    //     Vector3 dir = sampleHG(0.9f, e1, e2); //Sample a direction ~
    //     Henyey-Greenstein's phase function
    vec3 dir = sampleSphereUniform(
        e1, e2);  // Sample a direction ~ uniform phase function

    //     vec3 u, v;
    //     generateOrthoBasis(u, v, r.d);
    //     dir = u * dir.x + v * dir.y + r.d * dir.z;
    OrthonormalFrame onb(r.dir);
    dir = onb.toWorld(dir);

    //     dir.norm(); //by construction of the sample coordinates, dir is
    //     guaranteed to be unit
    sRay = Ray(x, dir);
    return 1.0f;
}

VolumetricProps::VolumetricProps()
    : sigma_s(vec3(1.09f, 1.59f, 1.79f)),
      sigma_a(vec3(0.013f, 0.070f, 0.145f)),
      //     scale(20)
      //     scale(10)
      //     scale(5)
      //     scale(2)
      scale(30)
{
    final_sigma_t = (sigma_s + sigma_a) * scale;
    final_sigma_s = (sigma_s)*scale;
    final_sigma_a = (sigma_a)*scale;
    albedo = div(sigma_s, sigma_s + sigma_a);

    printf("sigma_a:<%f, %f, %f>\n",
           final_sigma_a.x,
           final_sigma_a.y,
           final_sigma_a.z);
    printf("sigma_s:<%f, %f, %f>\n",
           final_sigma_s.x,
           final_sigma_s.y,
           final_sigma_s.z);
    printf("sigma_t:<%f, %f, %f>\n",
           final_sigma_t.x,
           final_sigma_t.y,
           final_sigma_t.z);
    printf("albedo<%f, %f, %f>\n", albedo.x, albedo.y, albedo.z);
    float w = 0.03f;
    //         real w=0.1;
    //         real w=0.3;
    //         real w=0.5;
    printf("specific to %f%%, depth should be <%d, %d, %d>\n",
           w * 100,
           int(logf(w) / logf(albedo.x)) + 1,
           int(logf(w) / logf(albedo.y)) + 1,
           int(logf(w) / logf(albedo.z)) + 1);

    maxDepths[0] = int(logf(w) / logf(albedo.x)) + 1;
    maxDepths[1] = int(logf(w) / logf(albedo.y)) + 1;
    maxDepths[2] = int(logf(w) / logf(albedo.z)) + 1;
}
