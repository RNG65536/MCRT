#ifndef volumetric_h__
#define volumetric_h__

#include "vectors.h"

class Ray;

float sampleSegment(float epsilon, float sigma);
vec3 sampleHG(float g /* mean cosine */, float e1, float e2);
float scatter(const Ray &r, Ray &sRay, float sigma_s, float &s, float e0, float e1, float e2);

struct VolumetricProps {
    VolumetricProps();

    float scale;
    vec3 sigma_s, sigma_a;

    vec3 final_sigma_t;//extinction coefficient
    vec3 final_sigma_s;
    vec3 final_sigma_a;
    vec3 albedo;//single scattering albedo

    int maxDepths[3];
};

// #define DEFAULT_IOR 1.49

// int maxDepths[3] = {102, 28, 16};
// int maxDepths[3] = {195, 54, 30};

extern VolumetricProps default_medium;
extern const int SSS_DEPTH;

#endif // volumetric_h__
