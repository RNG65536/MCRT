#ifndef bssrdf_h__
#define bssrdf_h__

#include "vectors.h"

class BSSRDF
{
public:
    virtual ~BSSRDF();
    virtual vec3 evaluate(const vec3& xi, const vec3& ni, const vec3& wi,
        const vec3& xo, const vec3& no, const vec3& wo) = 0;
    virtual float evaluate_channel(const vec3& xi, const vec3& ni, const vec3& wi,
        const vec3& xo, const vec3& no, const vec3& wo, int i) = 0;
};

class DirectionalDipoleBSSRDF
{
    vec3 sigma_s, sigma_a;
    float g;
    vec3 eta;
    vec3 Cp_norm, Cp, Ce, A; // Refractive
    vec3 Fdt;

public:
    DirectionalDipoleBSSRDF(const vec3& sigma_s, const vec3& sigma_a, float g, const vec3& eta);
    vec3 evaluate(const vec3& xi, const vec3& ni, const vec3& wi,
                  const vec3& xo, const vec3& no, const vec3& wo);
    float evaluate_channel(const vec3& xi, const vec3& ni, const vec3& wi,
                           const vec3& xo, const vec3& no, const vec3& wo, int i);
};

class StandardDipoleBSSRDF
{
    vec3 sigma_s, sigma_a;
    float g;
    vec3 eta;
    vec3 Cp, A;
    vec3 Fdt;

public:
    StandardDipoleBSSRDF(const vec3& sigma_s, const vec3& sigma_a, float g, const vec3& eta);
    vec3 evaluate(const vec3& xi, const vec3& ni, const vec3& wi,
                  const vec3& xo, const vec3& no, const vec3& wo);
    float evaluate_channel(const vec3& xi, const vec3& ni, const vec3& wi,
                           const vec3& xo, const vec3& no, const vec3& wo, int i);
};

#endif // bssrdf_h__
