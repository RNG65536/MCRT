#include "bsdf.h"
#include "constants.h"
#include "geometry.h"
#include "material.h"
#include "numeric.h"
#include "sample.h"

static int count;
inline float debug_f_max(float a, float b)
{
//     assert(a > 0);
//     if (a < -0.5)
//     {
//         //cout << a << endl;
//         count++;
//         if (0 == count % 10000)
//         {
//             cout << " ( " << count << " ) " << endl;
//         }
//     }
//     return fabs(a); // no fireflies
    return f_max(a, b); // fireflies without fudge
}

//////////////////////////////////////////////////////////////////////////

AbstractBSDF::~AbstractBSDF()
{

}

float AbstractBSDF::getContinuationProbability(const vec3& wi, const vec3& nl) const
{
    // ignore RR for the moment
    return 0.0f;
}

void AbstractBSDF::setParent(const Material *mat)
{
    parent = mat;
}

bool LambertianBSDF::isDelta() const
{
    return false;
}

//////////////////////////////////////////////////////////////////////////

vec3 LambertianBSDF::sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW, float *cos_wo, vec3 *bsdf_weight) const
{
    vec3 c = sampleHemisphereCosine(r1, r2);
    OrthonormalFrame onb(nl);
    vec3 wo = onb.toWorld(c);

    //     bsdf_weight = vec3(1.0 / M_PI) * m_color;
    float _cos = debug_f_max(dot(wo, nl), NUM_EPS_COSINE);
    if (forward_pdfW)
    {
        *forward_pdfW = _cos / NUM_PI; // wrt solid angle, i.e., dw is the integrand
        // *forward_pdfW = 1.0 / NUM_PI; // wrt projected solid angle, i.e. cos*dw is the integrand
    }
    if (cos_wo)
    {
        *cos_wo = _cos;
    }
    if (bsdf_weight)
    {
        *bsdf_weight = parent->color() * (1.0f / NUM_PI);
    }
    total_weight = parent->color();
    return wo;
}

float LambertianBSDF::calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const
{
    return debug_f_max(dot(wo, nl), NUM_EPS_COSINE) / NUM_PI; // not the same if wi and wo are swapped ??
}

vec3 LambertianBSDF::evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo, float *forward_pdfW /*= nullptr*/) const
{
    //         if (dot(wi, nl) * dot(wo, nl) < 0)
    //         {
    //             return vec3(0.0f);
    //         }

    float _cos = debug_f_max(dot(wo, nl), NUM_EPS_COSINE);
    if (cos_wo)
    {
        *cos_wo = _cos;
    }
    if (forward_pdfW)
    {
        *forward_pdfW = _cos / NUM_PI;
    }
    vec3 brdf = vec3(1.0 / NUM_PI);

    return brdf * parent->color();
}

LambertianBSDF * LambertianBSDF::clone() const
{
    return new LambertianBSDF(*this);
}

LambertianBSDF::~LambertianBSDF()
{

}

//////////////////////////////////////////////////////////////////////////

bool SpecularBSDF::isDelta() const
{
    return true;
}

vec3 SpecularBSDF::sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW, float *cos_wo, vec3 *bsdf_weight) const
{
    vec3 wo = reflectWorld(wi, nl);

    float _cos = debug_f_max(dot(wo, nl), NUM_EPS_COSINE);
    //     bsdf_weight = m_color / cos_wo;
    if (cos_wo)
    {
        *cos_wo = _cos;
    }
    if (forward_pdfW)
    {
        *forward_pdfW = 1; // wrt projected solid angle, otherwise needs to divide cos_wo
    }
    if (bsdf_weight)
    {
        *bsdf_weight = parent->color();
    }
    total_weight = parent->color();
    return wo;
}

float SpecularBSDF::calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const
{
    return 0;
}

vec3 SpecularBSDF::evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo, float *forward_pdfW /*= nullptr*/) const
{
    //         if (dot(wi, nl) * dot(wo, nl) < 0)
    //         {
    //             return vec3(0.0f);
    //         }

    if (cos_wo)
    {
        *cos_wo = debug_f_max(dot(wo, nl), NUM_EPS_COSINE);
    }
    if (forward_pdfW)
    {
        *forward_pdfW = 0;
    }
    return vec3(0.0f);
}

SpecularBSDF * SpecularBSDF::clone() const
{
    return new SpecularBSDF(*this);
}

SpecularBSDF::~SpecularBSDF()
{

}

//////////////////////////////////////////////////////////////////////////

bool GlossyBSDF::isDelta() const
{
    return false;
}

vec3 GlossyBSDF::sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW, float *cos_wo, vec3 *bsdf_weight) const
{
    float g = parent->glossiness();
    vec3 c = sampleHemisphereCosinePower(parent->glossiness(), r1, r2);
    vec3 refl = reflectWorld(wi, nl);
    OrthonormalFrame onb(refl);
    vec3 wo = onb.toWorld(c);
    if (dot(wo, nl) < 0)
    {
        assert(false);
        wo = vec3(0.0f);
    }

    float _cos = debug_f_max(dot(wo, nl), NUM_EPS_COSINE);
    if (cos_wo)
    {
        *cos_wo = _cos;
    }
    if (forward_pdfW)
    {
        // wrt solid angle
        *forward_pdfW = (g + 1.0) / (2.0 * NUM_PI) * pow(debug_f_max(dot(refl, wo), NUM_EPS_COSINE), g);
        // wrt projected solid angle
//         *forward_pdfW = (g + 1.0) / (2.0 * NUM_PI) * pow(f_max(dot(refl, wo), 0.0), g) / _cos;
    }

    float wo_n = dot(wo, nl);
    float wi_n = dot(wi, nl);
    float weight = (g + 2.0) / (g + 1.0) * _cos
        / debug_f_max(fabs(wi_n), fabs(wo_n));
    if (weight != weight) // handle NaN
    {
        weight = 0;
    }

    if (bsdf_weight)
    {
        *bsdf_weight = ((g + 2.0) / (2.0 * NUM_PI) * pow(debug_f_max(dot(refl, wo), NUM_EPS_COSINE), g)
            / debug_f_max(fabs(wi_n), fabs(wo_n))) * parent->color();
    }

    total_weight = weight * parent->color(); // bsdf_weight * cos_wo / forward_pdfW
    return wo;
}

float GlossyBSDF::calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const
{
    float g = parent->glossiness();
    vec3 refl = reflectWorld(wi, nl);
    return (g + 1.0) / (2.0 * NUM_PI) * pow(debug_f_max(dot(refl, wo), NUM_EPS_COSINE), g);
}

vec3 GlossyBSDF::evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo, float *forward_pdfW /*= nullptr*/) const
{
    float wo_n = dot(wo, nl);
    float wi_n = dot(wi, nl);
    vec3 refl = reflectWorld(wi, nl);
    float g = parent->glossiness();

    if (cos_wo)
    {
        *cos_wo = debug_f_max(wo_n, NUM_EPS_COSINE);
    }
    if (forward_pdfW)
    {
        *forward_pdfW = (g + 1.0) / (2.0 * NUM_PI) * pow(debug_f_max(dot(refl, wo), NUM_EPS_COSINE), g);
    }

    return float((g + 2.0) / (2.0 * NUM_PI) * pow(debug_f_max(dot(refl, wo), NUM_EPS_COSINE), g)
        / debug_f_max(fabs(wi_n), fabs(wo_n))) * parent->color();
}

GlossyBSDF * GlossyBSDF::clone() const
{
    return new GlossyBSDF(*this);
}

GlossyBSDF::~GlossyBSDF()
{

}

//////////////////////////////////////////////////////////////////////////

bool RefractiveBSDF::isDelta() const
{
    return true;
}

vec3 RefractiveBSDF::sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW, float *cos_wo, vec3 *bsdf_weight) const
{
    bool outside = dot(wi, nl) > 0; // Ray from outside going in?
    vec3 n = outside ? nl : -nl;

    vec3 refl_dir = normalize(reflectWorld(wi, n));     // Ideal dielectric REFRACTION
    float nc = 1, nt = parent->ior();

    // nnt is reciprocal relative index of refraction, ddn is negated cosine of incidence angle
    float nnt = outside ? nc / nt : nt / nc;
    float ddn = dot(-wi, n);
    float cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
    if (cos2t < 0)
    {
        // Total internal reflection
        total_weight = parent->color();
        if (cos_wo)
        {
            *cos_wo = debug_f_max(dot(refl_dir, n), NUM_EPS_COSINE);
        }
        if (forward_pdfW)
        {
            *forward_pdfW = 1.0; // wrt projected solid angle
        }
        if (bsdf_weight)
        {
            *bsdf_weight = parent->color();
        }
        return refl_dir;
    }
    else
    {
        vec3 tdir = normalize(-wi * nnt - n * (ddn * nnt + std::sqrt(cos2t)));
#if 0
        float a = nt - nc;
        float b = nt + nc;
        float R0 = a * a / (b * b);
        float c = 1 - (outside ? -ddn : dot(tdir, nl)); // R0 is the reflectance at normal incidence
        float Re = R0 + (1 - R0) * c * c * c * c * c; // Schlick's approximation for Fresnel reflectance
#else
        float Re = fresnelReflection(dot(wi, nl), nt / nc); // Schlick's approximation for Fresnel reflectance
#endif
        float Tr = 1 - Re;
        // float P = 0.25 + 0.5 * Re; // heuristic
        //float P = 0.5; // also works but not as well
        float P = Re; // also works but not as well

        total_weight = vec3(1.0f);

        float RP = Re / P;
        float TP = Tr / (1 - P); // P is used for importance sampling
        if (r1 < P)
        {
            total_weight *= RP * parent->color();
            if (cos_wo)
            {
                *cos_wo = debug_f_max(dot(refl_dir, n), NUM_EPS_COSINE);
            }
            if (forward_pdfW)
            {
                *forward_pdfW = P; // wrt projected solid angle
            }
            if (bsdf_weight)
            {
                *bsdf_weight = Re * parent->color();
            }
            return refl_dir;
        }
        else
        {
            total_weight *= TP * parent->color();
            if (cos_wo)
            {
                *cos_wo = debug_f_max(dot(tdir, -n), NUM_EPS_COSINE);
            }
            if (forward_pdfW)
            {
                *forward_pdfW = 1 - P; // wrt projected solid angle
            }
            if (bsdf_weight)
            {
                *bsdf_weight = Tr * parent->color();
            }
            return tdir;
        }
    }
}

float RefractiveBSDF::calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const
{
    return 0;
}

vec3 RefractiveBSDF::evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo, float *forward_pdfW /*= nullptr*/) const
{
    if (cos_wo)
    {
        *cos_wo = debug_f_max(dot(wo, nl), NUM_EPS_COSINE);
    }
    if (forward_pdfW)
    {
        *forward_pdfW = 0;
    }
    return vec3(0.0f);
}

RefractiveBSDF * RefractiveBSDF::clone() const
{
    return new RefractiveBSDF(*this);
}

RefractiveBSDF::~RefractiveBSDF()
{

}
