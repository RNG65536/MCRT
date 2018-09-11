#ifndef bsdf_h__
#define bsdf_h__

#include "vectors.h"

class Material;

class AbstractBSDF
{
protected:
    const Material *parent = nullptr;

public:
    AbstractBSDF() = default;
    virtual ~AbstractBSDF();
    virtual AbstractBSDF* clone() const = 0;
    void setParent(const Material *mat);

    // Evaluate, which takes input and output directions of light and a normal, and returns the BSDF weight, cosine of the angle of the input (??) direction, and color absorption of the scattering event. Evaluate can also optionally return the probability of the output direction given the input direction, with respect to solid angle.
    virtual vec3 evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo = nullptr, float *forward_pdfW = nullptr) const = 0;

    // CalculatePDFW, which takes the input and output directions of light and a normal, and returns the forward probability of the output direction given the input direction. In order to make the BSDF operate bidirectionally, this function also needs to be able to return the backwards probability if the input and output are reversed.
    virtual float calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const = 0;

    // Sample, which takes in an input direction, a normal, and a random number generator and returns an output direction, the BSDF weight, the forward probability of the output direction, and the cosine of the input (??) angle.
    virtual vec3 sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW = nullptr, float *cos_wo = nullptr, vec3 *bsdf_weight = nullptr) const = 0;

    // IsDelta, which returns true if the BSDF¡¯s probability distribution function is a Dirac delta function and false otherwise. This attribute is important for allowing BDPT and VCM to handle perfectly specular BSDFs correctly, since perfectly specular BSDFs are something of a special case.
    virtual bool isDelta() const = 0;

    // GetContinuationProbability, which takes in an input direction and normal and returns the probability of ending a ray path at this BSDF.This function is used for Russian Roulette early path termination.
    virtual float getContinuationProbability(const vec3& wi, const vec3& nl) const;
};

class LambertianBSDF : public AbstractBSDF
{
public:
    LambertianBSDF() = default;
    ~LambertianBSDF();
    LambertianBSDF *clone() const;

    vec3 evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo = nullptr, float *forward_pdfW = nullptr) const;
    float calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const;
    vec3 sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW = nullptr, float *cos_wo = nullptr, vec3 *bsdf_weight = nullptr) const;
    bool isDelta() const;
};

class SpecularBSDF : public AbstractBSDF
{
public:
    SpecularBSDF() = default;
    ~SpecularBSDF();
    SpecularBSDF *clone() const;
    vec3 evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo = nullptr, float *forward_pdfW = nullptr) const;
    float calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const;
    vec3 sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW = nullptr, float *cos_wo = nullptr, vec3 *bsdf_weight = nullptr) const;
    bool isDelta() const;
};

class GlossyBSDF : public AbstractBSDF
{
public:
    GlossyBSDF() = default;
    ~GlossyBSDF();
    GlossyBSDF * clone() const;

    vec3 evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo = nullptr, float *forward_pdfW = nullptr) const;
    float calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const;
    vec3 sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW = nullptr, float *cos_wo = nullptr, vec3 *bsdf_weight = nullptr) const;
    bool isDelta() const;
};

class RefractiveBSDF : public AbstractBSDF
{
public:
    RefractiveBSDF() = default;
    ~RefractiveBSDF();
    RefractiveBSDF *clone() const;

    vec3 evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo = nullptr, float *forward_pdfW = nullptr) const;
    float calculatePdfW(const vec3& wi, const vec3& nl, const vec3& wo) const;
    vec3 sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW = nullptr, float *cos_wo = nullptr, vec3 *bsdf_weight = nullptr) const;
    bool isDelta() const;
};

#endif // bsdf_h__