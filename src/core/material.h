#ifndef material_h__
#define material_h__

#include "vectors.h"

// materials (to be deprecated)

typedef enum{ LGHT, DIFF, SPEC, REFR, GLSY, PHNG, VOLM, SSSS } MaterialType;

class AbstractBSDF;

// todo : integrate bsdf into this class (esp. API)
class Material
{
private:
    vec3 m_color;
    float m_emission; // have no use, todo : deprecate this
    float m_ior = 1.33;
    float m_roughness = 1.0;
    float m_glossiness = 25.0;
    MaterialType m_type;
    AbstractBSDF *m_bsdf = nullptr;

public:
    vec3 color() const;
    float emission() const;
    float ior() const;
    float roughness() const;
    float glossiness() const;
    MaterialType type() const;
    Material() = default;
    Material(const Material& mat);
    Material& operator= (const Material& mat);
    ~Material();
    Material(const vec3& c, float e, MaterialType t);
    bool isLight() const;

    vec3 evaluate(const vec3& wi, const vec3& nl, const vec3& wo, float *cos_wo = nullptr, float *forward_pdfW = nullptr) const;
    float pdfW(const vec3& wi, const vec3& nl, const vec3& wo) const;
    vec3 sample(const vec3& wi, const vec3& nl, float r1, float r2, vec3& total_weight, float *forward_pdfW = nullptr, float *cos_wo = nullptr, vec3 *bsdf_weight = nullptr) const;
    bool isDelta() const;
    float getContinuationProbability(const vec3& wi, const vec3& nl) const;
};

int todo_addMaterial(const Material& mat);
const Material& todo_getMaterial(int mat_id);

#endif // material_h__
