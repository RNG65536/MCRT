#include <vector>
#include "bsdf.h"
#include "material.h"

Material::Material(const vec3& c, float e, MaterialType t)
    : m_color(c), m_emission(e), m_type(t)
{
    switch (m_type)
    {
        default:
        case DIFF:
            m_bsdf = new LambertianBSDF();
            break;
        case SSSS:
        case SPEC:
            m_bsdf = new SpecularBSDF();
            break;
        case REFR:
            m_bsdf = new RefractiveBSDF();
            break;
        case GLSY:
            m_bsdf = new GlossyBSDF();
            break;
    }

    m_bsdf->setParent(this);
}

Material::Material(const Material& mat)
{
    *this = mat;
}

MaterialType Material::type() const
{
    return m_type;
}

float Material::emission() const
{
    return m_emission;
}

vec3 Material::color() const
{
    return m_color;
}

Material::~Material()
{
    if (m_bsdf)
    {
        delete m_bsdf;
    }
}

Material& Material::operator=(const Material& mat)
{
    this->m_color = mat.m_color;
    this->m_emission = mat.m_emission;
    this->m_ior = mat.m_ior;
    this->m_roughness = mat.m_roughness;
    this->m_type = mat.m_type;
    if (mat.m_bsdf)
    {
        this->m_bsdf = mat.m_bsdf->clone();
        this->m_bsdf->setParent(this);
    }
    return *this;
}

float Material::roughness() const
{
    return m_roughness;
}

float Material::glossiness() const
{
    return m_glossiness;
}

float Material::ior() const
{
    return m_ior;
}

bool Material::isLight() const
{
    return m_type == LGHT || m_emission > 0;
}

float Material::getContinuationProbability(const vec3& wi, const vec3& nl) const
{
    return m_bsdf->getContinuationProbability(wi, nl);
}

bool Material::isDelta() const
{
    return m_bsdf->isDelta();
}

vec3 Material::sample(const vec3& wi,
                      const vec3& nl,
                      float       r1,
                      float       r2,
                      vec3&       total_weight,
                      float*      forward_pdfW,
                      float*      cos_wo,
                      vec3*       bsdf_weight) const
{
    return m_bsdf->sample(
        wi, nl, r1, r2, total_weight, forward_pdfW, cos_wo, bsdf_weight);
}

float Material::pdfW(const vec3& wi, const vec3& nl, const vec3& wo) const
{
    return m_bsdf->calculatePdfW(wi, nl, wo);
}

vec3 Material::evaluate(const vec3& wi,
                        const vec3& nl,
                        const vec3& wo,
                        float*      cos_wo,
                        float*      forward_pdfW) const
{
    return m_bsdf->evaluate(wi, nl, wo, cos_wo, forward_pdfW);
}

// TODO : implement material library
static std::vector<Material> mat_lib;
const Material&              todo_getMaterial(int mat_id)
{
    return mat_lib[mat_id];
}
int todo_addMaterial(const Material& mat)
{
    mat_lib.push_back(mat);
    return static_cast<int>(mat_lib.size()) - 1;
}
