#ifndef geometry_h__
#define geometry_h__

#include <string>
#include <glm/glm.hpp>
#include "vectors.h"

class Vert;
class Material;

class OrthonormalFrame
{
    vec3 m_x; // tangent
    vec3 m_y; // bitangent
    vec3 m_z; // normal

public:
    OrthonormalFrame();
    OrthonormalFrame(const vec3& _z);
    void setFromNormal(const vec3& _z);
    vec3 toWorld(const vec3& coords) const;
    vec3 toLocal(const vec3& coords) const;
    const vec3& getTangent() const;
    const vec3& getBitangent() const;
    const vec3& getNormal() const;
};

vec3 reflectWorld(const vec3& wi, const vec3& nl);
vec3 reflectLocal(const vec3& wi);

float cosineFactor(const vec3 wi, const vec3& nl);

float fresnelReflection(float cos_theta_i, float ior);

// local sampling PDFs and standard terms
float geometryTerm(const Vert& e0, const Vert& e1);
float geometryTerm(const vec3& p0, const vec3& n0, const vec3& p1, const vec3& n1);
// pdf conversion, r^2 falloff
float directionToArea(const Vert& current, const Vert& next);

//////////////////////////////////////////////////////////////////////////
// Utilities for converting PDF between Area (A) and Solid angle (W)
// WtoA = PdfW * cosine / distance_squared
// AtoW = PdfA * distance_squared / cosine
// -- cosThere = cosine on the emitter
// -- pdfW = 1 / Omega = r^2 / (Omega * r^2) = pdfA_proj * r^2 = pdfA / cos * r^2
// -- A_proj = A * cos, and pdf is reciprocal to A for the same sample space
float pdfWtoA(const vec3& v0, const vec3& v1, const vec3& n1);
float pdfAtoW(const vec3& v0, const vec3& v1, const vec3& n1);
float pdfWtoA(float pdfW, float dist, float cosineThere);
float pdfAtoW(float pdfA, float dist, float cosineThere);
float __pdfWtoA(float pdfW, float dist_sq, float cosineThere);
float __pdfAtoW(float pdfA, float dist_sq, float cosineThere);

bool isAligned(const vec3& a, const vec3& b);

class Scene;
void loadModel(Scene& scene, std::string filename, const glm::mat4& transform, const Material& mat);

#endif // geometry_h__
