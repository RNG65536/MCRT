#include "constants.h"
#include "light.h"
#include "numeric.h"
#include "sample.h"
#include "scene.h"

bool AreaLight::IsDelta() const {
    return false;
}

bool AreaLight::IsFinite() const
{
    return true;
}

vec3 AreaLight::GetRadiance(const SceneSphere &/*aSceneSphere*/, const vec3 &_rayDirection, const vec3 &_hitPoint, float *directPdfA_ /*= nullptr*/, float *emissionPdfW_ /*= nullptr*/) const
{
    const float cosOutL = f_max(0.0f, dot(m_frame.getNormal(), -_rayDirection));

    if (cosOutL == 0)
    {
        return vec3(0);
    }

    if (directPdfA_)
    {
        *directPdfA_ = m_invArea;
    }

    if (emissionPdfW_)
    {
        *emissionPdfW_ = sampleHemisphereCosinePdfW(m_frame.getNormal(), -_rayDirection);
        *emissionPdfW_ *= m_invArea;
    }

    return m_intensity;
}

vec3 AreaLight::Emit(const SceneSphere &/*aSceneSphere*/, const vec2 &_dirRndTuple, const vec2 &_posRndTuple, vec3 &position_, vec3 &direction_, float &emissionPdfW_, float *directPdfA_, float *cosThetaLight_) const
{
    const vec2 uv = sampleTriangleUniform(_posRndTuple.x, _posRndTuple.y);
    position_ = m_p0 + m_e1 * uv.x + m_e2 * uv.y;

    vec3 localDirOut = sampleHemisphereCosineW(_dirRndTuple.x, _dirRndTuple.y, &emissionPdfW_);

    emissionPdfW_ *= m_invArea;

    // cannot really not emit the particle, so just bias it to the correct angle
    localDirOut.z = f_max(localDirOut.z, NUM_EPS_COSINE);
    direction_ = m_frame.toWorld(localDirOut);

    if (directPdfA_)
    {
        *directPdfA_ = m_invArea;
    }

    if (cosThetaLight_)
    {
        *cosThetaLight_ = localDirOut.z;
    }

    return m_intensity * localDirOut.z;
}

vec3 AreaLight::Illuminate(const SceneSphere &/*aSceneSphere*/, const vec3 &_receivingPosition, const vec2 &_rndTuple, vec3 &directionToLight_, float &distance_, float &directPdfW_, float *emissionPdfW_ /*= nullptr*/, float *cosAtLight_ /*= nullptr*/) const
{
    const vec2 uv = sampleTriangleUniform(_rndTuple.x, _rndTuple.y);
    const vec3 lightPoint = m_p0 + m_e1 * uv.x + m_e2 * uv.y;

    directionToLight_ = lightPoint - _receivingPosition;
    const float distSqr = lengthSquared(directionToLight_);
    distance_ = std::sqrt(distSqr);
    directionToLight_ = directionToLight_ / distance_;

    const float cosNormalDir = dot(m_frame.getNormal(), -directionToLight_);

    // too close to, or under, tangent
    if (cosNormalDir < NUM_EPS_COSINE)
    {
        return vec3(0.0f);
    }

    directPdfW_ = m_invArea * distSqr / cosNormalDir;

    if (cosAtLight_)
    {
        *cosAtLight_ = cosNormalDir;
    }

    if (emissionPdfW_)
    {
        *emissionPdfW_ = m_invArea * cosNormalDir / NUM_PI;
    }

    return m_intensity;
}

AreaLight::AreaLight(const vec3 &_p0, const vec3 &_p1, const vec3 &_p2)
{
    m_p0 = _p0;
    m_e1 = _p1 - _p0;
    m_e2 = _p2 - _p0;

    vec3 normal = cross(m_e1, m_e2);
    float len = length(normal);
    m_invArea = 2.0f / len;
    m_frame.setFromNormal(normalize(normal));
}


bool DirectionalLight::IsDelta() const
{
    return true;
}

bool DirectionalLight::IsFinite() const
{
    return false;
}

vec3 DirectionalLight::GetRadiance(const SceneSphere &/*aSceneSphere*/, const vec3 &/*aRayDirection*/, const vec3 &/*aHitPoint*/, float *directPdfA_ /*= nullptr*/, float *emissionPdfW_ /*= nullptr*/) const
{
    return vec3(0);
}

vec3 DirectionalLight::Emit(const SceneSphere &_sceneSphere, const vec2 &/*aDirRndTuple*/, const vec2 &_posRndTuple, vec3 &position_, vec3 &direction_, float &emissionPdfW_, float *directPdfA_, float *cosThetaLight_) const
{
    const vec2 xy = sampleConcentricDiscA(_posRndTuple.x, _posRndTuple.y);

    position_ = _sceneSphere.m_cenetr + _sceneSphere.m_radius * m_frame.toWorld(vec3(xy.x, xy.y, -1.0f));

    direction_ = m_frame.getNormal();
    emissionPdfW_ = sampleConcentricDiscPdfA() * _sceneSphere.m_inv_radius_sqr;

    if (directPdfA_)
    {
        *directPdfA_ = 1.0f;
    }

    // Not used for infinite or delta lights
    if (cosThetaLight_)
    {
        *cosThetaLight_ = 1.0f;
    }

    return m_intensity;
}

vec3 DirectionalLight::Illuminate(const SceneSphere &_sceneSphere, const vec3 &/*aReceivingPosition*/, const vec2 &/*aRndTuple*/, vec3 &directionToLight_, float &distance_, float &directPdfW_, float *emissionPdfW_ /*= nullptr*/, float *cosAtLight_ /*= nullptr*/) const
{
    directionToLight_ = -m_frame.getNormal();
    distance_ = NUM_INFINITY;
    directPdfW_ = 1.0f;

    if (cosAtLight_)
    {
        *cosAtLight_ = 1.0f;
    }

    if (emissionPdfW_)
    {
        *emissionPdfW_ = sampleConcentricDiscPdfA() * _sceneSphere.m_inv_radius_sqr;
    }

    return m_intensity;
}

DirectionalLight::DirectionalLight(const vec3& aDirection)
{
    m_frame.setFromNormal(normalize(aDirection));
}

bool PointLight::IsDelta() const
{
    return true;
}

bool PointLight::IsFinite() const
{
    return true;
}

vec3 PointLight::GetRadiance(const SceneSphere &/*aSceneSphere*/, const vec3 &/*aRayDirection*/, const vec3 &/*aHitPoint*/, float *directPdfA_ /*= nullptr*/, float *emissionPdfW_ /*= nullptr*/) const
{
    return vec3(0);
}

vec3 PointLight::Emit(const SceneSphere &/*aSceneSphere*/, const vec2 &_dirRndTuple, const vec2 &/*aPosRndTuple*/, vec3 &position_, vec3 &direction_, float &emissionPdfW_, float *directPdfA_, float *oCosThetaLight) const
{
    position_ = m_position;
    direction_ = sampleSphereUniform(_dirRndTuple.x, _dirRndTuple.y, &emissionPdfW_);

    if (directPdfA_)
    {
        *directPdfA_ = 1.0f;
    }

    // Not used for infinite or delta lights
    if (oCosThetaLight)
    {
        *oCosThetaLight = 1.0f;
    }

    return m_intensity;
}

vec3 PointLight::Illuminate(const SceneSphere &/*aSceneSphere*/, const vec3 &_receivingPosition, const vec2 &_rndTuple, vec3 &directionToLight_, float &distance_, float &directPdfW_, float *emissionPdfW_ /*= nullptr*/, float *cosAtLight_ /*= nullptr*/) const
{
    directionToLight_ = m_position - _receivingPosition;
    const float distSqr = lengthSquared(directionToLight_);
    directPdfW_ = distSqr;
    distance_ = std::sqrt(distSqr);
    directionToLight_ = directionToLight_ / distance_;

    if (cosAtLight_)
    {
        *cosAtLight_ = 1.0f;
    }

    if (emissionPdfW_)
    {
        *emissionPdfW_ = sampleSphereUniformPdf();
    }

    return m_intensity;
}

PointLight::PointLight(const vec3& _position)
{
    m_position = _position;
}

bool BackgroundLight::IsDelta() const
{
    return false;
}

bool BackgroundLight::IsFinite() const
{
    return false;
}

vec3 BackgroundLight::GetRadiance(const SceneSphere &_sceneSphere, const vec3 &/*aRayDirection*/, const vec3 &/*aHitPoint*/, float *directPdfA_ /*= nullptr*/, float *emissionPdfW_ /*= nullptr*/) const
{
    // Replace this with image lookup (proper pdf and such)
    // use aRayDirection
    float directPdf = sampleSphereUniformPdf();
    vec3 radiance = m_backgroundColor * m_scale;

    const float positionPdf = sampleConcentricDiscPdfA() *
        _sceneSphere.m_inv_radius_sqr;

    if (directPdfA_)
    {
        *directPdfA_ = directPdf;
    }

    if (emissionPdfW_)
    {
        *emissionPdfW_ = directPdf * positionPdf;
    }

    return radiance;
}

vec3 BackgroundLight::Emit(const SceneSphere &_sceneSphere, const vec2 &_dirRndTuple, const vec2 &_posRndTuple, vec3 &position_, vec3 &direction_, float &emissionPdfW_, float *directPdfA_, float *cosThetaLight_) const
{
    float directPdf;

    // Replace these two lines with image sampling
    direction_ = sampleSphereUniform(_dirRndTuple.x, _dirRndTuple.y, &directPdf);
    //oDirection = -Vec3f(0.16123600f, -0.98195398f, 0.098840252f);
    vec3 radiance = m_backgroundColor * m_scale;

    // Stays even with image sampling
    const vec2 xy = sampleConcentricDiscA(_posRndTuple.x, _posRndTuple.y);

    OrthonormalFrame frame;
    frame.setFromNormal(direction_);

    position_ = _sceneSphere.m_cenetr + _sceneSphere.m_radius * frame.toWorld(vec3(xy.x, xy.y, -1.0f));

    //oPosition = Vec3f(-1.109054f, -2.15064538f, -1.087019148f);

    emissionPdfW_ = directPdf * sampleConcentricDiscPdfA() *
        _sceneSphere.m_inv_radius_sqr;

    // For background we lie about Pdf being in area measure
    if (directPdfA_)
    {
        *directPdfA_ = directPdf;
    }

    // Not used for infinite or delta lights
    if (cosThetaLight_)
    {
        *cosThetaLight_ = 1.0f;
    }

    return radiance;
}

vec3 BackgroundLight::Illuminate(const SceneSphere &_sceneSphere, const vec3 &_receivingPosition, const vec2 &_rndTuple, vec3 &directionToLight_, float &distance_, float &directPdfW_, float *emissionPdfW_ /*= nullptr*/, float *cosAtLight_ /*= nullptr*/) const
{
    // Replace these two lines with image sampling
    directionToLight_ = sampleSphereUniform(_rndTuple.x, _rndTuple.y, &directPdfW_);

    //oDirectionToLight = Vec3f(0.16123600f, -0.98195398f, 0.098840252f);
    vec3 radiance = m_backgroundColor * m_scale;

    // This stays even with image sampling
    distance_ = NUM_INFINITY;
    if (emissionPdfW_)
    {
        *emissionPdfW_ = directPdfW_ * sampleConcentricDiscPdfA() * _sceneSphere.m_inv_radius_sqr;
    }

    if (cosAtLight_)
    {
        *cosAtLight_ = 1.0f;
    }

    return radiance;
}

BackgroundLight::BackgroundLight()
{
    m_backgroundColor = vec3(135, 206, 250) / vec3(255.0f);
    m_scale = 1.0f;
}

bool SphericalLight::IsDelta() const
{
    return false;
}

bool SphericalLight::IsFinite() const
{
    return true;
}

vec3 SphericalLight::GetRadiance(const SceneSphere &/*aSceneSphere*/, const vec3 &_rayDirection, const vec3 &_hitPoint, float *directPdfA_ /*= nullptr*/, float *emissionPdfW_ /*= nullptr*/) const
{
    const vec3 normal = normalize(_hitPoint - m_pos);
    const float cosOutL = f_max(0.0f, dot(normal, -_rayDirection));

    if (cosOutL <= 0)
    {
        return vec3(0);
    }

    if (directPdfA_)
    {
        *directPdfA_ = m_invArea;
    }

    if (emissionPdfW_)
    {
        *emissionPdfW_ = sampleHemisphereCosinePdfW(normal, -_rayDirection);
        *emissionPdfW_ *= m_invArea;
    }

    return m_intensity;
}

vec3 SphericalLight::Emit(const SceneSphere &/*aSceneSphere*/, const vec2 &_dirRndTuple, const vec2 &_posRndTuple, vec3 &position_, vec3 &direction_, float &emissionPdfW_, float *directPdfA_, float *cosThetaLight_) const
{
    position_ = m_pos + m_radius * sampleSphereUniform(_posRndTuple.x, _posRndTuple.y);
    const vec3 normal = normalize(position_ - m_pos);
    OrthonormalFrame frame;
    frame.setFromNormal(normal);

    vec3 localDirOut = sampleHemisphereCosineW(_dirRndTuple.x, _dirRndTuple.y, &emissionPdfW_);

    emissionPdfW_ *= m_invArea;

    // cannot really not emit the particle, so just bias it to the correct angle
    localDirOut.z = f_max(localDirOut.z, NUM_EPS_COSINE);
    direction_ = frame.toWorld(localDirOut);

    if (directPdfA_)
    {
        *directPdfA_ = m_invArea;
    }

    if (cosThetaLight_)
    {
        *cosThetaLight_ = localDirOut.z;
    }

    return m_intensity * localDirOut.z;
}

vec3 SphericalLight::Illuminate(const SceneSphere &/*aSceneSphere*/, const vec3 &_receivingPosition, const vec2 &_rndTuple, vec3 &directionToLight_, float &distance_, float &directPdfW_, float *emissionPdfW_ /*= nullptr*/, float *cosAtLight_ /*= nullptr*/) const
{
    const vec3 lightPoint = m_pos + m_radius * sampleSphereUniform(_rndTuple.x, _rndTuple.y);
    const vec3 normal = normalize(lightPoint - m_pos);

    directionToLight_ = lightPoint - _receivingPosition;
    const float distSqr = lengthSquared(directionToLight_);
    distance_ = std::sqrt(distSqr);
    directionToLight_ = directionToLight_ / distance_;

    const float cosNormalDir = dot(normal, -directionToLight_);

    // too close to, or under, tangent
    if (cosNormalDir < NUM_EPS_COSINE)
    {
        return vec3(0.0f);
    }

    directPdfW_ = m_invArea * distSqr / cosNormalDir; // inlined PdfAtoW

    if (cosAtLight_)
        *cosAtLight_ = cosNormalDir;

    if (emissionPdfW_)
        *emissionPdfW_ = m_invArea * cosNormalDir / NUM_PI;

    return m_intensity;
}

SphericalLight::SphericalLight(const vec3 &_pos, const float _rad)
{
    m_pos = _pos;
    m_radius = _rad;

    m_invArea = 1.0f / (4.0f * NUM_PI * sq(m_radius));
}
