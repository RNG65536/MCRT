#ifndef light_h__
#define light_h__

#include "geometry.h"
#include "vectors.h"

class SceneSphere;

// grabbed from smallvcm
// http://www.smallvcm.com/
class AbstractLight
{
public:
    /* \brief Illuminates a given point in the scene.
     *
     * Given a point and two random samples (e.g., for position on area lights),
     * this method returns direction from point to light, distance,
     * pdf of having chosen this direction (e.g., 1 / area).
     * Optionally also returns pdf of emitting particle in this direction,
     * and cosine from lights normal (helps with PDF of hitting the light,
     * but set to 1 for point lights).
     *
     * Returns radiance.
     */
    virtual vec3 Illuminate(const SceneSphere& _sceneSphere,
                            const vec3&        _receivingPosition,
                            const vec2&        _RndTuple,
                            vec3&              directionToLight_,
                            float&             distance_,
                            float&             directPdfW_,
                            float*             emissionPdfW_ = nullptr,
                            float*             cosAtLight_ = nullptr) const = 0;

    /* \brief Emits particle from the light.
     *
     * Given two sets of random numbers (e.g., position and direction on area
     * light), this method generates a position and direction for light
     * particle, along with the pdf.
     *
     * Can also supply pdf (w.r.t. area) of choosing this position when calling
     * Illuminate. Also provides cosine on the light (this is 1 for point lights
     * etc.).
     *
     * Returns "energy" that particle carries
     */
    virtual vec3 Emit(const SceneSphere& _sceneSphere,
                      const vec2&        _dirRndTuple,
                      const vec2&        _posRndTuple,
                      vec3&              position_,
                      vec3&              direction_,
                      float&             emissionPdfW_,
                      float*             directPdfA_,
                      float*             cosThetaLight_) const = 0;

    /* \brief Returns radiance for ray randomly hitting the light
     *
     * Given ray direction and hitpoint, it returns radiance.
     * Can also provide area pdf of sampling hitpoint in Illuminate,
     * and of emitting particle along the ray (in opposite direction).
     */
    virtual vec3 GetRadiance(const SceneSphere& _sceneSphere,
                             const vec3&        _rayDirection,
                             const vec3&        _hitPoint,
                             float*             directPdfA_ = nullptr,
                             float* emissionPdfW_ = nullptr) const = 0;

    // Whether the light has a finite extent (area, point) or not (directional,
    // env. map)
    virtual bool IsFinite() const = 0;

    // Whether the light has delta function (point, directional) or not (area)
    virtual bool IsDelta() const = 0;
};

//////////////////////////////////////////////////////////////////////////
class AreaLight : public AbstractLight
{
public:
    vec3             m_p0, m_e1, m_e2;
    OrthonormalFrame m_frame;
    vec3             m_intensity;
    float            m_invArea;

public:
    AreaLight(const vec3& _p0, const vec3& _p1, const vec3& _p2);

    virtual vec3 Illuminate(const SceneSphere& /*aSceneSphere*/,
                            const vec3& _receivingPosition,
                            const vec2& _rndTuple,
                            vec3&       directionToLight_,
                            float&      distance_,
                            float&      directPdfW_,
                            float*      emissionPdfW_ = nullptr,
                            float*      cosAtLight_ = nullptr) const;

    virtual vec3 Emit(const SceneSphere& /*aSceneSphere*/,
                      const vec2& _dirRndTuple,
                      const vec2& _posRndTuple,
                      vec3&       position_,
                      vec3&       direction_,
                      float&      emissionPdfW_,
                      float*      directPdfA_,
                      float*      cosThetaLight_) const;

    virtual vec3 GetRadiance(const SceneSphere& /*aSceneSphere*/,
                             const vec3& _rayDirection,
                             const vec3& _hitPoint,
                             float*      directPdfA_ = nullptr,
                             float*      emissionPdfW_ = nullptr) const;
    // Whether the light has a finite extent (area, point) or not (directional,
    // env. map)
    virtual bool IsFinite() const;

    // Whether the light has delta function (point, directional) or not (area)
    virtual bool IsDelta() const;
};

//////////////////////////////////////////////////////////////////////////
class SphericalLight : public AbstractLight
{
public:
    SphericalLight(const vec3& _pos, const float _rad);

    virtual vec3 Illuminate(const SceneSphere& /*aSceneSphere*/,
                            const vec3& _receivingPosition,
                            const vec2& _rndTuple,
                            vec3&       directionToLight_,
                            float&      distance_,
                            float&      directPdfW_,
                            float*      emissionPdfW_ = nullptr,
                            float*      cosAtLight_ = nullptr) const;

    virtual vec3 Emit(const SceneSphere& /*aSceneSphere*/,
                      const vec2& _dirRndTuple,
                      const vec2& _posRndTuple,
                      vec3&       position_,
                      vec3&       direction_,
                      float&      emissionPdfW_,
                      float*      directPdfA_,
                      float*      cosThetaLight_) const;

    virtual vec3 GetRadiance(const SceneSphere& /*aSceneSphere*/,
                             const vec3& _rayDirection,
                             const vec3& _hitPoint,
                             float*      directPdfA_ = nullptr,
                             float*      emissionPdfW_ = nullptr) const;
    // Whether the light has a finite extent (area, point) or not (directional,
    // env. map)
    virtual bool IsFinite() const;

    // Whether the light has delta function (point, directional) or not (area)
    virtual bool IsDelta() const;

public:
    vec3  m_pos;
    float m_radius;
    float m_invArea;
    vec3  m_intensity;
};

//////////////////////////////////////////////////////////////////////////
class DirectionalLight : public AbstractLight
{
public:
    OrthonormalFrame m_frame;
    vec3             m_intensity;

public:
    DirectionalLight(const vec3& aDirection);

    virtual vec3 Illuminate(const SceneSphere& _sceneSphere,
                            const vec3& /*aReceivingPosition*/,
                            const vec2& /*aRndTuple*/,
                            vec3&  directionToLight_,
                            float& distance_,
                            float& directPdfW_,
                            float* emissionPdfW_ = nullptr,
                            float* cosAtLight_ = nullptr) const;

    virtual vec3 Emit(const SceneSphere& _sceneSphere,
                      const vec2& /*aDirRndTuple*/,
                      const vec2& _posRndTuple,
                      vec3&       position_,
                      vec3&       direction_,
                      float&      emissionPdfW_,
                      float*      directPdfA_,
                      float*      cosThetaLight_) const;

    virtual vec3 GetRadiance(const SceneSphere& /*aSceneSphere*/,
                             const vec3& /*aRayDirection*/,
                             const vec3& /*aHitPoint*/,
                             float* directPdfA_ = nullptr,
                             float* emissionPdfW_ = nullptr) const;

    // Whether the light has a finite extent (area, point) or not (directional,
    // env. map)
    virtual bool IsFinite() const;

    // Whether the light has delta function (point, directional) or not (area)
    virtual bool IsDelta() const;
};

//////////////////////////////////////////////////////////////////////////
class PointLight : public AbstractLight
{
public:
    vec3 m_position;
    vec3 m_intensity;

public:
    PointLight(const vec3& _position);

    virtual vec3 Illuminate(const SceneSphere& /*aSceneSphere*/,
                            const vec3& _receivingPosition,
                            const vec2& _rndTuple,
                            vec3&       directionToLight_,
                            float&      distance_,
                            float&      directPdfW_,
                            float*      emissionPdfW_ = nullptr,
                            float*      cosAtLight_ = nullptr) const;

    virtual vec3 Emit(const SceneSphere& /*aSceneSphere*/,
                      const vec2& _dirRndTuple,
                      const vec2& /*aPosRndTuple*/,
                      vec3&  position_,
                      vec3&  direction_,
                      float& emissionPdfW_,
                      float* directPdfA_,
                      float* oCosThetaLight) const;

    virtual vec3 GetRadiance(const SceneSphere& /*aSceneSphere*/,
                             const vec3& /*aRayDirection*/,
                             const vec3& /*aHitPoint*/,
                             float* directPdfA_ = nullptr,
                             float* emissionPdfW_ = nullptr) const;

    // Whether the light has a finite extent (area, point) or not (directional,
    // env. map)
    virtual bool IsFinite() const;

    // Whether the light has delta function (point, directional) or not (area)
    virtual bool IsDelta() const;
};

//////////////////////////////////////////////////////////////////////////
class BackgroundLight : public AbstractLight
{
public:
    vec3  m_backgroundColor;
    float m_scale;

public:
    BackgroundLight();

    virtual vec3 Illuminate(const SceneSphere& _sceneSphere,
                            const vec3&        _receivingPosition,
                            const vec2&        _rndTuple,
                            vec3&              directionToLight_,
                            float&             distance_,
                            float&             directPdfW_,
                            float*             emissionPdfW_ = nullptr,
                            float*             cosAtLight_ = nullptr) const;

    virtual vec3 Emit(const SceneSphere& _sceneSphere,
                      const vec2&        _dirRndTuple,
                      const vec2&        _posRndTuple,
                      vec3&              position_,
                      vec3&              direction_,
                      float&             emissionPdfW_,
                      float*             directPdfA_,
                      float*             cosThetaLight_) const;

    virtual vec3 GetRadiance(const SceneSphere& _sceneSphere,
                             const vec3& /*aRayDirection*/,
                             const vec3& /*aHitPoint*/,
                             float* directPdfA_ = nullptr,
                             float* emissionPdfW_ = nullptr) const;

    // Whether the light has a finite extent (area, point) or not (directional,
    // env. map)
    virtual bool IsFinite() const;

    // Whether the light has delta function (point, directional) or not (area)
    virtual bool IsDelta() const;
};

#endif  // light_h__
