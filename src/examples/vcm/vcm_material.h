// Disclaimer: This demo is adapted from smallvcm http://www.smallvcm.com/

#pragma once

#include "vectors.h"
#include "geometry.h"
#include "ray.h"
#include "sample.h"

class Material_debug
{
public:
    Material_debug()
    {
        Reset();
    }

    void Reset()
    {
        m_diffuseReflectance = vec3(0);
        m_phongReflectance = vec3(0);
        m_phongExponent = 1.0f;
        m_mirrorReflectance = vec3(0);
        m_IOR = -1.0f;
    }

    // diffuse is simply added to the others
    vec3 m_diffuseReflectance;

    // Phong is simply added to the others
    vec3 m_phongReflectance;
    float m_phongExponent;

    // mirror can be either simply added, or mixed using Fresnel term
    // this is governed by mIOR, if it is >= 0, fresnel is used, otherwise
    // it is not
    vec3 m_mirrorReflectance;

    // When mIOR >= 0, we also transmit (just clear glass)
    float m_IOR;
};

//////////////////////////////////////////////////////////////////////////
// BSDF, most magic happens here
//
// One of important conventions is prefixing direction with World when they
// are in world coordinates and with Local when they are in local frame,
// i.e., mFrame.
//
// Another important convention is suffix Fix and Gen.
// For PDF computation, we need to know which direction is given (Fix),
// and which is the generated (Gen) direction. This is important even
// when simply evaluating BSDF.
// In BPT, we call Evaluate() when directly connecting to light/camera.
// This gives us both directions required for evaluating BSDF.
// However, for MIS we also need to know probabilities of having sampled
// this path via BSDF sampling, and we need that for both possible directions.
// The Fix/Gen convention (along with Direct and Reverse for PDF) clearly
// establishes which PDF is which.
//
// The BSDF is also templated by direction of tracing, whether from camera
// (BSDF<false>) or from light (BSDF<true>). This is identical to Veach's
// Adjoint BSDF (except the name is more straightforward).
// For us this is only used when refracting.

#define EPS_PHONG 1e-3f

template<bool FixIsLight>
class BSDF
{
private:
    struct ComponentProbabilities
    {
        float m_diffProb;
        float m_phongProb;
        float m_reflProb;
        float m_refrProb;
    };

private:
    const Material_debug *m_material = nullptr;  //!< * of scene material, nullptr ~ invalid
    OrthonormalFrame m_frame;            //!< Local frame of reference
    vec3 m_localDirFix;      //!< Incoming (fixed) direction, in local
    bool  m_isDelta;          //!< True when material is purely specular
    ComponentProbabilities m_probabilities; //!< Sampling probabilities
    float m_continuationProb; //!< Russian roulette probability
    float m_reflectCoeff;     //!< Fresnel reflection coefficient (for glass)

public:
    enum Events
    {
        kNONE = 0,
        kDiffuse = 1,
        kPhong = 2,
        kReflect = 4,
        kRefract = 8,
        kSpecular = (kReflect | kRefract),
        kNonSpecular = (kDiffuse | kPhong),
        kAll = (kSpecular | kNonSpecular)
    };

public:
    BSDF() {};

    BSDF(
        const Ray&      _ray,
        const vec3&     _normal,
        const Material_debug  *_mat)
    {
        m_frame.setFromNormal(_normal);
        m_localDirFix = m_frame.toLocal(-_ray.dir);

        // reject rays that are too parallel with tangent plane
        if (std::abs(m_localDirFix.z) < NUM_EPS_COSINE)
        {
            return;
        }

        m_material = _mat;
        const Material_debug &mat = *m_material;
        GetComponentProbabilities(mat, m_probabilities);

        m_isDelta = (m_probabilities.m_diffProb == 0) && (m_probabilities.m_phongProb == 0);
    }

    /* \brief Given a direction, evaluates BSDF
    *
    * Returns value of BSDF, as well as cosine for the
    * aWorldDirGen direction.
    * Can return probability (w.r.t. solid angle W),
    * of having sampled aWorldDirGen given mLocalDirFix (oDirectPdfW),
    * and of having sampled mLocalDirFix given aWorldDirGen (oReversePdfW).
    *
    */
    vec3 Evaluate(
        const vec3&     _worldDirGen,
        float&          cosThetaGen_,
        float           *directPdfW_ = nullptr,
        float           *reversePdfW_ = nullptr) const
    {
        vec3 result(0);

        if (directPdfW_)
        {
            *directPdfW_ = 0;
        }
        if (reversePdfW_)
        {
            *reversePdfW_ = 0;
        }

        const vec3 localDirGen = m_frame.toLocal(_worldDirGen);

        if (localDirGen.z * m_localDirFix.z < 0)
        {
            return result;
        }

        cosThetaGen_ = std::abs(localDirGen.z);

        const Material_debug &mat = *m_material;

        result += EvaluateDiffuse(mat, localDirGen, directPdfW_, reversePdfW_);
        result += EvaluatePhong(mat, localDirGen, directPdfW_, reversePdfW_);

        return result;
    }

    /* \brief Given a direction, evaluates Pdf
    *
    * By default returns PDF with which would be aWorldDirGen
    * generated from mLocalDirFix. When aEvalRevPdf == true,
    * it provides PDF for the reverse direction.
    */
    float Pdf(
        const vec3&     _worldDirGen,
        const bool      _evalRevPdf = false) const
    {
        const vec3 localDirGen = m_frame.toLocal(_worldDirGen);

        if (localDirGen.z * m_localDirFix.z < 0)
        {
            return 0;
        }

        const Material_debug &mat = *m_material;

        float directPdfW = 0;
        float reversePdfW = 0;

        PdfDiffuse(mat, localDirGen, &directPdfW, &reversePdfW);
        PdfPhong(mat, localDirGen, &directPdfW, &reversePdfW);

        return _evalRevPdf ? reversePdfW : directPdfW;
    }

    /* \brief Given 3 random numbers, samples new direction from BSDF.
    *
    * Uses z component of random triplet to pick BSDF component from
    * which it will sample direction. If non-specular component is chosen,
    * it will also evaluate the other (non-specular) BSDF components.
    * Return BSDF factor for given direction, as well as PDF choosing that direction.
    * Can return event which has been sampled.
    * If result is Vec3f(0,0,0), then the sample should be discarded.
    */
    vec3 Sample(
        const vec3&     _rndTriplet,
        vec3&           worldDirGen_,
        float&          pdfW_,
        float&          cosThetaGen_,
        uint32_t        *sampledEvent_ = nullptr) const
    {
        uint32_t sampledEvent;

        if (_rndTriplet.z < m_probabilities.m_diffProb)
        {
            sampledEvent = kDiffuse;
        }
        else if (_rndTriplet.z < m_probabilities.m_diffProb + m_probabilities.m_phongProb)
        {
            sampledEvent = kPhong;
        }
        else if (_rndTriplet.z < m_probabilities.m_diffProb + m_probabilities.m_phongProb + m_probabilities.m_reflProb)
        {
            sampledEvent = kReflect;
        }
        else
        {
            sampledEvent = kRefract;
        }

        if (sampledEvent_)
        {
            *sampledEvent_ = sampledEvent;
        }

        const Material_debug &mat = *m_material;

        pdfW_ = 0;
        vec3 result(0);
        vec3 localDirGen;

        if (sampledEvent == kDiffuse)
        {
            result += SampleDiffuse(mat, vec2(_rndTriplet.x, _rndTriplet.y, 0), localDirGen, pdfW_);

            if (isZero(result))
            {
                return vec3(0);
            }

            result += EvaluatePhong(mat, localDirGen, &pdfW_);
        }
        else if (sampledEvent == kPhong)
        {
            result += SamplePhong(mat, vec2(_rndTriplet.x, _rndTriplet.y, 0), localDirGen, pdfW_);

            if (isZero(result))
            {
                return vec3(0);
            }

            result += EvaluateDiffuse(mat, localDirGen, &pdfW_);
        }
        else if (sampledEvent == kReflect)
        {
            result += SampleReflect(mat, vec2(_rndTriplet.x, _rndTriplet.y, 0), localDirGen, pdfW_);

            if (isZero(result))
            {
                return vec3(0);
            }
        }
        else
        {
            result += SampleRefract(mat, vec2(_rndTriplet.x, _rndTriplet.y, 0), localDirGen, pdfW_);
            if (isZero(result))
            {
                return vec3(0);
            }
        }

        cosThetaGen_ = std::abs(localDirGen.z);
        if (cosThetaGen_ < NUM_EPS_COSINE)
        {
            return vec3(0.0f);
        }

        worldDirGen_ = m_frame.toWorld(localDirGen);
        return result;
    }


    bool  IsValid() const           { return nullptr != m_material; }
    bool  IsDelta() const           { return m_isDelta; }
    float ContinuationProb() const  { return m_continuationProb; }
    float CosThetaFix() const       { return m_localDirFix.z; }
    vec3 WorldDirFix() const       { return m_frame.toWorld(m_localDirFix); }

private:

    ////////////////////////////////////////////////////////////////////////////
    // Sampling methods
    // All sampling methods take material, 2 random numbers [0-1[,
    // and return BSDF factor, generated direction in local coordinates, and PDF
    ////////////////////////////////////////////////////////////////////////////

    vec3 SampleDiffuse(
        const Material_debug& _material,
        const vec2&     _rndTuple,
        vec3&           localDirGen_,
        float&          pdfW_) const
    {
        if (m_localDirFix.z < NUM_EPS_COSINE)
            return vec3(0);

        float unweightedPdfW;
        localDirGen_ = sampleHemisphereCosineW(_rndTuple.x, _rndTuple.y, &unweightedPdfW);
        pdfW_ += unweightedPdfW * m_probabilities.m_diffProb;

        return _material.m_diffuseReflectance / NUM_PI;
    }

    vec3 SamplePhong(
        const Material_debug& _material,
        const vec2&     _rndTuple,
        vec3&           localDirGen_,
        float&          pdfW_) const
    {
        localDirGen_ = sampleHemisphereCosinePowerW(_rndTuple.x, _rndTuple.y, _material.m_phongExponent, nullptr);

        // Due to numeric issues in MIS, we actually need to compute all pdfs
        // exactly the same way all the time!!!
        const vec3 reflLocalDirFixed = reflectLocal(m_localDirFix);
        {
            OrthonormalFrame frame;
            frame.setFromNormal(reflLocalDirFixed);
            localDirGen_ = frame.toWorld(localDirGen_);
        }

        const float dot_R_Wi = dot(reflLocalDirFixed, localDirGen_);

        if (dot_R_Wi <= EPS_PHONG)
        {
            return vec3(0.0f);
        }

        PdfPhong(_material, localDirGen_, &pdfW_);

        const vec3 rho = _material.m_phongReflectance *
            (_material.m_phongExponent + 2.0f) * 0.5f / NUM_PI;

        return rho * std::pow(dot_R_Wi, _material.m_phongExponent);
    }

    vec3 SampleReflect(
        const Material_debug& _material,
        const vec2&     _rndTuple,
        vec3&           localDirGen_,
        float&          pdfW_) const
    {
        localDirGen_ = reflectLocal(m_localDirFix);

        pdfW_ += m_probabilities.m_reflProb;
        // BSDF is multiplied (outside) by cosine (oLocalDirGen.z),
        // for mirror this shouldn't be done, so we pre-divide here instead
        return m_reflectCoeff * _material.m_mirrorReflectance /
            std::abs(localDirGen_.z);
    }

    vec3 SampleRefract(
        const Material_debug& _material,
        const vec2&     _rndTuple,
        vec3&           localDirGen_,
        float&          pdfW_) const
    {
        if (_material.m_IOR < 0)
        {
            return vec3(0);
        }

        float cosI = m_localDirFix.z;

        float cosT;
        float etaIncOverEtaTrans;

        if (cosI < 0.0f) // hit from inside
        {
            etaIncOverEtaTrans = _material.m_IOR;
            cosI = -cosI;
            cosT = 1.0f;
        }
        else
        {
            etaIncOverEtaTrans = 1.0f / _material.m_IOR;
            cosT = -1.0f;
        }

        const float sinI2 = 1.0f - cosI * cosI;
        const float sinT2 = sq(etaIncOverEtaTrans) * sinI2;

        if (sinT2 < 1.0f) // no total internal reflection
        {
            cosT *= std::sqrt(f_max(0.0f, 1.0f - sinT2));

            localDirGen_ = vec3(
                -etaIncOverEtaTrans * m_localDirFix.x,
                -etaIncOverEtaTrans * m_localDirFix.y,
                cosT);

            pdfW_ += m_probabilities.m_refrProb;

            const float refractCoeff = 1.0f - m_reflectCoeff;
            // only camera paths are multiplied by this factor, and etas
            // are swapped because radiance flows in the opposite direction
            if (!FixIsLight)
            {
                return vec3(refractCoeff * sq(etaIncOverEtaTrans) / std::abs(cosT));
            }
            else
            {
                return vec3(refractCoeff / std::abs(cosT));
            }
        }
        //else total internal reflection, do nothing
        // because it is sampled as reflection ??

        pdfW_ += 0.0f;
        return vec3(0.0f);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Evaluation methods
    ////////////////////////////////////////////////////////////////////////////

    vec3 EvaluateDiffuse(
        const Material_debug& _material,
        const vec3&     _localDirGen,
        float           *directPdfW_ = nullptr,
        float           *reversePdfW_ = nullptr) const
    {
        if (m_probabilities.m_diffProb == 0)
        {
            return vec3(0);
        }

        if (m_localDirFix.z < NUM_EPS_COSINE || _localDirGen.z < NUM_EPS_COSINE)
        {
            return vec3(0);
        }

        if (directPdfW_)
        {
            *directPdfW_ += m_probabilities.m_diffProb * f_max(0.0f, _localDirGen.z / NUM_PI);
        }

        if (reversePdfW_)
        {
            *reversePdfW_ += m_probabilities.m_diffProb * f_max(0.0f, m_localDirFix.z / NUM_PI);
        }

        return _material.m_diffuseReflectance / NUM_PI;
    }

    vec3 EvaluatePhong(
        const Material_debug& _material,
        const vec3&     _localDirGen,
        float           *directPdfW_ = nullptr,
        float           *reversePdfW_ = nullptr) const
    {
        if (m_probabilities.m_phongProb == 0)
        {
            return vec3(0);
        }

        if (m_localDirFix.z < NUM_EPS_COSINE || _localDirGen.z < NUM_EPS_COSINE)
        {
            return vec3(0);
        }

        // assumes this is never called when rejectShadingCos(oLocalDirGen.z) is true
        const vec3 reflLocalDirIn = reflectLocal(m_localDirFix);
        const float dot_R_Wi = dot(reflLocalDirIn, _localDirGen);

        if (dot_R_Wi <= EPS_PHONG)
        {
            return vec3(0.0f);
        }

        if (directPdfW_ || reversePdfW_)
        {
            // the sampling is symmetric
            const float pdfW = m_probabilities.m_phongProb *
                sampleHemisphereCosinePowerPdfW(reflLocalDirIn, _localDirGen, _material.m_phongExponent);

            if (directPdfW_)
            {
                *directPdfW_ += pdfW;
            }

            if (reversePdfW_)
            {
                *reversePdfW_ += pdfW;
            }
        }

        const vec3 rho = _material.m_phongReflectance *
            (_material.m_phongExponent + 2.0f) * 0.5f / NUM_PI;

        return rho * std::pow(dot_R_Wi, _material.m_phongExponent);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Pdf methods
    ////////////////////////////////////////////////////////////////////////////

    void PdfDiffuse(
        const Material_debug& _material,
        const vec3&     _localDirGen,
        float           *directPdfW_ = nullptr,
        float           *reversePdfW_ = nullptr) const
    {
        if (m_probabilities.m_diffProb == 0)
        {
            return;
        }

        if (directPdfW_)
        {
            *directPdfW_ += m_probabilities.m_diffProb *
                f_max(0.0f, _localDirGen.z / NUM_PI);
        }

        if (reversePdfW_)
        {
            *reversePdfW_ += m_probabilities.m_diffProb *
                f_max(0.0f, m_localDirFix.z / NUM_PI);
        }
    }

    void PdfPhong(
        const Material_debug& _material,
        const vec3&     _localDirGen,
        float           *directPdfW_ = nullptr,
        float           *reversePdfW_ = nullptr) const
    {
        if (m_probabilities.m_phongProb == 0)
        {
            return;
        }

        // assumes this is never called when rejectShadingCos(oLocalDirGen.z) is true
        const vec3 reflLocalDirIn = reflectLocal(m_localDirFix);
        const float dot_R_Wi = dot(reflLocalDirIn, _localDirGen);

        if (dot_R_Wi <= EPS_PHONG)
        {
            return;
        }

        if (directPdfW_ || reversePdfW_)
        {
            // the sampling is symmetric
            const float pdfW = sampleHemisphereCosinePowerPdfW(reflLocalDirIn, _localDirGen,
                _material.m_phongExponent) * m_probabilities.m_phongProb;

            if (directPdfW_)
            {
                *directPdfW_ += pdfW;
            }

            if (reversePdfW_)
            {
                *reversePdfW_ += pdfW;
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Albedo methods
    ////////////////////////////////////////////////////////////////////////////

    float AlbedoDiffuse(const Material_debug& _material) const
    {
        return luminance(_material.m_diffuseReflectance);
    }

    float AlbedoPhong(const Material_debug& _material) const
    {
        return luminance(_material.m_phongReflectance);
    }

    float AlbedoReflect(const Material_debug& _material) const
    {
        return luminance(_material.m_mirrorReflectance);
    }

    float AlbedoRefract(const Material_debug& _material) const
    {
        return _material.m_IOR > 0.0f ? 1.0f : 0.0f;
    }

    void GetComponentProbabilities(
        const Material_debug         &_material,
        ComponentProbabilities &probabilities_)
    {
        m_reflectCoeff = fresnelReflectance(m_localDirFix.z, _material.m_IOR);

        const float albedoDiffuse = AlbedoDiffuse(_material);
        const float albedoPhong = AlbedoPhong(_material);
        const float albedoReflect = m_reflectCoeff         * AlbedoReflect(_material);
        const float albedoRefract = (1.0f - m_reflectCoeff) * AlbedoRefract(_material);

        const float totalAlbedo = albedoDiffuse + albedoPhong + albedoReflect + albedoRefract;

        if (totalAlbedo < 1e-9f)
        {
            probabilities_.m_diffProb = 0.0f;
            probabilities_.m_phongProb = 0.0f;
            probabilities_.m_reflProb = 0.0f;
            probabilities_.m_refrProb = 0.0f;
            m_continuationProb = 0.0f;
        }
        else
        {
            probabilities_.m_diffProb = albedoDiffuse / totalAlbedo;
            probabilities_.m_phongProb = albedoPhong / totalAlbedo;
            probabilities_.m_reflProb = albedoReflect / totalAlbedo;
            probabilities_.m_refrProb = albedoRefract / totalAlbedo;
            // The continuation probability is max component from reflectance.
            // That way the weight of sample will never rise.
            // Luminance is another very valid option.
            m_continuationProb =
                maxComponent(_material.m_diffuseReflectance +
                _material.m_phongReflectance +
                m_reflectCoeff * _material.m_mirrorReflectance) +
                (1.0f - m_reflectCoeff);

            m_continuationProb = clampf(m_continuationProb, 0.0f, 1.0f);
        }
    }
};
