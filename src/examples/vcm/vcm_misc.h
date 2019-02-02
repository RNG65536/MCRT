// Disclaimer: This demo is adapted from smallvcm http://www.smallvcm.com/

#pragma once

#include "film.h"

class EyeLight
{
public:
    uint32_t            m_maxPathLength;
    uint32_t            m_minPathLength;

    uint32_t            m_iterations;
    FrameBuffer         m_framebuffer;
    const Scene_debug&  m_scene;

public:

    EyeLight(const Scene_debug& _scene) : m_scene(_scene)
    {
        m_minPathLength = 0;
        m_maxPathLength = 10;
        m_iterations = 0;
        m_framebuffer.resize(m_scene.m_camera->filmWidth(), m_scene.m_camera->filmHeight());
    }

    // final output of this rendering work thread
    void GetFramebuffer(FrameBuffer& framebuffer_)
    {
        framebuffer_ = m_framebuffer;

        if (m_iterations > 0)
        {
            framebuffer_.scale(1.0f / m_iterations);
        }
    }

    virtual void RunIteration(int _iteration)
    {
        const int resX = int(m_scene.m_camera->filmWidth());
        const int resY = int(m_scene.m_camera->filmHeight());

        for (int pixID = 0; pixID < resX * resY; pixID++)
        {
            //////////////////////////////////////////////////////////////////////////
            // Generate ray
            const int x = pixID % resX;
            const int y = pixID / resX;

            const vec2 sample = (vec2(float(x), float(y), 0) + vec2(randf(), randf(), 0)) / vec3(float(resX), float(resY), 1);

            Ray   ray = m_scene.m_camera->makeRay(sample.x, sample.y);

            HitInfo_debug hit = m_scene.intersect(ray);
            if (hit)
            {
                float dotLN = dot(hit.shadingNormal(), -ray.dir);

                if (dotLN > 0)
                {
                    m_framebuffer.accumulatePixel(x, y, vec3(dotLN));
                }
                else
                {
                    m_framebuffer.accumulatePixel(x, y, vec3(-dotLN, 0, 0));
                }
            }
        }

        m_iterations++;
    }
};

class PathTracer
{
public:
    uint32_t            m_maxPathLength;
    uint32_t            m_minPathLength;

    uint32_t            m_iterations;
    FrameBuffer         m_framebuffer;
    const Scene_debug&  m_scene;

private:
    // Mis power (1 for balance heuristic)
    float Mis(float _pdf) const
    {
        return _pdf;
    }

    // Mis weight for 2 pdfs
    float Mis2(
        float _samplePdf,
        float _otherPdf) const
    {
        return Mis(_samplePdf) / (Mis(_samplePdf) + Mis(_otherPdf));
    }

public:
    PathTracer(const Scene_debug& _scene) : m_scene(_scene)
    {
        m_minPathLength = 0;
        m_maxPathLength = 10; // 2;
        m_iterations = 0;
        m_framebuffer.resize(m_scene.m_camera->filmWidth(), m_scene.m_camera->filmHeight());
    }

    // final output of this rendering work thread
    void GetFramebuffer(FrameBuffer& framebuffer_)
    {
        framebuffer_ = m_framebuffer;

        if (m_iterations > 0)
        {
            framebuffer_.scale(1.0f / m_iterations);
        }
    }

    void RunIteration(int _iteration)
    {
        // We sample lights uniformly
        const int   lightCount = m_scene.GetLightCount();
        const float lightPickProb = 1.0f / lightCount;

        const int resX = int(m_scene.m_camera->filmWidth());
        const int resY = int(m_scene.m_camera->filmHeight());

        for (int pixID = 0; pixID < resX * resY; pixID++)
        {
            const int x = pixID % resX;
            const int y = pixID / resX;

            const vec2 sample = (vec2(float(x), float(y), 0) + vec2(randf(), randf(), 0)) / vec3(float(resX), float(resY), 1);

            Ray   ray = m_scene.m_camera->makeRay(sample.x, sample.y);

            vec3 pathWeight(1.0f); // throughput
            vec3 color(0.0f);
            uint32_t  pathLength = 1;
            bool  lastSpecular = true;
            float lastPdfW = 1;

            for (;; ++pathLength)
            {
#pragma region
                HitInfo_debug hit = m_scene.intersect(ray);
                if (!hit)
                {
                    if (pathLength < m_minPathLength)
                    {
                        break;
                    }

                    const BackgroundLight* background = m_scene.GetBackground();
                    if (!background)
                    {
                        break;
                    }
                    // For background we cheat with the A/W suffixes,
                    // and GetRadiance actually returns W instead of A
                    float directPdfW;
                    vec3 contrib = background->GetRadiance(m_scene.m_sceneSphere,
                        ray.dir, vec3(0), &directPdfW);
                    if (isZero(contrib))
                    {
                        break;
                    }

                    float misWeight = 1.0f;
                    if (pathLength > 1 && !lastSpecular)
                    {
                        misWeight = Mis2(lastPdfW, directPdfW * lightPickProb);
                    }

                    color += pathWeight * misWeight * contrib;
                    break;
                }
#pragma endregion background

                vec3 hitPoint = ray.at(hit.distance());
                const Material_debug& mat = m_scene.m_materials[hit.materialID()];

                vec3 nhit = hit.shadingNormal();
//                 if (dot(ray.d, nhit) > 0)
//                 {
//                     nhit = -nhit;
//                 }
                BSDF<false> bsdf(ray, nhit, &mat);
                bool isDelta = bsdf.IsDelta();


                if (!bsdf.IsValid())
                {
                    break;
                }

#pragma region
                int light_id = m_scene.toLightID(hit.materialID());

                // directly hit some light, lights do not reflect
                if (light_id >= 0)
                {
                    if (pathLength < m_minPathLength)
                    {
                        break;
                    }

                    const AbstractLight *light = m_scene.GetLightPtr(light_id);
                    float directPdfA;
                    vec3 contrib = light->GetRadiance(m_scene.m_sceneSphere,
                        ray.dir, hitPoint, &directPdfA);
                    if (isZero(contrib))
                    {
                        break;
                    }

                    float misWeight = 1.0f;
                    if (pathLength > 1 && !lastSpecular) // no passive evaluation for primary rays, also no specular
                        // on last vertex because it's fully accounted for by passive sampling (misWeight = 1)
                    {
                        const float directPdfW = pdfAtoW(directPdfA, hit.distance(),
                            bsdf.CosThetaFix()); // the pdf of having sampled the passive hit direction
                        // in direct lighting, so using cosine on the light
                        misWeight = Mis2(lastPdfW, directPdfW * lightPickProb); // multi-sample MIS,
                        // since NEE on the last vertex has taken care of direct light sampling
                    }

                    color += pathWeight * misWeight * contrib;
                    break;
                }
#pragma endregion passive light hit, terminates path

                if (pathLength >= m_maxPathLength)
                {
                    break;
                }

                if (bsdf.ContinuationProb() == 0) // at this point NEE and passive lighting are both evaluated
                {
                    break;
                }

#pragma region
                // next event estimation
                if (!bsdf.IsDelta() && pathLength + 1 >= m_minPathLength) // delta bsdf is accounted for entirely by passive lighting
                {
                    int lightID = int(randf() * lightCount);
                    const AbstractLight *light = m_scene.GetLightPtr(lightID);

                    vec3 directionToLight;
                    float distance, directPdfW;
                    vec3 radiance = light->Illuminate(m_scene.m_sceneSphere, hitPoint,
                        vec2(randf(), randf(), 0), directionToLight, distance, directPdfW);

                    if (!isZero(radiance))
                    {
                        float bsdfPdfW, cosThetaOut;
                        const vec3 factor = bsdf.Evaluate(
                            directionToLight, cosThetaOut, &bsdfPdfW);

                        if (!isZero(factor))
                        {
                            float weight = 1.0f;
                            if (!light->IsDelta()) // if light is delta, does not use MIS,
                                // since there's no chance for passive light hit
                            {
                                const float contProb = bsdf.ContinuationProb();
                                bsdfPdfW *= contProb; // Russian Roulette
                                weight = Mis2(directPdfW * lightPickProb, bsdfPdfW);
                            }

                            vec3 contrib = (weight * cosThetaOut / (lightPickProb * directPdfW)) *
                                (radiance * factor);

                            if (!m_scene.Occluded(hitPoint, directionToLight, distance)) // visibility test
                            {
                                color += pathWeight * contrib;
                            }
                        }
                    }
                }
#pragma endregion direct lighting / next event estimation

#pragma region
                // continue random walk
                {
                    vec3 rndTriplet = vec3(randf(), randf(), randf());
                    float pdf, cosThetaOut;
                    uint32_t  sampledEvent;

                    vec3 factor = bsdf.Sample(rndTriplet, ray.dir,
                        pdf, cosThetaOut, &sampledEvent);

                    if (isZero(factor))
                    {
                        break;
                    }

                    // Russian roulette
                    const float contProb = bsdf.ContinuationProb();

                    lastSpecular = (sampledEvent & BSDF<true>::kSpecular) != 0;
                    lastPdfW = pdf * contProb;

                    if (contProb < 1.0f)
                    {
                        if (randf() > contProb)
                        {
                            break;
                        }
                        pdf *= contProb;
                    }

                    // pdf == lastPdfW at this point
                    pathWeight = pathWeight * (factor * (cosThetaOut / pdf));

                    // We offset ray origin instead of setting tmin due to numeric
                    // issues in ray-sphere intersection. The isect.dist has to be
                    // extended by this EPS_RAY after hitpoint is determined
                    ray.orig = hitPoint + NUM_EPS_RAY * ray.dir;
                }
#pragma endregion continue random walk
            }
            m_framebuffer.accumulatePixel(x, y, color);
        }

        m_iterations++;
    }
};

