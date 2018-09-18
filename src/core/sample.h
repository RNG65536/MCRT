#pragma once

#include <random>
#include <stack>
#include <vector>
#include "vectors.h"

class RandomNumberGenerator
{
    std::uniform_real_distribution<float> _dist;
    std::default_random_engine            _rng;

public:
    RandomNumberGenerator();
    float randf();
};

class DiscreteSampler
{
    std::vector<float> m_cdf;
    std::vector<float> m_pdf;

public:
    DiscreteSampler();
    DiscreteSampler(std::vector<float>& data);
    void update(std::vector<float>& data);

    int   sample();
    int   sample(float r);
    float pdf(int idx) const;
};

class Sampler
{
protected:
    RandomNumberGenerator _random;

public:
    virtual ~Sampler();
    virtual float next() = 0;
    float         randf();
};

class StandardSampler : public Sampler
{
public:
    StandardSampler();
    float next();
};

struct PrimarySample
{
    int   modify_time = 0;
    float value = 0.0f;

    PrimarySample();
    PrimarySample(float val);
};

struct PathSample
{
    int   x, y;
    vec3  F;
    float weight;
};

class MetropolisSampler : public Sampler
{
    float mutate(float x);

    int  used_rands;
    int  global_time;
    int  large_step_time;
    bool large_step;

    std::vector<PrimarySample> primary_samples;
    std::stack<PrimarySample>  backup;

public:
    MetropolisSampler();
    void  reset(bool large);
    float next();
    void  accept();
    void  reject();
};

// return in tangent, bitangent, normal order
vec3 sampleHemisphereCosine(float x, float y);
vec3 sampleHemisphereCosinePower(float g, float x, float y);
vec3 sampleHemisphereUniform(float x, float y, float* pdf = nullptr);
vec3 sampleSphereUniform(float x, float y, float* pdf = nullptr);
vec3 sampleSolidAngleUniform(float  x,
                             float  y,
                             float  cos_theta_min,
                             float* pdf = nullptr);

// sample triangle, returns barycentric coordinates
vec2 sampleTriangleUniform(float x, float y);

/// Sample direction in the upper hemisphere with cosine-proportional pdf
/** The returned PDF is with respect to solid angle measure */
vec3  sampleHemisphereCosineW(float x, float y, float* pdfW = nullptr);
float sampleHemisphereCosinePdfW(const vec3& normal, const vec3& direction);

// Cosine lobe hemisphere sampling
vec3  sampleHemisphereCosinePowerW(float  x,
                                   float  y,
                                   float  power,
                                   float* pdfW = nullptr);
float sampleHemisphereCosinePowerPdfW(const vec3& normal,
                                      const vec3& direction,
                                      float       power);

// Disk sampling, map from square to disk
vec2  sampleConcentricDiscA(float x, float y);
float sampleConcentricDiscPdfA();

// Sphere sampling
// vec3 sampleSphereUniformW(float x, float y, float *pdfW);
float sampleSphereUniformPdf();
