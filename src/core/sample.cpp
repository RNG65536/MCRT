#include <cassert>
#include <functional>
#include <numeric>
#include <algorithm>
#include "constants.h"
#include "numeric.h"
#include "sample.h"

// TODO : use pi-related constants

int DiscreteSampler::sample(float r)
{
#if 0  // wrong ??
    int a = 0, b = cdf.size() - 1;
    int debug_count = 0;
    while (b - a > 1) {
        int n = (a + b) / 2;
        if (cdf[n] > r) {
            b = n;
        }
        else {
            a = n;
        }
        debug_count++;
    }
    return b;
#elif 0
    int count;
    for (int n = 0; n < m_cdf.size(); n++)
    {
        count = n;
        if (m_cdf[n] > u)
        {
            break;
        }
    }
    return count;
#else
    int begin = 0;
    int end = m_cdf.size() - 1;
    while (end > begin)
    {
        int   mid = begin + (end - begin) / 2;
        float c = m_cdf[mid];
        if (c >= r)
        {
            end = mid;
        }
        else
        {
            begin = mid + 1;
        }
    }
    return begin;
#endif
}

int DiscreteSampler::sample()
{
    // a : begin
    // b : end
    float r = randf();
    int   begin = 0;
    int   end = m_cdf.size() - 1;
    int   debug_count = 0;

    //         while (cdf[b] > r && r >= cdf[a]) {
    while (end > begin)
    {
        int mid = begin + (end - begin) / 2;
        //             printf("%d\n", n);
        if (m_cdf[mid] >= r)
        {
            end = mid;
        }
        else
        {
            begin = mid + 1;
        }
        debug_count++;
    }

    //         printf("took %d iters, a(%d), b(%d), b-a(%d)\n\tcdf[a](%f) <=
    //         r(%f) < cdf[b](%f)\n",
    //             debug_count, a, b, b-a, cdf[a], r, cdf[b]);
    return end;
}

float DiscreteSampler::pdf(int idx) const
{
    return m_pdf[idx];
}

DiscreteSampler::DiscreteSampler(std::vector<float>& data)
{
    update(data);
}

DiscreteSampler::DiscreteSampler()
{
}

void DiscreteSampler::update(std::vector<float>& data)
{
    m_pdf = data;  // pdf is not needed for sampling purpose
    m_cdf.resize(data.size());
    float norm = std::accumulate(m_pdf.begin(), m_pdf.end(), 0.0f);
    std::transform(m_pdf.begin(),
                   m_pdf.end(),
                   m_pdf.begin(),
                   std::bind2nd(std::multiplies<float>(), 1.0f / norm));
    float sum = 0;
    for (int n = 0; n < m_pdf.size(); n++)
    {
        sum += m_pdf[n];
        m_cdf[n] = sum;
    }
}

// see tungsten\src\core\sampling\distribution1d
#if 0  // WITH BUGS
int sample() {
    return sample(randf());
}
int sample(float r) {
    return std::distance(cdf.begin(), std::upper_bound(cdf.begin(), cdf.end(), r)) - 1;
}
#endif

//////////////////////////////////////////////////////////////////////////

vec3 sampleHemisphereCosine(float x, float y)
{
    x *= 2 * NUM_PI;
    float r2s = std::sqrt(y);
    return vec3(cos(x) * r2s, sin(x) * r2s, std::sqrt(1 - y));
}

vec3 sampleHemisphereCosinePower(float g, float x, float y)
{
    x *= 2 * NUM_PI;
    float r2p = pow(y, 1.0 / (g + 1.0));
    float s = sin(x), c = cos(x), t = std::sqrt(1.0 - r2p * r2p);
    return vec3(s * t, c * t, r2p);
}

vec3 sampleHemisphereUniform(float x, float y, float* pdf)
{
    y *= 2 * NUM_PI;
    float r = std::sqrt(1 - x * x);

    if (pdf)
    {
        *pdf = 0.5f / NUM_PI;
    }

    return vec3(cos(y) * r, sin(y) * r, x);
}

vec3 sampleSphereUniform(float x, float y, float* pdf)
{
    y *= 2 * NUM_PI;
    float cos_theta = 1.0f - 2.0f * x;
    assert(cos_theta >= -1 && cos_theta <= 1);
    float sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);

    if (pdf)
    {
        *pdf = 0.25f / NUM_PI;
    }

    return vec3(cosf(y) * sin_theta, sinf(y) * sin_theta, cos_theta);
}

vec3 sampleSolidAngleUniform(float x, float y, float cos_theta_min, float* pdf)
{
    y *= 2 * NUM_PI;
    float cos_theta = (1 - x) + x * cos_theta_min;  // uniform sample cos_theta
    assert(cos_theta >= -1 && cos_theta <= 1);
    float sin_theta = std::sqrt(1 - cos_theta * cos_theta);

    if (pdf)
    {
        *pdf = 0.5 / (NUM_PI * (1 - cos_theta_min));
    }

    return vec3(cos(y) * sin_theta, sin(y) * sin_theta, cos_theta);
}

float sampleSphereUniformPdf()
{
    return 0.25f / NUM_PI;
}

// vec3 sampleSphereUniformW(float x, float y, float *pdfW)
// {
//     const float term1 = 2.0f * NUM_PI * x;
//     const float term2 = 2.0f * std::sqrt(y - y * y);
//
//     const vec3 ret(
//         std::cos(term1) * term2,
//         std::sin(term1) * term2,
//         1.0f - 2.0f * y);
//
//     if (pdfW)
//     {
//         *pdfW = 0.25f / NUM_PI;
//     }
//
//     return ret;
// }

//////////////////////////////////////////////////////////////////////////

vec2 sampleTriangleUniform(float x, float y)
{
    float t = std::sqrt(x);
    return vec3(1.0f - t, y * t, 0.0f);
}

// vec2 sampleTriangleUniform(float x, float y)
// {
//     vec2 ret;
//     float u_sqrt = sqrt(x);
//     ret.x = 1 - u_sqrt;
//     ret.y = (1 - y) * u_sqrt;
//     return ret;
// }

vec3 sampleHemisphereCosineW(float x, float y, float* pdfW)
{
    const float term1 = 2.0f * NUM_PI * x;
    const float term2 = std::sqrt(1.0f - y);

    const vec3 ret(
        std::cos(term1) * term2, std::sin(term1) * term2, std::sqrt(y));

    if (pdfW)
    {
        *pdfW = ret.z / NUM_PI;
    }

    return ret;
}

float sampleHemisphereCosinePdfW(const vec3& normal, const vec3& direction)
{
    return f_max(0.0f, dot(normal, direction)) / NUM_PI;
}

vec2 sampleConcentricDiscA(float x, float y)
{
    float phi, r;

    float a = 2 * x - 1; /* (a,b) is now on [-1,1]^2 */
    float b = 2 * y - 1;

    if (a > -b) /* region 1 or 2 */
    {
        if (a > b) /* region 1, also |a| > |b| */
        {
            r = a;
            phi = (NUM_PI / 4.0f) * (b / a);
        }
        else /* region 2, also |b| > |a| */
        {
            r = b;
            phi = (NUM_PI / 4.0f) * (2.0f - (a / b));
        }
    }
    else /* region 3 or 4 */
    {
        if (a < b) /* region 3, also |a| >= |b|, a != 0 */
        {
            r = -a;
            phi = (NUM_PI / 4.0f) * (4.0f + (b / a));
        }
        else /* region 4, |b| >= |a|, but a==0 and b==0 could occur. */
        {
            r = -b;

            if (b != 0)
            {
                phi = (NUM_PI / 4.0f) * (6.0f - (a / b));
            }
            else
            {
                phi = 0;
            }
        }
    }

    vec2 res;
    res.x = r * std::cos(phi);
    res.y = r * std::sin(phi);
    return res;
}

float sampleConcentricDiscPdfA()
{
    return 1 / NUM_PI;
}

vec3 sampleHemisphereCosinePowerW(float x, float y, float power, float* pdfW)
{
    const float term1 = 2.0f * NUM_PI * x;
    const float term2 = std::pow(y, 1.0f / (power + 1.0f));
    const float term3 = std::sqrt(1.0f - term2 * term2);

    if (pdfW)
    {
        *pdfW = (power + 1.0f) * std::pow(term2, power) * (0.5f / NUM_PI);
    }

    return vec3(std::cos(term1) * term3, std::sin(term1) * term3, term2);
}

float sampleHemisphereCosinePowerPdfW(const vec3& normal,
                                      const vec3& direction,
                                      float       power)
{
    const float cosTheta = f_max(0.0f, dot(normal, direction));

    return (power + 1.0f) * std::pow(cosTheta, power) * (0.5f / NUM_PI);
}

Sampler::~Sampler()
{
}

float Sampler::randf()
{
    return _random.randf();
}

float StandardSampler::next()
{
    return _random.randf();
}

StandardSampler::StandardSampler()
{
}

void MetropolisSampler::reject()
{
    int idx = used_rands - 1;
    while (!backup.empty())
    {
        primary_samples[idx--] = backup.top();
        backup.pop();
    }
}

void MetropolisSampler::accept()
{
    if (large_step)
    {
        large_step_time = global_time;
    }

    global_time++;

    std::stack<PrimarySample> empty;
    backup.swap(empty);
}

float MetropolisSampler::next()
{
    if (primary_samples.size() <= used_rands)
    {
        int old_size = primary_samples.size();
        primary_samples.resize(primary_samples.size() * 2);
    }

    PrimarySample& s = primary_samples[used_rands];
    used_rands++;

    if (s.modify_time < global_time)
    {
        if (large_step)
        {
            backup.push(s);
            s.modify_time = global_time;
            s.value = _random.randf();
        }
        else
        {
            if (s.modify_time < large_step_time)
            {
                s.modify_time = large_step_time;
                s.value = _random.randf();
            }

            while (s.modify_time < global_time - 1)
            {
                s.value = mutate(s.value);
                s.modify_time++;
            }
            backup.push(s);
            s.value = mutate(s.value);
            s.modify_time = global_time;
        }
    }

    return s.value;
}

void MetropolisSampler::reset(bool large)
{
    used_rands = 0;
    large_step = large;
}

MetropolisSampler::MetropolisSampler()
{
    global_time = 0;
    large_step_time = 0;
    used_rands = 0;
    large_step = false;

    primary_samples.resize(128);
}

float MetropolisSampler::mutate(float x)
{
    float r = _random.randf();
    float s1 = 1.0 / 512.0, s2 = 1.0 / 16.0;
    //     float s1 = 1.0 / 1024.0, s2 = 1.0 / 32.0;
    float dx = s1 / (s1 / s2 + std::abs(2.0 * r - 1.0)) - s1 / (s1 / s2 + 1.0);
    if (r < 0.5)
    {
        float x1 = x + dx;
        return x1 < 1.0 ? x1 : x1 - 1.0;
    }
    else
    {
        float x1 = x - dx;
        return x1 < 0.0 ? x1 + 1.0 : x1;
    }
}

float RandomNumberGenerator::randf()
{
    return _dist(_rng);
}

RandomNumberGenerator::RandomNumberGenerator() : _dist(0.0f, 1.0f)
{
    union {
        float        f;
        unsigned int i;
    } u;
    u.f = ::randf();
    unsigned int seed = u.i;
    _rng.seed(seed);
    _dist.reset();
}

PrimarySample::PrimarySample(float val)
{
    //     value = randf();
}

PrimarySample::PrimarySample()
{
    //     value = randf();
}
