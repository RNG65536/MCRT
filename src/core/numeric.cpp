#include <mutex>
#include <random>
#include "constants.h"
#include "numeric.h"
#include "vectors.h"

static std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);
static std::default_random_engine rng;
static std::mutex _mutex;

float randf()
{
//     std::unique_lock<std::mutex> lock(_mutex); // extremely slow

    //     return (float(rand()) + float(rand()) * float(RAND_MAX)) / (float(RAND_MAX) * float(RAND_MAX));

    return uniform_dist(rng);
}

float f_min(float a, float b)
{
    return a < b ? a : b;
}

float f_max(float a, float b)
{
    return a > b ? a : b;
}

int i_min(int a, int b)
{
    return a < b ? a : b;
}

int i_max(int a, int b)
{
    return a > b ? a : b;
}

float clamp(float x, float a, float b)
{
    return x<a ? a : x>b ? b : x;
}

// float log2(float x)
// {
//     return std::log(x) * 1.4426950408889634f;
// }

float fastlog2(float x)
{
    union {
        float f;
        uint32_t i;
    } vx = { x };
    union {
        uint32_t i;
        float f;
    } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;

    return y - 124.22551499f
        - 1.498030302f * mx.f
        - 1.72587999f / (0.3520887068f + mx.f);
}

float fasterlog2(float x)
{
    union {
        float f;
        uint32_t i;
    } vx = { x };
    float y = vx.i;
    y *= 1.1920928955078125e-7f;
    return y - 126.94269504f;
}

float sq(float x)
{
    return x * x;
}

float linearRatio(float x, float a, float b)
{
    return (x - a) / (b - a);
}

float toRadians(float deg)
{
    return deg / 180.0f * NUM_PI;
}

float toDegrees(float rad)
{
    return rad / NUM_PI * 180.0f;
}

float luminance(const vec3& rgb)
{
    return 0.212671f * rgb.x +
        0.715160f * rgb.y +
        0.072169f * rgb.z;

// return 0.299 * c.x + 0.587 * c.y + 0.114 * c.z;
}


void TKahanAdder::add(const double b)
{
    y = b - carry;
    const double t = sum + y;
    carry = (t - sum) - y;
    sum = t;
}

TKahanAdder::TKahanAdder(const double b /*= 0.0*/)
{
    sum = b;
    carry = 0.0;
    y = 0.0;
}

// int DiscreteSampler::sample(float rnd)
// {
//     int count;
// 
//     // TODO: use binary search
// #if 0
//     for (int n = 0; n < m_cdf.size(); n++) {
//         count = n;
//         if (cdf[n] > u) {
//             break;
//         }
//     }
// #else
//         {
//             int begin = 0;
//             int end = m_cdf.size() - 1;
//             while (end > begin) {
//                 int mid = begin + (end - begin) / 2;
//                 float c = m_cdf[mid];
//                 if (c >= u) {
//                     end = mid;
//                 }
//                 else {
//                     begin = mid + 1;
//                 }
//             }
//             count = begin;
//         }
// #endif
// 
//         return count;
// }
// 
// void DiscreteSampler::update(std::vector<float> pdf)
// {
//     m_cdf.clear();
//     m_cdf.resize(pdf.size());
// 
//     float acc = 0;
//     for (int n = 0; n < pdf.size(); n++)
//     {
//         acc += pdf[n];
//         m_cdf[n] = acc;
//     }
// }
// 
// DiscreteSampler::DiscreteSampler()
// {
// 
// }
