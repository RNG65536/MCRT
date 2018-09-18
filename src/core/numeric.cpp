#include <mutex>
#include <random>
#include "constants.h"
#include "numeric.h"
#include "vectors.h"

static std::uniform_real_distribution<float> uniform_dist(0.0f, 1.0f);
static std::default_random_engine            rng;
static std::mutex                            _mutex;

float randf()
{
    //     std::unique_lock<std::mutex> lock(_mutex); // extremely slow
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
    return x < a ? a : x > b ? b : x;
}

float sq(float x)
{
    return x * x;
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
    // return 0.212671f * rgb.x + 0.715160f * rgb.y + 0.072169f * rgb.z;
    return 0.299 * rgb.x + 0.587 * rgb.y + 0.114 * rgb.z;
}

void TKahanAdder::add(const double b)
{
    y = b - carry;
    const double t = sum + y;
    carry = (t - sum) - y;
    sum = t;
}

TKahanAdder::TKahanAdder(const double b)
{
    sum = b;
    carry = 0.0;
    y = 0.0;
}
