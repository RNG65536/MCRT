#include <limits>
#include "constants.h"

static const double const_pi = 3.1415926535897932384626433832795;

const float NUM_INFINITY = std::numeric_limits<float>::infinity();
const float NUM_PI = static_cast<float>(const_pi);

const float NUM_EPS = 1e-10f;
//const float NUM_EPS_DET = 1e-10f;
const float NUM_EPS_DET = 1e-6f;
// const float M_EPS_RAY = 1e-5f; //++++++++++++++++++++++ better than 1e-6??
const float NUM_EPS_RAY = 1e-4f;  // works for BDPT
// const float NUM_EPS_RAY = 1e-3f; // works for BDPT
const float NUM_EPS_COSINE = 1e-6f;
// #define eps_ray 1e-4f
// #define eps_ray 1e-6f
// #define eps 1e-6f // or 1e-4, for all

const float M_PI = static_cast<float>(const_pi);
const float M_PI_2 = static_cast<float>(const_pi / 2.0);
const float M_1_PI = static_cast<float>(1.0 / const_pi);
const float M_2PI = static_cast<float>(2.0 * const_pi);
const float M_4PI = static_cast<float>(4.0 * const_pi);
const float M_1_4PI = static_cast<float>(1.0 / (4.0 * const_pi));
const float M_1_4PIPI = static_cast<float>(1.0 / (4.0 * const_pi * const_pi));
