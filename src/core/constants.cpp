#include <limits>
#include "constants.h"

const float NUM_INFINITY = std::numeric_limits<float>::infinity();
const float NUM_PI = 4.0 * atan(1.0);

const float NUM_EPS = 1e-10f;
const float NUM_EPS_DET = 1e-10f;
// const float M_EPS_RAY = 1e-5f; //++++++++++++++++++++++ better than 1e-6??
const float NUM_EPS_RAY = 1e-4f; // works for BDPT
// const float NUM_EPS_RAY = 1e-3f; // works for BDPT
const float NUM_EPS_COSINE = 1e-6f;
// #define eps_ray 1e-4f
// #define eps_ray 1e-6f
// #define eps 1e-6f // or 1e-4, for all

const float M_PI = 3.1415926535897932384626433832795;
const float M_1_PI = 1.0 / M_PI;
const float M_4PI = 4.0 * M_PI;
const float M_1_4PI = 1.0 / M_4PI;
const float M_1_4PIPI = M_1_4PI * M_1_PI;

