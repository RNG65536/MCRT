#pragma once

extern const float NUM_INFINITY;
extern const float NUM_PI;

extern const float NUM_EPS;
extern const float NUM_EPS_DET;
extern const float NUM_EPS_RAY;
extern const float NUM_EPS_COSINE;

// $TODO : rename to avoid conflict
#ifdef M_PI
#undef M_PI
#endif
#ifdef M_PI_2
#undef M_PI_2
#endif
#ifdef M_1_PI
#undef M_1_PI
#endif
#ifdef M_4PI
#undef M_4PI
#endif
#ifdef M_1_4PI
#undef M_1_4PI
#endif
#ifdef M_1_4PIPI
#undef M_1_4PIPI
#endif

extern const float M_PI;
extern const float M_PI_2;
extern const float M_1_PI;
extern const float M_4PI;
extern const float M_1_4PI;
extern const float M_1_4PIPI;
