#ifndef numeric_h__
#define numeric_h__

// utilities

class vec3;

// care not to use with multiple threads
float randf();

float f_min(float a, float b);
float f_max(float a, float b);
int i_min(int a, int b);
int i_max(int a, int b);
float clamp(float x, float a, float b);
float sq(float x);
float linearRatio(float x, float a, float b);

// float log2(float x);
float fastlog2(float x);
float fasterlog2(float x);

// compensated summation
struct TKahanAdder
{
    double sum, carry, y;
    TKahanAdder(const double b = 0.0);
    void add(const double b);
};

float toRadians(float deg);
float toDegrees(float rad);

// sRGB luminance
float luminance(const vec3& rgb);

#endif // numeric_h__
