#pragma once

// utilities

class vec3;

// not thread-safe
float randf();
float randfFast();

float f_min(float a, float b);
float f_max(float a, float b);
int   i_min(int a, int b);
int   i_max(int a, int b);
float clampf(float x, float a, float b);

// square
float sq(float x);

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
