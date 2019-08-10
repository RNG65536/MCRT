#pragma once
#include <cstdint>

// custom vector and operators

class vec3
{
public:
    friend vec3 operator*(float a, const vec3& b);

    vec3();
    explicit vec3(float a);
    vec3(float x, float y, float z);
    vec3         operator+(const float& v) const;
    vec3         operator+(const vec3& v) const;
    vec3         operator-(const float& v) const;
    vec3         operator-(const vec3& v) const;
    vec3&        operator+=(const vec3& v);
    vec3&        operator-=(const vec3& v);
    vec3         operator*(float f) const;
    vec3         operator*(const vec3& v) const;
    vec3         operator/(const vec3& v) const;
    vec3&        operator*=(float f);
    vec3&        operator*=(const vec3& v);
    vec3         operator-() const;
    vec3         operator/(float f) const;
    vec3&        operator/=(float f);
    const float& operator[](int i) const;
    float&       operator[](int i);
    vec3         normalize() const;
    float        lengthSquared() const;
    float        length() const;
    bool         operator==(const vec3& v) const;
    bool         operator!=(const vec3& v) const;

    bool isNormalized() const;

public:
    float x, y, z;
};

typedef vec3 float3;

void debugPrint(const vec3& a);

vec3  operator+(const float a, const vec3& v);
vec3  operator-(const float a, const vec3& v);
vec3  mult(const vec3& v1, const vec3& v2);
vec3  div(const vec3& v1, const vec3& v2);
float dot(const vec3& v1, const vec3& v2);
vec3  cross(const vec3& v1, const vec3& v2);
float length(const vec3& v);
float lengthSquared(const vec3& v);
vec3  normalize(const vec3& v);

vec3 exp(vec3 a);
vec3 sqrt(vec3 a);
vec3 sq(vec3 a);

vec3 f_min(const vec3& a, const vec3& b);
vec3 f_max(const vec3& a, const vec3& b);

float maxComponent(const vec3& v);
float minComponent(const vec3& v);
bool  isZero(const vec3& v);

class int3
{
public:
    int x, y, z;
    int3();
    int3(int a, int b, int c);
};

class byte3
{
public:
    int8_t x, y, z;
    byte3();
    byte3(int8_t x, int8_t y, int8_t z);
};

class vec3d
{
public:
    double x, y, z;

    vec3d();
    explicit vec3d(double a);
    vec3d(double x, double y, double z);
    explicit vec3d(const vec3& v);
    vec3  toFloat() const;
    vec3d operator+(const vec3d& v) const;
    vec3d operator-(const vec3d& v) const;
    vec3d operator*(double f) const;
};

double length(const vec3d& v);
vec3d  normalize(const vec3d& v);
double dot(const vec3d& v1, const vec3d& v2);

typedef vec3 vec2;  // TODO : implement vec2
typedef int3 int2;
