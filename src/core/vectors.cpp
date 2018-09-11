#include <cmath>
#include "vectors.h"
#include "numeric.h"

vec3::vec3(float a) : x(a), y(a), z(a)
{

}

vec3::vec3() : x(0), y(0), z(0)
{

}

vec3::vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {    }

vec3 vec3::operator+(const float &v) const{ return vec3(x + v, y + v, z + v); }
vec3 vec3::operator-(const float &v) const{ return vec3(x - v, y - v, z - v); }
vec3 vec3::operator+(const vec3 &v) const{ return vec3(x + v.x, y + v.y, z + v.z); }
vec3 vec3::operator-(const vec3 &v) const{ return vec3(x - v.x, y - v.y, z - v.z); }
vec3& vec3::operator+=(const vec3 &v) { x += v.x; y += v.y; z += v.z;        return *this; }
vec3& vec3::operator-=(const vec3 &v) { x -= v.x; y -= v.y; z -= v.z;        return *this; }
vec3 vec3::operator*(float f) const{ return vec3(x*f, y*f, z*f); }

vec3 vec3::operator*(const vec3& v) const
{
    return vec3(x * v.x, y * v.y, z * v.z);
}

vec3& vec3::operator*=(float f) { x *= f;        y *= f;        z *= f;        return *this; }
vec3& vec3::operator*=(const vec3& v) { x *= v.x;        y *= v.y;        z *= v.z;        return *this; }
vec3 vec3::operator-() const{ return vec3(-x, -y, -z); }
vec3 vec3::operator/(float f) const{ float r = 1.0f / f;        return vec3(x * r, y * r, z * r); }

vec3 vec3::operator/(const vec3& v) const
{
    return vec3(x / v.x, y / v.y, z / v.z);
}

vec3& vec3::operator/=(float f) { float r = 1.0f / f;        x *= r;        y *= r;        z *= r;        return *this; }
// float vec3::operator % (const vec3 &b) const{ return x*b.x + y*b.y + z*b.z; }
float vec3::operator[](int i) const { return (&x)[i]; }
float& vec3::operator[](int i) { return (&x)[i]; }
// vec3 vec3::norm() _CONST_ {        real len = length() + 1e-9;        return vec3(x / len, y / len, z / len);    }
vec3 vec3::norm() const{ float len = length()/* + 1e-30*/;        return vec3(x / len, y / len, z / len); }
// vec3 vec3::norm() _CONST_ {        real inv_len = 1.0f / length();        return vec3(x *inv_len, y *inv_len, z *inv_len);    }
// vec3 vec3::norm() _CONST_ {        real len = length();   if (length() < 1e-9) return vec3(1,0,0);    return vec3(x / len, y / len, z / len);    }
vec3 vec3::normalize() const{ float invLen = 1.0f / length();        return vec3(x * invLen, y * invLen, z * invLen); }
float vec3::lengthSquared() const{ return x*x + y*y + z*z; }
float vec3::length() const{ return std::sqrt(lengthSquared()); }
bool vec3::operator==(const vec3 &v) const{ return (v.x == x && v.y == y && v.z == z); }
bool vec3::operator!=(const vec3 &v) const{ return !operator==(v); }
// vec3 vec3::mult(const vec3 &v2) const{ return vec3(x*v2.x, y*v2.y, z*v2.z); }
// vec3 vec3::operator ^ (const vec3 &b) const{ return vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }//cross product

vec3 operator*(float a, const vec3 &b){ return b * a; }

vec3 operator + (const float a, const vec3& v){ return v + a; }
vec3 operator - (const float a, const vec3& v){ return -(v - a); }
vec3 mult(const vec3 &v1, const vec3 &v2) { return vec3(v1.x*v2.x, v1.y*v2.y, v1.z*v2.z); }
vec3 div(const vec3 &v1, const vec3 &v2){ return vec3(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z); }
float dot(const vec3 &v1, const vec3 &v2) { return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z; }

double dot(const vec3d &v1, const vec3d &v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

vec3 cross(const vec3 &v1, const vec3 &v2) {
    return vec3(
        (v1.y * v2.z) - (v1.z * v2.y),
        (v1.z * v2.x) - (v1.x * v2.z),
        (v1.x * v2.y) - (v1.y * v2.x)
        );
}

vec3 exp(vec3 a)
{
    return vec3(expf(a.x), expf(a.y), expf(a.z));
}

vec3 sqrt(vec3 a)
{
    return vec3(std::sqrtf(a.x), std::sqrtf(a.y), std::sqrtf(a.z));
}

vec3 sq(vec3 a)
{
    return vec3(a.x*a.x, a.y*a.y, a.z*a.z);
}

vec3 f_min(const vec3& a, const vec3& b)
{
    return vec3(f_min(a.x, b.x), f_min(a.y, b.y), f_min(a.z, b.z));
}

vec3 f_max(const vec3& a, const vec3& b)
{
    return vec3(f_max(a.x, b.x), f_max(a.y, b.y), f_max(a.z, b.z));
}

float length(const vec3& v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

double length(const vec3d &v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

float lengthSquared(const vec3& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

vec3 normalize(const vec3& v)
{
    float invLen = 1.0f / length(v);
    return vec3(v.x * invLen, v.y * invLen, v.z * invLen);
}

vec3d normalize(const vec3d &v)
{
    double invLen = 1.0 / length(v);
    return vec3d(v.x * invLen, v.y * invLen, v.z * invLen);
}

float maxComponent(const vec3& v)
{
    return f_max(f_max(v.x, v.y), v.z);
}

float minComponent(const vec3& v)
{
    return f_min(f_min(v.x, v.y), v.z);
}

bool isZero(const vec3& v)
{
    return v.x <= 0 && v.y <= 0 && v.z <= 0;
}

int3::int3(int a, int b, int c) : x(a), y(b), z(c)
{

}

int3::int3() : x(0), y(0), z(0)
{

}

vec3d vec3d::operator*(double f) const
{
    return vec3d(x * f, y * f, z * f);
}

vec3d vec3d::operator-(const vec3d &v) const
{
    return vec3d(x - v.x, y - v.y, z - v.z);
}

vec3d vec3d::operator+(const vec3d &v) const
{
    return vec3d(x + v.x, y + v.y, z + v.z);
}

vec3 vec3d::toFloat() const
{
    return vec3(x, y, z);
}

vec3d::vec3d(const vec3& v) : x(v.x), y(v.y), z(v.z)
{

}

vec3d::vec3d(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
{

}

vec3d::vec3d(double a) : x(a), y(a), z(a)
{

}

vec3d::vec3d() : x(0), y(0), z(0)
{

}
