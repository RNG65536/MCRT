#include <cmath>
#include "aabb.h"
#include "constants.h"
#include "numeric.h"
#include "triangle.h"

TriangleObject::TriangleObject()
{
}

TriangleObject::TriangleObject(vec3 a, vec3 b, vec3 c)
    : a(a),
      b(b),
      c(c),
      n(normalize(cross((a - c), (b - c)))),
      na(n),
      nb(n),
      nc(n)
{
}

TriangleObject::TriangleObject(
    vec3 a, vec3 b, vec3 c, vec3 na, vec3 nb, vec3 nc)
    : a(a),
      b(b),
      c(c),
      n(normalize(cross((a - c), (b - c)))),
      na(na),
      nb(nb),
      nc(nc)
{
}

vec3 TriangleObject::normal(float u, float v) const
{
    return normalize((1.0f - u - v) * na + u * nb + v * nc);
}

vec3 TriangleObject::point(float u, float v) const
{
    return (1.0f - u - v) * a + u * b + v * c;
}

vec3 TriangleObject::samplePoint(float& u, float& v, float r0, float r1) const
{
    r0 = sqrt(r0);
    u = 1.0f - r0;
    v = r1 * r0;  // rand * (1 - u)
    return a * (1.0f - u - v) + b * u + c * v;
}

vec3 TriangleObject::samplePoint(float& u, float& v) const
{
    float r0 = sqrt(randf());
    u = 1.0f - r0;
    v = randf() * r0;  // rand * (1 - u)
    return a * (1.0f - u - v) + b * u + c * v;
}

vec3 TriangleObject::samplePoint() const
{
    float r0 = sqrt(randf());
    float u = 1.0f - r0;
    float v = randf() * r0;  // rand * (1 - u)
    return a * (1.0f - u - v) + b * u + c * v;
}

vec3 TriangleObject::normal() const
{
    return n;
}

float TriangleObject::area() const
{
    return length(cross((b - a), (c - a))) * 0.5f;
}

vec3 TriangleObject::boxCenter() const
{
    return (a + b + c) * (1.0f / 3.0f);
}

void TriangleObject::setMaterialID(int id)
{
    m_matID = id;
}

int TriangleObject::materialID() const
{
    return m_matID;
}

vec3 TriangleObject::consistentNormal(
    const vec3& wi, float u, float v, vec3* refl, float* cos_wi) const
{
    const vec3 np = normal(u, v);
    //     if (dot(wi, np) > 0)
    //     {
    //         return np;
    //     }

    // for alpha = 0, consistent normal becomes shading normal
    float alpha = ((1.0f - u - v) * m_vertex_alpha.x + u * m_vertex_alpha.y +
                   v * m_vertex_alpha.z);

    float q = sq(1.0f - (2.0f / NUM_PI) * alpha) /
              (1.0f + 2.0f * (1.0f - 2.0f / NUM_PI) * alpha);
    float      b = dot(wi, np);
    float      g = 1.0f + q * (b - 1.0f);
    float      rho = sqrt(q * (1.0f + g) / (1.0f + b));
    const vec3 r = (g + rho * b) * np - rho * wi;

    // reflected by nc
    if (refl)
    {
        *refl = r;
    }

    // dot(wi, nc)
    if (cos_wi)
    {
        *cos_wi = sqrt((1.0f + (g + rho * b) * b - rho) * 0.5f);
    }

    return normalize(wi + r);
}

void TriangleObject::setVertexAlpha(const vec3& alpha)
{
    m_vertex_alpha = alpha;
}

AABB TriangleObject::boundingBox() const
{
    AABB ret(f_min(a, b), f_max(a, b));
    ret.enclose(c);
    return ret;
}
