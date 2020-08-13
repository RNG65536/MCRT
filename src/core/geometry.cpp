#include <cassert>
#include "geometry.h"
#include "lightpath.h"
#include "material.h"
#include "numeric.h"
#include "scene.h"
#include "triangle.h"

float cosineFactor(const vec3 wi, const vec3& nl)
{
    return f_max(dot(wi, nl), 0.0f);
}

float fresnelReflectance(float cos_theta_i, float ior)
{
    if (ior < 0)
    {
        return 1;
    }

    float eta;  // eta_i / eta_t
    if (cos_theta_i < 0)
    {
        cos_theta_i = -cos_theta_i;
        eta = ior;
    }
    else
    {
        eta = 1.0f / ior;
    }

    const float sin_theta_t_sqr =
        f_min(1.0f, (eta * eta) * (1.0f - cos_theta_i * cos_theta_i));
    const float cos_theta_t = std::sqrt(1.0f - sin_theta_t_sqr);

    float       t1 = eta * cos_theta_t;
    const float r_s = (cos_theta_i - t1) / (cos_theta_i + t1);

    float       t2 = eta * cos_theta_i;
    const float r_p = (t2 - cos_theta_t) / (t2 + cos_theta_t);

    return (r_s * r_s + r_p * r_p) * 0.5f;
}

vec3 reflectWorld(const vec3& wi, const vec3& nl)
{
#if 1
    //     assert(dot(wi, nl) > 0);
    return nl * (2.0f * dot(nl, wi)) - wi;
#else
    float cos_wi = dot(wi, nl);
    //     assert(cos_wi > 0);
    return nl * (2.0f * cos_wi) - wi;
#endif
}

vec3 reflectLocal(const vec3& wi)
{
    return vec3(-wi.x, -wi.y, wi.z);
}

OrthonormalFrame::OrthonormalFrame(const vec3& nl)
{
    setFromNormal(nl);
}

OrthonormalFrame::OrthonormalFrame()
{
    m_x = vec3(1, 0, 0);
    m_y = vec3(0, 1, 0);
    m_z = vec3(0, 0, 1);
}

void OrthonormalFrame::setFromNormal(const vec3& nl)
{
    assert(std::abs(length(nl) - 1.0f) < 1e-6f);

#if 0
    vec3 u, w, v = nl;
    if (nl.z < -0.9999999)
    {
        u = vec3(0.0, -1.0, 0.0);
        w = vec3(-1.0, 0.0, 0.0);
    }
    else
    {
        const float a = 1.0f / (1.0f + nl.z);
        const float b = -nl.x * nl.y * a;
        u = vec3(1.0f - nl.x * nl.x * a, b, -nl.x);
        w = vec3(b, 1.0f - nl.y * nl.y * a, -nl.y);
    }
    m_x = u;
    m_y = w;
    m_z = v;
#else
    // Building an Orthonormal Basis, Revisited
    // http://jcgt.org/published/0006/01/01/
    float       sign = copysignf(1.0f, nl.z);
    const float a = -1.0f / (sign + nl.z);
    const float b = nl.x * nl.y * a;
    m_x = vec3(1.0f + sign * nl.x * nl.x * a, sign * b, -sign * nl.x);
    m_y = vec3(b, sign + nl.y * nl.y * a, -nl.y);
    m_z = nl;
#endif
}

vec3 OrthonormalFrame::toWorld(const vec3& coords) const
{
    return m_x * coords.x + m_y * coords.y + m_z * coords.z;
}

vec3 OrthonormalFrame::toLocal(const vec3& coords) const
{
    return vec3(dot(coords, m_x), dot(coords, m_y), dot(coords, m_z));
}

const vec3& OrthonormalFrame::getNormal() const
{
    return m_z;
}

const vec3& OrthonormalFrame::getBitangent() const
{
    return m_y;
}

const vec3& OrthonormalFrame::getTangent() const
{
    return m_x;
}

float geometryTerm(const Vert& e0, const Vert& e1)
{
    vec3        dv = e1.p - e0.p;
    const float d2 = dot(dv, dv);
    dv = normalize(dv);
    float g = 1.0f;
    //     if (e0.obj && !
    //     todo_getMaterial(e0.obj->getMaterialID()).getBSDF()->isDelta())
    {
        //         g *= f_max(dot(e0.n, dv), 0); // this dims the bdpt result,
        //         todo : inspect the cause
        g *= dot(e0.n, dv);
    }
    //     if (e1.obj && !
    //     todo_getMaterial(e1.obj->getMaterialID()).getBSDF()->isDelta())
    {
        //         g *= f_max(-dot(e1.n, dv), 0); // this dims the bdpt result,
        //         todo : inspect the cause
        g *= dot(e1.n, dv);
    }
    return fabs(g) / (d2);

    // this is not ok since cosine factor is included depending on the materials
    // this has caused fireflies when specular material is involved
    //     return fabs(g) / (d2 * d2); // one d2 for the distance square,
    //     another for normalizing the two dv
}

float geometryTerm(const vec3& p0,
                   const vec3& n0,
                   const vec3& p1,
                   const vec3& n1)
{
    vec3        dv = p1 - p0;
    const float d2 = dot(dv, dv);
    dv = normalize(dv);
    float g = 1.0f;
    {
        g *= dot(n0, dv);
    }
    {
        g *= dot(n1, -dv);
    }
    return std::abs(g) / (d2);
}

float pdfWtoA(const vec3& v0, const vec3& v1, const vec3& n1)
{
    const vec3  dv = v0 - v1;
    const float d2 = dot(dv, dv);
    return f_max(dot(n1, dv), 0.0) /
           (d2 * std::sqrt(d2));  // sqrt(d2) for normalizing the vector
}

float pdfAtoW(const vec3& v0, const vec3& v1, const vec3& n1)
{
    const vec3  dv = v0 - v1;
    const float d2 = dot(dv, dv);
    return (d2 * std::sqrt(d2)) /
           f_max(dot(n1, dv), 0.0);  // sqrt(d2) for normalizing the vector
}

float directionToArea(const Vert& current, const Vert& next)
{
#if 1
    const vec3  dv = next.p - current.p;
    const float d2 = dot(dv, dv);
    return fabs(dot(next.n, dv)) /
           (d2 * std::sqrt(d2));  // sqrt(d2) for normalizing the vector
#else
    return pdfWtoA(1.0,
                   length(current.p - next.p),
                   fabs(dot(next.n, normalize(current.p - next.p))));
#endif
}

float pdfWtoA(float pdfW, float dist, float cosineThere)
{
    return pdfW * std::abs(cosineThere) / sq(dist);
}

float pdfAtoW(float pdfA, float dist, float cosineThere)
{
    return pdfA * sq(dist) / std::abs(cosineThere);
}

float __pdfWtoA(float pdfW, float dist_sq, float cosineThere)
{
    return pdfW * fabs(cosineThere) / dist_sq;
}
float __pdfAtoW(float pdfA, float dist_sq, float cosineThere)
{
    return pdfA * dist_sq / fabs(cosineThere);
}

bool isAligned(const vec3& a, const vec3& b)
{
    return dot(a, b) > 0.95;
}
