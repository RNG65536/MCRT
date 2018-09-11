#include <cmath>
#include "aabb.h"
#include "constants.h"
#include "intersection.h"
#include "numeric.h"
#include "ray.h"
#include "triangle.h"

vec3 TriangleObject::normal(float u, float v) const {
    return normalize((1 - u - v) * na + u * nb + v * nc);
}

vec3 TriangleObject::point(float u, float v) const
{
    return (1 - u - v) * a + u * b + v * c;
}

// vec3 TriangleObject::normal(const vec3& p0) const
// {
//     if (na.length() < NUM_EPS && nb.length() < NUM_EPS && nc.length() < NUM_EPS)
//         return n;
//     //         Real area = ((a-c)^(b-c)).length();
//     float area_a = length(cross((b - p0), (c - p0)));
//     float area_b = length(cross((c - p0), (a - p0)));
//     float area_c = length(cross((a - p0), (b - p0)));
//     return (area_a * na + area_b * nb + area_c * nc).norm();
// }

vec3 TriangleObject::samplePoint(float& u, float& v) const
{
    u = randf();
    v = randf();
    if (u + v > 1)
    {
        u = 1 - u;
        v = 1 - v;
    }
    return a + (b - a) * u + (c - a) * v;
}

vec3 TriangleObject::samplePoint() const
{
    float u = randf();
    float v = randf();
    if (u + v > 1)
    {
        u = 1 - u;
        v = 1 - v;
    }
    return a + (b - a) * u + (c - a) * v;
}

vec3 TriangleObject::getFaceNormal() const
{
    return n;
}

float TriangleObject::getArea() const
{
    return 0.5 * length(cross((b - a), (c - a)));
}

vec3 TriangleObject::getCenter() const
{
    return (a + b + c) * 0.3333333333333333f;
}

HitInfo TriangleObject::intersect(const HitInfo &hit, const Ray &ray) const
{
#if 1 //Moller-Trumbore (not much improvement though)
    const vec3& v0 = a;
    const vec3& v1 = b;
    const vec3& v2 = c;
    vec3 edge1 = v1 - v0;
    vec3 edge2 = v2 - v0;
    vec3 pvec = cross(ray.d, edge2);
    float det = dot(edge1, pvec);
    if (fabs(det) < NUM_EPS_DET) return hit;//do not try to cull back for transparent objects
    float invDet = 1.0f / det;
    vec3 tvec = ray.o - v0;
    float u = dot(tvec, pvec) * invDet;
    if (u < 0 || u > 1) return hit;
    vec3 qvec = cross(tvec, edge1);
    float v = dot(ray.d, qvec) * invDet;
    if (v < 0 || u + v > 1) return hit;
    float t = dot(edge2, qvec) * invDet;
    //         if(t < eps || t > hit.t) return hit;
    if (t < NUM_EPS_RAY || t > hit.getDistance()) return hit;
    //         return Hit_CUDA(t,normal(ray(t)),this);
    //         return Hit_CUDA(t,normal(u*edge1+v*edge2),this);
    //         return Hit_CUDA(t,normal(ray(t)),this);
//     return HitInfo(t, normal(u, v), this);
    return HitInfo(t, u, v, normal(u, v), this);
#else
    //first find intersection with plane
    float d0 = dot(n, ray.d);
    bool front = d0 < 0;
#ifdef BACK_CULLING
    if (!front) return hit;
#endif
    float t = dot(n, a - ray.o) / d0;
    // 		if (t < eps*1e3 || t>hit.t) return hit;//for BVH, a miss-hit/inverse-hit doesn't take place
    // 		if (t < eps || t>hit.t) return hit;//for BVH, a miss-hit/inverse-hit doesn't take place
    if (t < NUM_EPS || t>hit.t) return hit;//for BVH, a miss-hit/inverse-hit doesn't take place
    //second, test if intersection is within triangle
#if FIX_FIREFLY
    if (getArea() < NUM_EPS) return hit;//::fix bug elegantly (the fire flies)
#endif
    //         if(getArea()<eps*0.001) return hit;//::fix bug elegantly (the fire flies) //using a much smaller threshold makes no significant difference, 
    //so the dark edges must be related to the following routine that interpolates the normal
    vec3 po(ray(t) - ray.o), ao(a - ray.o), bo(b - ray.o), co(c - ray.o);
#if 1
    vec3 n1(bo^ao);        if (dot(po, n1)*(front ? 1 : -1) < 0) return hit;
    vec3 n2(co^bo);        if (dot(po, n2)*(front ? 1 : -1) < 0) return hit;
    vec3 n3(ao^co);        if (dot(po, n3)*(front ? 1 : -1) < 0) return hit;
#else//wireframe
    vec3 n1(bo^ao);        if (dot(po, n1)*(front ? 1 : -1) < 0) return hit;
    vec3 n2(co^bo);        if (dot(po, n2)*(front ? 1 : -1) < 0) return hit;
    vec3 n3(ao^co);        if (dot(po, n3)*(front ? 1 : -1) < 0) return hit;
    vec3 p0(ray(t));
    Real area_a = ((b - p0) ^ (c - p0)).length();
    Real area_b = ((c - p0) ^ (a - p0)).length();
    Real area_c = ((a - p0) ^ (b - p0)).length();
    Real area = ((a - c) ^ (b - c)).length();
    if (!(area_a / area < .08 ||
        area_b / area < .08 ||
        area_c / area < .08)) return hit;
#endif
    //         return Hit(t,n,this);
    return HitInfo(t, normal(ray(t)), this);
#endif
}

TriangleObject::TriangleObject(vec3 a_, vec3 b_, vec3 c_, vec3 na_, vec3 nb_, vec3 nc_)
    : a(a_), b(b_), c(c_), n(normalize(cross((a - c), (b - c)))), na(na_), nb(nb_), nc(nc_)
{

}

TriangleObject::TriangleObject(vec3 a_, vec3 b_, vec3 c_)
    : a(a_), b(b_), c(c_), n(normalize(cross((a - c), (b - c)))), na(n), nb(n), nc(n)
{

}

TriangleObject::TriangleObject()
{

}

// void TriangleObject::setMat(const Material& m_)
// {
//     mat = m_;
// }
// 
// void TriangleObject::setMat(vec3 cl_ /*= vec3(1, 1, 1)*/, float emission_ /*= 0*/, MaterialType type_ /*= REFR*/)
// {
//     mat = Material(cl_, emission_, type_);
// }

// const Material& TriangleObject::getMaterial() const
// {
//     return mat;
// }

// void TriangleObject::setMaterial(const Material& _mat)
// {
//     mat = _mat;
// }

void TriangleObject::setMaterialID(int id)
{
    mat_id = id;
}

int TriangleObject::materialID() const
{
    return mat_id;
}

vec3 TriangleObject::consistentNormal(const vec3& wi, float u, float v,
    vec3* refl, float *cos_wi) const
{
    const vec3 np = normal(u, v);
//     if (dot(wi, np) > 0)
//     {
//         return np;
//     }

    // for alpha = 0, consistent normal becomes shading normal
    float alpha = ((1.0 - u - v) * m_vertex_alpha.x
        + u * m_vertex_alpha.y + v * m_vertex_alpha.z);

    float q = sq(1.0 - (2.0 / NUM_PI) * alpha)
        / (1.0 + 2.0 * (1.0 - 2.0 / NUM_PI) * alpha);
    float b = dot(wi, np);
    float g = 1.0 + q * (b - 1.0);
    float rho = sqrt(q * (1.0 + g) / (1.0 + b));
    const vec3 r = (g + rho * b) * np - rho * wi;

    // reflected by nc
    if (refl)
    {
        *refl = r;
    }

    // dot(wi, nc)
    if (cos_wi)
    {
        *cos_wi = sqrt((1.0 + (g + rho * b) * b - rho) * 0.5);
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

