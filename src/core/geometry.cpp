#include "geometry.h"
#include "lightpath.h"
#include "material.h"
#include "mesh.h"
#include "numeric.h"
#include "scene.h"
#include "triangle.h"

float cosineFactor(const vec3 wi, const vec3& nl) {
    return f_max(dot(wi, nl), 0.0f);
}

float fresnelReflection(float cos_theta_i, float ior)
{
    if (ior < 0)
    {
        return 1;
    }

    float eta; // eta_i / eta_t
    if (cos_theta_i < 0)
    {
        cos_theta_i = -cos_theta_i;
        eta = ior;
    }
    else
    {
        eta = 1.0 / ior;
    }

    const float sin_theta_t_sqr = f_min(1.0f, (eta * eta) * (1.0 - cos_theta_i * cos_theta_i));
    const float cos_theta_t = std::sqrtf(1.0 - sin_theta_t_sqr);

    float t1 = eta * cos_theta_t;
    const float r_s = (cos_theta_i - t1) / (cos_theta_i + t1);

    float t2 = eta * cos_theta_i;
    const float r_p = (t2 - cos_theta_t) / (t2 + cos_theta_t);

    return (r_s * r_s + r_p * r_p) * 0.5;
}

vec3 reflectWorld(const vec3& wi, const vec3& nl)
{
#if 1
//     assert(dot(wi, nl) > 0);
    return nl * (2.0 * dot(nl, wi)) - wi;
#else
    float cos_wi = dot(wi, nl);
//     assert(cos_wi > 0);
    return nl * (2.0 * cos_wi) - wi;
#endif
}

vec3 reflectLocal(const vec3& wi)
{
    return vec3(-wi.x, -wi.y, wi.z);
}

OrthonormalFrame::OrthonormalFrame(const vec3& _z)
{
    setFromNormal(_z);
}

OrthonormalFrame::OrthonormalFrame()
{
    m_x = vec3(1, 0, 0);
    m_y = vec3(0, 1, 0);
    m_z = vec3(0, 0, 1);
}

// void generateOrthoBasis(vec3 &u, vec3 &v, vec3 w)
// {
//     //w is normal
//     vec3 coVec;
//     if (fabs(w.x) <= fabs(w.y))
//         if (fabs(w.x) <= fabs(w.z)) coVec = vec3(0, -w.z, w.y);
//         else coVec = vec3(-w.y, w.x, 0);
//     else if (fabs(w.y) <= fabs(w.z)) coVec = vec3(-w.z, 0, w.x);
//     else coVec = vec3(-w.y, w.x, 0);
//     coVec.norm();
//     u = cross(w, coVec);  //bug caused by inconsistent convention of operator overloading and silent conversion
//     v = cross(w, u);      //bug caused by inconsistent convention of operator overloading and silent conversion
// }

void OrthonormalFrame::setFromNormal(const vec3& _z)
{
    assert(std::abs(length(_z) - 1) < 1e-6);

#if 0 // causing visual bug, todo : investigate the cause
    m_z = _z;

    //w is normal
    vec3 coVec;
    if (fabs(_z.x) <= fabs(_z.y))
    {
        if (fabs(_z.x) <= fabs(_z.z))
        {
            coVec = vec3(0, -_z.z, _z.y);
        }
        else
        {
            coVec = vec3(-_z.y, _z.x, 0);
        }
    }
    else
    {
        if (fabs(_z.y) <= fabs(_z.z))
        {
            coVec = vec3(-_z.z, 0, _z.x);
        }
        else
        {
            coVec = vec3(-_z.y, _z.x, 0);
        }
    }
    coVec.norm();
    m_x = normalize(cross(_z, coVec));  //bug caused by inconsistent convention of operator overloading and silent conversion
    m_y = normalize(cross(_z, m_x));      //bug caused by inconsistent convention of operator overloading and silent conversion
#else
    vec3 u, w, v = _z;
    if (_z.z < -0.9999999)
    {
        u = vec3(0.0, -1.0, 0.0);
        w = vec3(-1.0, 0.0, 0.0);
    }
    else
    {
        const float a = 1.0 / (1.0 + _z.z);
        const float b = -_z.x * _z.y * a;
        u = vec3(1.0 - _z.x * _z.x * a, b, -_z.x);
        w = vec3(b, 1.0 - _z.y * _z.y * a, -_z.y);
    }
    m_x = u;
    m_y = w;
    m_z = v;
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
    vec3 dv = e1.p - e0.p;
    const float d2 = dot(dv, dv);
    dv = normalize(dv);
    float g = 1.0f;
//     if (e0.obj && ! todo_getMaterial(e0.obj->getMaterialID()).getBSDF()->isDelta())
    {
//         g *= f_max(dot(e0.n, dv), 0); // this dims the bdpt result, todo : inspect the cause
        g *= dot(e0.n, dv);
    }
//     if (e1.obj && ! todo_getMaterial(e1.obj->getMaterialID()).getBSDF()->isDelta())
    {
//         g *= f_max(-dot(e1.n, dv), 0); // this dims the bdpt result, todo : inspect the cause
        g *= dot(e1.n, dv);
    }
    return fabs(g) / (d2);

    // this is not ok since cosine factor is included depending on the materials
    // this has caused fireflies when specular material is involved
//     return fabs(g) / (d2 * d2); // one d2 for the distance square, another for normalizing the two dv
}

float geometryTerm(const vec3& p0, const vec3& n0, const vec3& p1, const vec3& n1)
{
    vec3 dv = p1 - p0;
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
    const vec3 dv = v0 - v1;
    const float d2 = dot(dv, dv);
    return f_max(dot(n1, dv), 0.0) / (d2 * std::sqrt(d2)); // sqrt(d2) for normalizing the vector
}

float pdfAtoW(const vec3& v0, const vec3& v1, const vec3& n1)
{
    const vec3 dv = v0 - v1;
    const float d2 = dot(dv, dv);
    return (d2 * std::sqrt(d2)) / f_max(dot(n1, dv), 0.0); // sqrt(d2) for normalizing the vector
}

float directionToArea(const Vert& current, const Vert& next)
{
#if 1
    const vec3 dv = next.p - current.p;
    const float d2 = dot(dv, dv);
    return fabs(dot(next.n, dv)) / (d2 * std::sqrt(d2)); // sqrt(d2) for normalizing the vector
#else
    return pdfWtoA(1.0, length(current.p - next.p), fabs(dot(next.n, normalize(current.p - next.p))));
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

void loadModel(Scene& scene, std::string filename, const glm::mat4& transform, const Material& mat)
{
    CMesh obj(filename.c_str());
    {
        for (int n = 0; n < obj.m_verts.size(); n++)
        {
            glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
            v = transform * v;
            obj.m_verts[n] = vec3(v.x, v.y, v.z);
        }
        obj.genNormals();
    }

    int id = todo_addMaterial(mat);
    for (int i = 0; i < obj.m_faces.size(); i++)
    {
        int a = obj.m_faces[i].x;
        int b = obj.m_faces[i].y;
        int c = obj.m_faces[i].z;

        TriangleObject tri(
            vec3(obj.m_verts[a].x, (obj.m_verts[a].y), obj.m_verts[a].z),
            vec3(obj.m_verts[b].x, (obj.m_verts[b].y), obj.m_verts[b].z),
            vec3(obj.m_verts[c].x, (obj.m_verts[c].y), obj.m_verts[c].z),
            vec3(obj.m_vertex_normals[a].x, (obj.m_vertex_normals[a].y), obj.m_vertex_normals[a].z),
            vec3(obj.m_vertex_normals[b].x, (obj.m_vertex_normals[b].y), obj.m_vertex_normals[b].z),
            vec3(obj.m_vertex_normals[c].x, (obj.m_vertex_normals[c].y), obj.m_vertex_normals[c].z)
            );
        tri.setMaterialID(id);

        scene.add(tri);
    }
}
