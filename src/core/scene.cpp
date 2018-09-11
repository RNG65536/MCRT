#include <iostream>
#include <numeric>
#include "aabb.h"
#include "camera.h"
#include "envmap.h"
#include "geometry.h"
#include "intersection.h"
#include "material.h"
#include "numeric.h"
#include "ray.h"
#include "sample.h"
#include "scene.h"
#include "scenebvh.h"

EnvmapLoader* Scene::getEnvmap() const {
    return m_envmap;
}

Scene::~Scene()
{
    clear();
}

Scene::Scene()
{
//     m_scene_bvh = new SceneBVH(m_triangles);
//     m_envmap = new EnvmapLoader;
}

void Scene::sampleLightAreaAndDirection(vec3& position, vec3& normal, float& pdfA, vec3& direction, float& pdfW, const TriangleObject *& obj) const
{
    // sample the triangles proportional to area
    int i = m_light_sampler->sample();
    const TriangleObject& tri = *m_lights[i];
    obj = &tri;

    // sample inside the triangle uniformly
    float u, v;
    vec3 p = tri.samplePoint(u, v);
    vec3 n = tri.normal(u, v);
    vec3 fn = tri.getFaceNormal();
    position = p;
    normal = n;
#if 0
    pdfA = 1.0 / getLightArea();
#else
    pdfA = (1.0 / tri.getArea()) * m_light_sampler->pdf(i);
#endif

#if 1
    // diffuse emission
    const vec3 c = sampleHemisphereCosine(randf(), randf());
    OrthonormalFrame onb(n);
    direction = onb.toWorld(c);

    // avoid inward emission
    if (dot(direction, fn) < 0)
    {
        OrthonormalFrame onb(fn);
        direction = onb.toWorld(c);
        normal = fn;
    }

    pdfW = sampleHemisphereCosinePdfW(normal, direction);
#else
    // not diffuse emission, but not much difference
    return Ray(p, fn);
#endif
}

HitInfo Scene::intersect(const Ray& ray) const
{
    return m_scene_bvh->intersect(ray);
}

bool Scene::occluded(const Ray& ray, float tfar) const
{
    return m_scene_bvh->occluded(ray, tfar);
}

void Scene::clear()
{
    delete m_scene_bvh;
    delete m_light_sampler;
    delete m_envmap;
    delete m_camera;
    m_triangles.clear();
    m_lights.clear();
    m_light_area = 0;
    m_bounding_sphere = SceneSphere();
}

const TriangleObject& Scene::getTriangle(size_t n) const
{
    return m_triangles[n];
}

TriangleObject& Scene::getTriangle(size_t n)
{
    return m_triangles[n];
}

void Scene::add(const TriangleObject& object)
{
    m_triangles.push_back(object);
}

float Scene::getTriangleArea(size_t n) const
{
    return m_triangles[n].getArea();
}

size_t Scene::numTrangles() const
{
    return m_triangles.size();
}

vec3 Scene::getCenter() const
{
    vec3 bmin(1e10, 1e10, 1e10);
    vec3 bmax(-1e10, -1e10, -1e10);

    size_t numtris = m_triangles.size();
    for (int n = 0; n < numtris; n++){
        bmin = f_min(bmin, m_triangles[n].a);
        bmin = f_min(bmin, m_triangles[n].b);
        bmin = f_min(bmin, m_triangles[n].c);
        bmax = f_max(bmax, m_triangles[n].a);
        bmax = f_max(bmax, m_triangles[n].b);
        bmax = f_max(bmax, m_triangles[n].c);
    }

    printf("scene bounding box:\n");
    printf("%f, %f, %f\n", bmin.x, bmin.y, bmin.z);
    printf("%f, %f, %f\n", bmax.x, bmax.y, bmax.z);

    return (bmin + bmax) * 0.5;
}

void Scene::finalize()
{
    //build bounding sphere
    AABB scene_bbox;
    for (auto& t : m_triangles)
    {
        AABB tri_bbox = getAABB(t);
        scene_bbox.enclose(tri_bbox);
    }
    m_bounding_sphere.m_cenetr = scene_bbox.getCenter();
    m_bounding_sphere.m_radius = scene_bbox.diagonalLength() * 0.5f;
    m_bounding_sphere.m_inv_radius_sqr = 1.0f / sq(m_bounding_sphere.m_radius);
    std::cout << "scene sphere : ( " << m_bounding_sphere.m_cenetr.x << ", "
                                << m_bounding_sphere.m_cenetr.y << ", "
                                << m_bounding_sphere.m_cenetr.z << " ), "
                                << m_bounding_sphere.m_radius << std::endl;

    m_scene_bvh = new SceneBVH(m_triangles);

    int n = m_triangles.size();

    m_lights.clear();
    int tri_id = 0, light_id = 0;
    for (auto& t : m_triangles)
    {
        if (todo_getMaterial(t.materialID()).type() == LGHT)
        {
            m_lights.push_back(&t);
            _primitiveID_to_lightID.insert(std::make_pair(tri_id, light_id));
            ++light_id;
        }
        ++tri_id;
    }

    std::vector<float> areas(m_lights.size());
    for (int i = 0; i < areas.size(); i++)
    {
        areas[i] = m_lights[i]->getArea();
    }
//     TKahanAdder sum_areas(0.0);
//     for (auto t : areas)
//     {
//         sum_areas.add(t);
//     }
//     m_light_area = sum_areas.sum;
    m_light_area = std::accumulate(areas.begin(), areas.end(), 0.0f);
    std::cout << "light area = " << m_light_area << std::endl;

    m_light_sampler = new DiscreteSampler(areas);
}

float Scene::getLightArea() const
{
    return m_light_area;
}

CCamera * Scene::getCamera() const
{
    return m_camera;
}

void Scene::setEnvmap(const std::string filename)
{
    if (m_envmap)
    {
        delete m_envmap;
    }
    m_envmap = new EnvmapLoader(filename);
}

void Scene::setEnvmap(const std::string filename, int width, int height)
{
    if (m_envmap)
    {
        delete m_envmap;
    }
    m_envmap = new EnvmapLoader(filename, width, height);
}

void Scene::setCamera(const vec3& origin, const vec3& target, int width, int height, float fov_y)
{
    if (m_camera)
    {
        delete m_camera;
    }
    m_camera = new CCamera(origin, target, vec3(0, 1, 0), width, height, fov_y);
}

const SceneSphere& Scene::sceneSphere() const
{
    return m_bounding_sphere;
}

float Scene::lightPdfW(const vec3& pos, const vec3& nl, const vec3& dir) const
{
    Ray shadow_ray(pos, dir);
    HitInfo hit = intersect(shadow_ray);
    if (!hit) // not hitting anything, pdfW_light = 0
    {
        return 0;
    }

    const TriangleObject *tri = hit.getObject();
    if (todo_getMaterial(tri->materialID()).type() != LGHT) //  not hitting a light surface, pdfW_light = 0
    {
        return 0;
    }

    float cos_light = dot(hit.getShadingNormal(), -dir);
    if (cos_light < 1e-6f) // not visible, pdfW_light = 0
    {
        return 0;
    }

    int tri_id = hit.primitiveID();
    int light_id = _primitiveID_to_lightID.at(tri_id);
    float pdfA = (1.0 / tri->getArea()) * m_light_sampler->pdf(light_id);
    float pdfW = pdfA * sq(hit.getDistance()) / cos_light; // TODO : use pdfAtoW

    return pdfW;
}

void Scene::sampleLight(float& pdfLight, int& light_id, vec3& sample_pos, vec3& sample_normal,
    float rnd0, float rnd1, float rnd2) const
{
    // sample the triangles proportional to area
    int i = m_light_sampler->sample(rnd0);
    const TriangleObject& tri = *m_lights[i];
    light_id = i;

    // sample inside the triangle uniformly
    vec2 uv = sampleTriangleUniform(rnd1, rnd2);
    vec3 p = tri.point(uv.x, uv.y);
    vec3 n = tri.normal(uv.x, uv.y);

//     float u, v;
//     vec3 p = tri.samplePoint(u, v);
//     vec3 n = tri.normal(u, v);


    vec3 fn = tri.getFaceNormal();
    sample_pos = p;
    sample_normal = n;

    pdfLight = (1.0 / tri.getArea()) * m_light_sampler->pdf(i);
}

vec3 Scene::sampleLightArea(const vec3& position, /*const vec3& normal,*/ float& pdfW, int& light_id, vec3 *sample_pos, vec3 *sample_normal) const
{
    // sample the triangles proportional to area
    int i = m_light_sampler->sample();
    const TriangleObject& tri = *m_lights[i];
    light_id = i;

    // sample inside the triangle uniformly
    float u, v;
    vec3 p = tri.samplePoint(u, v);
    vec3 n = tri.normal(u, v);
    vec3 fn = tri.getFaceNormal();

    if (sample_pos)
    {
        *sample_pos = p;
    }
    if (sample_normal)
    {
        *sample_normal = n;
    }

#if 0
    pdfA = 1.0 / getLightArea();
#else
    float pdfA = (1.0 / tri.getArea()) * m_light_sampler->pdf(i);
//     cout << pdfA << ", " << 1.0 / getLightArea() << endl;
#endif

    vec3 wo = p - position;
    float dist2 = wo.lengthSquared();
    wo = normalize(wo);
    float cos_light = dot(n, -wo);

    pdfW = __pdfAtoW(pdfA, dist2, cos_light);
    if (cos_light <= 0)
    {
        pdfW = 0;
//         light_id = -1;
    }
    return wo;
}

vec3 Scene::lightRadiance(const vec3& pos, /*const vec3& nl,*/ const vec3& dir, const int light_id /*= -1*/) const
{
    HitInfo hit = intersect(Ray(pos, dir));

    if (!hit
        || todo_getMaterial(hit.materialID()).type() != LGHT // not hit
        || dot(hit.getShadingNormal(), -dir) < 0 // or not visible
        || (light_id >= 0 && light_id != _primitiveID_to_lightID.at(hit.primitiveID())) // this is for active light sampling
        )
    {
        return vec3(0, 0, 0);
    }

    else
    {
        return todo_getMaterial(hit.getObject()->materialID()).color();
    }
}

