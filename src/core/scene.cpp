#include <iostream>
#include <numeric>
#include "aabb.h"
#include "camera.h"
#include "envmap.h"
#include "geometry.h"
#include "intersection.h"
#include "material.h"
#include "mesh.h"
#include "numeric.h"
#include "ray.h"
#include "sample.h"
#include "scene.h"
#include "scenebvh.h"

EnvmapLoader* Scene::envmap() const
{
    return m_envmap;
}

Scene::~Scene()
{
    clear();
}

Scene::Scene()
{
}

void Scene::sampleLightAreaAndDirection(vec3&                  position,
                                        vec3&                  normal,
                                        float&                 pdfA,
                                        vec3&                  direction,
                                        float&                 pdfW,
                                        const TriangleObject*& obj,
                                        StandardSampler&       sampler) const
{
    // sample the triangles proportional to area
    int                   i = m_light_sampler->sample(sampler.next());
    const TriangleObject& tri = *m_lights[i];
    obj = &tri;

    // sample inside the triangle uniformly
    float u, v;
    vec3  p = tri.samplePoint(u, v, sampler.next(), sampler.next());
    vec3  n = tri.normal(u, v);
    vec3  fn = tri.normal();
    position = p;
    normal = n;
#if 0
    pdfA = 1.0 / lightArea();
#else
    pdfA = (1.0f / tri.area()) * m_light_sampler->pdf(i);
#endif

#if 1
    // diffuse emission
    const vec3       c = sampleHemisphereCosine(sampler.next(), sampler.next());
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

const TriangleObject& Scene::triangle(size_t n) const
{
    return m_triangles[n];
}

TriangleObject& Scene::triangle(size_t n)
{
    return m_triangles[n];
}

void Scene::add(const TriangleObject& object)
{
    m_triangles.push_back(object);
}

float Scene::triangleArea(size_t n) const
{
    return m_triangles[n].area();
}

size_t Scene::numTrangles() const
{
    return m_triangles.size();
}

vec3 Scene::boxCenter() const
{
    return m_box_center;
}

void Scene::finalize()
{
    // build bounding sphere
    AABB scene_bbox;
    for (auto& t : m_triangles)
    {
        AABB tri_bbox = getAABB(t);
        scene_bbox.expandBy(tri_bbox);
    }
    m_bounding_sphere.m_cenetr = scene_bbox.center();
    m_bounding_sphere.m_radius = scene_bbox.diagonal().length() * 0.5f;
    m_bounding_sphere.m_inv_radius_sqr = 1.0f / sq(m_bounding_sphere.m_radius);
    std::cout << "scene sphere : ( " << m_bounding_sphere.m_cenetr.x << ", "
              << m_bounding_sphere.m_cenetr.y << ", "
              << m_bounding_sphere.m_cenetr.z << " ), "
              << m_bounding_sphere.m_radius << std::endl;

    assert(nullptr == m_scene_bvh);
    m_scene_bvh = new SceneBVH(m_triangles);

    int n = static_cast<int>(m_triangles.size());

    m_lights.clear();
    int tri_id = 0, light_id = 0;
    for (auto& t : m_triangles)
    {
        if (todo_getMaterial(t.materialID()).type() == LGHT)
        {
            m_lights.push_back(&t);
            m_primitiveID_to_lightID.insert(std::make_pair(tri_id, light_id));
            ++light_id;
        }
        ++tri_id;
    }

    std::vector<float> areas(m_lights.size());
    for (int i = 0; i < areas.size(); i++)
    {
        areas[i] = m_lights[i]->area();
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

    AABB   bbox;
    size_t num_triangles = m_triangles.size();
    for (int n = 0; n < num_triangles; n++)
    {
        bbox.expandBy(m_triangles[n].a);
        bbox.expandBy(m_triangles[n].b);
        bbox.expandBy(m_triangles[n].c);
    }

    printf("scene bounding box:\n");
    printf("%f, %f, %f\n", bbox.min().x, bbox.min().y, bbox.min().z);
    printf("%f, %f, %f\n", bbox.max().x, bbox.max().y, bbox.max().z);

    m_box_center = bbox.center();
}

float Scene::lightArea() const
{
    return m_light_area;
}

Camera* Scene::camera() const
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

void Scene::setCamera(
    const vec3& origin, const vec3& target, int width, int height, float fov_y)
{
    if (m_camera)
    {
        delete m_camera;
    }
    m_camera = new Camera(origin, target, vec3(0, 1, 0), width, height, fov_y);
}

const SceneSphere& Scene::boundingSphere() const
{
    return m_bounding_sphere;
}

float Scene::lightPdfW(const vec3& pos, const vec3& nl, const vec3& dir) const
{
    Ray     shadow_ray(pos, dir);
    HitInfo hit = intersect(shadow_ray);
    if (!hit)  // not hitting anything, pdfW_light = 0
    {
        return 0;
    }

    const TriangleObject* tri = hit.triangleObject();
    if (todo_getMaterial(tri->materialID()).type() !=
        LGHT)  //  not hitting a light surface, pdfW_light = 0
    {
        return 0;
    }

    float cos_light = dot(hit.shadingNormal(), -dir);
    if (cos_light < 1e-6f)  // not visible, pdfW_light = 0
    {
        return 0;
    }

    int   tri_id = hit.primitiveID();
    int   light_id = m_primitiveID_to_lightID.at(tri_id);
    float pdfA = (1.0f / tri->area()) * m_light_sampler->pdf(light_id);
    float pdfW = pdfA * sq(hit.distance()) / cos_light;  // TODO : use pdfAtoW

    return pdfW;
}

void Scene::sampleLight(float& pdfLight,
                        int&   light_id,
                        vec3&  sample_pos,
                        vec3&  sample_normal,
                        float  rnd0,
                        float  rnd1,
                        float  rnd2) const
{
    // sample the triangles proportional to area
    int                   i = m_light_sampler->sample(rnd0);
    const TriangleObject& tri = *m_lights[i];
    light_id = i;

    // sample inside the triangle uniformly
    vec2 uv = sampleTriangleUniform(rnd1, rnd2);
    vec3 p = tri.point(uv.x, uv.y);
    vec3 n = tri.normal(uv.x, uv.y);

    //     float u, v;
    //     vec3 p = tri.samplePoint(u, v);
    //     vec3 n = tri.normal(u, v);

    vec3 fn = tri.normal();
    sample_pos = p;
    sample_normal = n;

    pdfLight = (1.0f / tri.area()) * m_light_sampler->pdf(i);
}

vec3 Scene::sampleLightArea(const vec3&                    position,
                            /*const vec3& normal,*/ float& pdfW,
                            int&                           light_id,
                            vec3*                          sample_pos,
                            vec3*                          sample_normal) const
{
    // sample the triangles proportional to area
    int                   i = m_light_sampler->sample();
    const TriangleObject& tri = *m_lights[i];
    light_id = i;

    // sample inside the triangle uniformly
    float u, v;
    vec3  p = tri.samplePoint(u, v);
    vec3  n = tri.normal(u, v);
    vec3  fn = tri.normal();

    if (sample_pos)
    {
        *sample_pos = p;
    }
    if (sample_normal)
    {
        *sample_normal = n;
    }

#if 0
    pdfA = 1.0 / lightArea();
#else
    float pdfA = (1.0f / tri.area()) * m_light_sampler->pdf(i);
//     cout << pdfA << ", " << 1.0 / lightArea() << endl;
#endif

    vec3  wo = p - position;
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

vec3 Scene::lightRadiance(const vec3&                     pos,
                          /*const vec3& nl,*/ const vec3& dir,
                          const int light_id /*= -1*/) const
{
    HitInfo hit = intersect(Ray(pos, dir));

    if (!hit || todo_getMaterial(hit.materialID()).type() != LGHT  // not hit
        || dot(hit.shadingNormal(), -dir) < 0  // or not visible
        || (light_id >= 0 &&
            light_id !=
                m_primitiveID_to_lightID.at(
                    hit.primitiveID()))  // this is for active light sampling
    )
    {
        return vec3(0, 0, 0);
    }

    else
    {
        return todo_getMaterial(hit.triangleObject()->materialID()).color();
    }
}

void Scene::loadModel(const std::string& filename,
                      const glm::mat4&   transform,
                      const Material&    mat)
{
    Mesh obj(filename.c_str());
    {
        for (int n = 0; n < obj.m_verts.size(); n++)
        {
            glm::vec4 v = glm::vec4(
                obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
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
            vec3(obj.m_vertex_normals[a].x,
                 (obj.m_vertex_normals[a].y),
                 obj.m_vertex_normals[a].z),
            vec3(obj.m_vertex_normals[b].x,
                 (obj.m_vertex_normals[b].y),
                 obj.m_vertex_normals[b].z),
            vec3(obj.m_vertex_normals[c].x,
                 (obj.m_vertex_normals[c].y),
                 obj.m_vertex_normals[c].z));
        tri.setMaterialID(id);

        this->add(tri);
    }
}
