#pragma once

#include <glm/glm.hpp>
#include <unordered_map>
#include "camera.h"
#include "constants.h"
#include "envmap.h"
#include "intersection.h"
#include "material.h"
#include "ray.h"
#include "sample.h"
#include "triangle.h"
#include "vectors.h"

// scene management

class SceneBVH;
class DiscreteSampler;

class SceneSphere
{
public:
    vec3  m_cenetr;
    float m_radius;
    float m_inv_radius_sqr;
};

class Scene
{
public:
    Scene();
    ~Scene();

    void clear();
    void finalize();
    void add(const TriangleObject& object);
    void setCamera(const vec3& origin,
                   const vec3& target,
                   int         width,
                   int         height,
                   float       fov_y);

    // environment map
    void          setEnvmap(const std::string filename);
    void          setEnvmap(const std::string filename, int width, int height);
    EnvmapLoader* envmap() const;

    // other lights
    float lightArea() const;
    float lightPdfW(const vec3& pos, const vec3& nl, const vec3& dir) const;
    vec3  lightRadiance(const vec3& pos,
                        const vec3& dir,
                        const int   light_id = -1) const;
    void  sampleLightAreaAndDirection(
         vec3&                  position,
         vec3&                  normal,
         float&                 pdfA,
         vec3&                  direction,
         float&                 pdfW,
         const TriangleObject*& obj,
         StandardSampler&       sampler) const;  // assuming diffuse emission
    vec3 sampleLightArea(const vec3& position,
                         float&      pdfW,
                         int&        light_id,
                         vec3*       sample_pos = nullptr,
                         vec3*       sample_normal = nullptr) const;
    void sampleLight(float& pdfLight,
                     int&   light_id,
                     vec3&  sample_pos,
                     vec3&  sample_normal,
                     float  rnd0,
                     float  rnd1,
                     float  rnd2) const;

    // intersection
    HitInfo intersect(const Ray& ray) const;
    bool    occluded(const Ray& ray, float tfar = NUM_INFINITY) const;

    // scene objects
    TriangleObject&       triangle(size_t n);
    const TriangleObject& triangle(size_t n) const;
    float                 triangleArea(size_t n) const;
    Camera*               camera() const;
    size_t                numTrangles() const;
    vec3                  boxCenter() const;
    const SceneSphere&    boundingSphere() const;

    void loadModel(const std::string& filename,
                   const glm::mat4&   transform,
                   const Material&    mat);

private:
    SceneBVH*                    m_scene_bvh = nullptr;
    std::vector<TriangleObject>  m_triangles;
    std::vector<TriangleObject*> m_lights;
    float                        m_light_area;
    DiscreteSampler*             m_light_sampler = nullptr;
    EnvmapLoader*                m_envmap = nullptr;
    Camera*                      m_camera = nullptr;
    SceneSphere                  m_bounding_sphere;
    vec3                         m_box_center;
    std::unordered_map<int, int> m_primitiveID_to_lightID;
};
