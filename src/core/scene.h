#ifndef scene_h__
#define scene_h__

#include <unordered_map>
#include "constants.h"
#include "vectors.h"
#include "triangle.h"

class SceneBVH;
class EnvmapLoader;
class DiscreteSampler;
class CCamera;

class SceneSphere
{
public:
    vec3 m_cenetr;
    float m_radius;
    float m_inv_radius_sqr;
};

class Scene
{
private:
    SceneBVH *m_scene_bvh = nullptr;
    std::vector<TriangleObject> m_triangles;
    std::vector<TriangleObject*> m_lights;
    float m_light_area;
    DiscreteSampler *m_light_sampler = nullptr;
    EnvmapLoader *m_envmap = nullptr;
    CCamera *m_camera = nullptr;
    SceneSphere m_bounding_sphere;
    std::unordered_map<int, int> _primitiveID_to_lightID;

public:
    Scene();
    ~Scene();

    void clear();
    void finalize();
    void add(const TriangleObject& object);
    void setCamera(const vec3& origin, const vec3& target, int width, int height, float fov_y);

    void setEnvmap(const std::string filename);
    void setEnvmap(const std::string filename, int width, int height);

    EnvmapLoader *getEnvmap() const;
    void sampleLightAreaAndDirection(vec3& position, vec3& normal, float& pdfA, vec3& direction, float& pdfW, const TriangleObject *& obj) const; // assuming diffuse emission
    vec3 sampleLightArea(const vec3& position, /*const vec3& normal,*/ float& pdfW, int& light_id, vec3 *sample_pos = nullptr, vec3 *sample_normal = nullptr) const;
    float getLightArea() const;
    float lightPdfW(const vec3& pos, const vec3& nl, const vec3& dir) const;
    vec3 lightRadiance(const vec3& pos, /*const vec3& nl,*/ const vec3& dir, const int light_id = -1) const;

    HitInfo intersect(const Ray& ray) const;
    bool occluded(const Ray& ray, float tfar = NUM_INFINITY) const;

    TriangleObject& getTriangle(size_t n);
    const TriangleObject& getTriangle(size_t n) const;
    float getTriangleArea(size_t n) const;
    vec3 getCenter() const;
    size_t numTrangles() const;
    CCamera *getCamera() const;
    const SceneSphere& sceneSphere() const;
    void sampleLight(float& pdfLight, int& light_id, vec3& sample_pos, vec3& sample_normal, float rnd0, float rnd1, float rnd2) const;
};

#endif // scene_h__