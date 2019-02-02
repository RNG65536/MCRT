// Disclaimer: This demo is adapted from smallvcm http://www.smallvcm.com/

#pragma once

#include <glm/gtc/matrix_transform.hpp>
#include "intersection.h"
#include "aabb.h"
#include "scene.h"
#include "scenebvh.h"
#include "vectors.h"
#include "mesh.h"
#include "numeric.h"
#include "camera.h"
#include "light.h"

typedef HitInfo HitInfo_debug;

#if 0
class API_Scene
{
    class Sphere_debug
    {
    public:
        float rad;
        vec3 p;

        Sphere_debug(float r_, vec3 p_)
            : rad(r_), p(p_)
        {
        }

        float intersect(const Ray &r) const
        {
            // returns distance
            vec3d op = vec3d(p) - vec3d(r.o);
            double t, b = dot(op, vec3d(r.d)), det = b * b - dot(op, op) + rad * rad;

            if (det < 0)
            {
                return 1e20;
            }
            else
            {
                det = sqrt(det);
            }

            return (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);
        }

        vec3 sample() const
        {
            return p + sampleSphereUniform(randf(), randf()) * (rad * 2);
        }
    };

    std::vector<Sphere_debug>       m_geometry;
    std::vector<int> material_ids;

public:

    void finalize(SceneSphere& _sceneSphere)
    {
        AABB scene_bbox;
        for (auto& s : m_geometry)
        {
            AABB tri_bbox = AABB(s.p - vec3(s.rad), s.p + vec3(s.rad));
            scene_bbox.enclose(tri_bbox);
        }
        _sceneSphere.m_cenetr = scene_bbox.getCenter();
        _sceneSphere.m_radius = scene_bbox.diagonal().length() * 0.5f;
        _sceneSphere.m_inv_radius_sqr = 1.0f / sq(_sceneSphere.m_radius);
    }

    HitInfo_debug intersect(const Ray& ray) const
    {
        float t;
        int id = -1;

        //ray-sphere intersect.
        int n = m_geometry.size();
        float d, inf = 1e20;
        t = inf;
        for (int i = 0; i < n; i++)
        {
            d = m_geometry[i].intersect(ray);
            if (d < t)
            {
                t = d;
                id = i;
            }
        }

        HitInfo_debug ret;
        ret.setDistance(t);
        ret.setShadingNormal(normalize(ray.at(t) - m_geometry[id].p));
        ret.setPrimitiveID(id);
        ret.setMaterialID(id >= 0 ? material_ids[id] : -1);
        return ret;
    }

    bool occluded(const Ray& ray, float tfar) const
    {
        //ray-sphere intersect.
        int n = m_geometry.size();
        for (int i = 0; i < n; i++)
        {
            float d = m_geometry[i].intersect(ray);
            if (d < tfar)
            {
                return true;
            }
        }
        return false;
    }

    void addSphere(float radius, const vec3& center, int mat_id)
    {
        m_geometry.push_back(Sphere_debug(radius, center));
        material_ids.push_back(mat_id);
    }
};
#else
class API_Scene
{
    SceneBVH *bvh = nullptr;
    std::vector<TriangleObject> m_triangles;
    std::vector<int> material_ids;

public:
    API_Scene()
    {

    }
    ~API_Scene()
    {
        delete bvh;
    }
    void finalize(SceneSphere& _sceneSphere)
    {
        bvh = new SceneBVH(m_triangles);

        AABB scene_bbox;
        for (auto& s : m_triangles)
        {
            AABB tri_bbox = s.boundingBox();
            scene_bbox.enclose(tri_bbox);
        }
        _sceneSphere.m_cenetr = scene_bbox.center();
        _sceneSphere.m_radius = scene_bbox.diagonal().length() * 0.5f;
        _sceneSphere.m_inv_radius_sqr = 1.0f / sq(_sceneSphere.m_radius);
    }

    HitInfo_debug intersect(const Ray& ray) const
    {
        HitInfo_debug ret = bvh->intersect(ray);
        int pid = ret.primitiveID();
        ret.setMaterialID(pid >= 0 ? material_ids[pid] : -1);
        return ret;
    }

    bool occluded(const Ray& ray, float tfar) const
    {
        return bvh->occluded(ray, tfar);
    }

    void addSphere(float radius, const vec3& center, int mat_id)
    {
        Mesh obj("unitsph.obj");
        {
            glm::mat4 translation = glm::translate(glm::mat4(1.0f), glm::vec3(center.x, center.y, center.z));
            glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(radius));

            for (int n = 0; n < obj.m_verts.size(); n++)
            {
                glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
                v = translation * scaling * v;
                obj.m_verts[n] = vec3(v.x, v.y, v.z);
            }
            obj.genNormals();
        }

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
            tri.setMaterialID(mat_id);
            m_triangles.push_back(tri);
            material_ids.push_back(mat_id);
        }
    }

    void addTriangle(const vec3& a, const vec3& b, const vec3& c, float scale, int mat_id)
    {
        vec3 n = normalize(cross(b - a, c - a));
        TriangleObject tri(a * scale, b * scale, c * scale, n, n, n);
        tri.setMaterialID(mat_id);
        m_triangles.push_back(tri);
        material_ids.push_back(mat_id);
    }
};
#endif

class Scene_debug
{
public:
    API_Scene                       m_scene;
    Camera                         *m_camera = nullptr;

    std::vector<Material_debug>     m_materials;
    std::unordered_map<int, int>    m_material2Light;

    SceneSphere                     m_sceneSphere;
    std::vector<AbstractLight*>     m_lights;
    BackgroundLight*                m_background = nullptr;

public:
    Scene_debug()
    {
        build_scene_1();

        std::cout << "lights : " << m_lights.size() << std::endl;

        BuildSceneSphere();
    }

    void addLight(const vec3& a, const vec3& b, const vec3& c, float scale, int mat_id, float brightness)
    {
        AreaLight *l = new AreaLight(a * scale, b * scale, c * scale);
        l->m_intensity = vec3(brightness);
        m_lights.push_back(l);
        m_material2Light.insert(std::make_pair(mat_id, (int)m_lights.size() - 1));
        m_scene.addTriangle(a, b, c, scale, mat_id);
    }

    void build_scene_0()
    {
        // care the material id
        m_scene.addSphere(6.0f, vec3(10, 70, 51.6f), 2);// light
        //         m_scene.addSphere(1e5, vec3(1e5 + 1, 40.8, 81.6), 2);// walls
        //         m_scene.addSphere(1e5, vec3(-1e5 + 99, 40.8, 81.6), 3);
        //         m_scene.addSphere(1e5, vec3(50, 40.8, 1e5), 4);
        //         m_scene.addSphere(1e5, vec3(50, 40.8, -1e5 + 350), 4);
        //         m_scene.addSphere(1e5, vec3(50, 1e5, 81.6), 4);
        //         m_scene.addSphere(1e5, vec3(50, -1e5 + 81.6, 81.6), 4);
        if (1)
        {
            m_scene.addSphere(1e5f, vec3(1e5f + 1 - 2e5f, 40.8f, 81.6f), 2);// walls
            m_scene.addSphere(1e5f, vec3(-1e5 + 99 + 2e5f, 40.8f, 81.6f), 3);
            m_scene.addSphere(1e5f, vec3(50, 40.8f, 1e5f - 2e5f), 4);
            m_scene.addSphere(1e5f, vec3(50, 40.8f, -1e5f + 350 + 2e5f), 4);
            m_scene.addSphere(1e5f, vec3(50, 1e5f - 2e5f, 81.6f), 4);
            m_scene.addSphere(1e5f, vec3(50, -1e5f + 81.6f + 2e5f, 81.6f), 4);
        }
        m_scene.addSphere(20, vec3(50, 20, 50), 5);// objects
        m_scene.addSphere(16.5f, vec3(19, 16.5f, 25), 6);
        m_scene.addSphere(16.5f, vec3(77, 16.5f, 78), 7);

        // Camera
        m_camera = new Camera(vec3(50.0f, 40.8f, 220.0f), vec3(50.0f, 40.8f, 0.0f), vec3(0, 1, 0), 512, 512, 40.0f);

        //////////////////////////////////////////////////////////////////////////
        // Materials
        Material_debug mat;
        // 0) light1, will only emit
        m_materials.push_back(mat);
        // 1) light2, will only emit
        m_materials.push_back(mat);

        // 2) red diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0.75, 0.25, 0.25);
        //         mat.m_phongReflectance = vec3(0.0f);
        //         mat.m_phongExponent = 90.f;
        m_materials.push_back(mat);

        // 3) blue diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0.25, 0.25, 0.75);
        m_materials.push_back(mat);

        // 4) white diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0.75, 0.75, 0.75);
        m_materials.push_back(mat);

        // 5) green diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0.25, 0.75, 0.25);
        m_materials.push_back(mat);

        // 6) mirror ball
        mat.Reset();
        mat.m_mirrorReflectance = vec3(1.0f);
        m_materials.push_back(mat);

        // 7) glass ball
        mat.Reset();
        mat.m_mirrorReflectance = vec3(1.0f);
        mat.m_IOR = 1.6f;
        m_materials.push_back(mat);

        // 8) spherical light
        mat.Reset();
        m_materials.push_back(mat);

        //////////////////////////////////////////////////////////////////////////
        // Lights
        if (1)
        {
            //             // With light box
            //             m_lights.resize(2);
            //             AreaLight *l = new AreaLight(m_geometry[0].sample(), m_geometry[0].sample(), m_geometry[0].sample());
            //             //l->mIntensity = Vec3f(0.95492965f);
            //             l->m_intensity = vec3(1.03329895614464f);
            //             m_lights[0] = l;
            //             m_material2Light.insert(std::make_pair(0, 0));
            // 
            //             l = new AreaLight(m_geometry[0].sample(), m_geometry[0].sample(), m_geometry[0].sample());
            //             //l->mIntensity = Vec3f(0.95492965f);
            //             l->m_intensity = vec3(1.03329895614464f);
            //             m_lights[1] = l;
            //             m_material2Light.insert(std::make_pair(1, 1));
        }

        if (0)
        {
            SphericalLight *l = new SphericalLight(vec3(10, 70, 51.6f), 6.0f);
            l->m_intensity = vec3(1, 1, 1) / NUM_PI * sq(l->m_radius);
            m_lights.push_back(l);
            m_material2Light.insert(std::make_pair(8, (int)m_lights.size() - 1));
        }

        if (0)
        {
            DirectionalLight *l = new DirectionalLight(vec3(-1.0f, 1.5f, -1.0f));
            l->m_intensity = vec3(0.5f, 0.2f, 0.0f) * 20.0f;
            m_lights.push_back(l);
        }

        if (1)
        {
            PointLight *l = new PointLight(vec3(50, 60, 51.6f));
            l->m_intensity = vec3(8000.0f * (0.25f / NUM_PI));
            m_lights.push_back(l);
        }

        if (0)
        {
            BackgroundLight *l = new BackgroundLight;
            l->m_scale = 1.0f;
            m_lights.push_back(l);
            m_background = l;
        }
    }

    void build_scene_1()
    {
        const float scale = 0.1f;
        const float brightness = 10;// 20; //

#if 1
        // care the material id
        vec3 cbox_luminaire[] = {
            vec3(343, 548.79999f - 0.5f, 227),
            vec3(343, 548.79999f - 0.5f, 332),
            vec3(213, 548.79999f - 0.5f, 332),
            vec3(213, 548.79999f - 0.5f, 227) };
         addLight(cbox_luminaire[0], cbox_luminaire[1], cbox_luminaire[2], scale, 8, brightness);
         addLight(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[3], scale, 9, brightness);

//         vec3 cbox_luminaire[] = {
//             vec3(343, 548.79999 - 0.5 - 100, 227),
//             vec3(343, 548.79999 - 0.5 - 100, 332),
//             vec3(213, 548.79999 - 0.5 - 100, 332),
//             vec3(213, 548.79999 - 0.5 - 100, 227) };
//         addLight(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[1], scale, 8, brightness);
//         addLight(cbox_luminaire[0], cbox_luminaire[3], cbox_luminaire[2], scale, 9, brightness);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_redwall[] = {
            vec3(556, 0, 0),
            vec3(556, 0, 559.20001f),
            vec3(556, 548.79999f, 559.20001f),
            vec3(556, 548.79999f, 0) };

        m_scene.addTriangle(cbox_redwall[0], cbox_redwall[1], cbox_redwall[2], scale, 1);
        m_scene.addTriangle(cbox_redwall[0], cbox_redwall[2], cbox_redwall[3], scale, 1);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_greenwall[] = {
            vec3(0, 0, 559.20001f),
            vec3(0, 0, 0),
            vec3(0, 548.79999f, 0),
            vec3(0, 548.79999f, 559.20001f) };

        m_scene.addTriangle(cbox_greenwall[0], cbox_greenwall[1], cbox_greenwall[2], scale, 2);
        m_scene.addTriangle(cbox_greenwall[0], cbox_greenwall[2], cbox_greenwall[3], scale, 2);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_floor[] = {
            vec3(556, 0, 0),
            vec3(0, 0, 0),
            vec3(0, 0, 559.20001f),
            vec3(556, 0, 559.20001f) };

        m_scene.addTriangle(cbox_floor[0], cbox_floor[1], cbox_floor[2], scale, 3);
        m_scene.addTriangle(cbox_floor[0], cbox_floor[2], cbox_floor[3], scale, 3);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_back[] = {
            vec3(556, 0, 559.20001f),
            vec3(0, 0, 559.20001f),
            vec3(0, 548.79999f, 559.20001f),
            vec3(556, 548.79999f, 559.20001f) };

        m_scene.addTriangle(cbox_back[0], cbox_back[1], cbox_back[2], scale, 4);
        m_scene.addTriangle(cbox_back[0], cbox_back[2], cbox_back[3], scale, 4);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_ceiling[] = {
            vec3(556, 548.79999f, 0),
            vec3(556, 548.79999f, 559.20001f),
            vec3(0, 548.79999f, 559.20001f),
            vec3(0, 548.79999f, 0) };

        m_scene.addTriangle(cbox_ceiling[0], cbox_ceiling[1], cbox_ceiling[2], scale, 5);
        m_scene.addTriangle(cbox_ceiling[0], cbox_ceiling[2], cbox_ceiling[3], scale, 5);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_smallbox[] = {
            vec3(130.000000, 165.000000, 65.000000),
            vec3(82.000000, 165.000000, 225.000000),
            vec3(240.000000, 165.000000, 272.000000),
            vec3(290.000000, 165.000000, 114.000000),
            vec3(290.000000, 0.000000, 114.000000),
            vec3(290.000000, 165.000000, 114.000000),
            vec3(240.000000, 165.000000, 272.000000),
            vec3(240.000000, 0.000000, 272.000000),
            vec3(130.000000, 0.000000, 65.000000),
            vec3(130.000000, 165.000000, 65.000000),
            vec3(290.000000, 165.000000, 114.000000),
            vec3(290.000000, 0.000000, 114.000000),
            vec3(82.000000, 0.000000, 225.000000),
            vec3(82.000000, 165.000000, 225.000000),
            vec3(130.000000, 165.000000, 65.000000),
            vec3(130.000000, 0.000000, 65.000000),
            vec3(240.000000, 0.000000, 272.000000),
            vec3(240.000000, 165.000000, 272.000000),
            vec3(82.000000, 165.000000, 225.000000),
            vec3(82.000000, 0.000000, 225.000000),
            vec3(290.000000, 0.000000, 114.000000),
            vec3(240.000000, 0.000000, 272.000000),
            vec3(82.000000, 0.000000, 225.000000),
            vec3(130.000000, 0.000000, 65.000000) };

        m_scene.addTriangle(cbox_smallbox[ 0], cbox_smallbox[ 1], cbox_smallbox[ 2], scale, 6);
        m_scene.addTriangle(cbox_smallbox[ 0], cbox_smallbox[ 2], cbox_smallbox[ 3], scale, 6);
        m_scene.addTriangle(cbox_smallbox[ 4], cbox_smallbox[ 5], cbox_smallbox[ 6], scale, 6);
        m_scene.addTriangle(cbox_smallbox[ 4], cbox_smallbox[ 6], cbox_smallbox[ 7], scale, 6);
        m_scene.addTriangle(cbox_smallbox[ 8], cbox_smallbox[ 9], cbox_smallbox[10], scale, 6);
        m_scene.addTriangle(cbox_smallbox[ 8], cbox_smallbox[10], cbox_smallbox[11], scale, 6);
        m_scene.addTriangle(cbox_smallbox[12], cbox_smallbox[13], cbox_smallbox[14], scale, 6);
        m_scene.addTriangle(cbox_smallbox[12], cbox_smallbox[14], cbox_smallbox[15], scale, 6);
        m_scene.addTriangle(cbox_smallbox[16], cbox_smallbox[17], cbox_smallbox[18], scale, 6);
        m_scene.addTriangle(cbox_smallbox[16], cbox_smallbox[18], cbox_smallbox[19], scale, 6);
        m_scene.addTriangle(cbox_smallbox[20], cbox_smallbox[21], cbox_smallbox[22], scale, 6);
        m_scene.addTriangle(cbox_smallbox[20], cbox_smallbox[22], cbox_smallbox[23], scale, 6);

        //////////////////////////////////////////////////////////////////////////
        vec3 cbox_largebox[] = {
            vec3(423.000000, 330.000000, 247.000000),
            vec3(265.000000, 330.000000, 296.000000),
            vec3(314.000000, 330.000000, 456.000000),
            vec3(472.000000, 330.000000, 406.000000),
            vec3(423.000000, 0.000000, 247.000000),
            vec3(423.000000, 330.000000, 247.000000),
            vec3(472.000000, 330.000000, 406.000000),
            vec3(472.000000, 0.000000, 406.000000),
            vec3(472.000000, 0.000000, 406.000000),
            vec3(472.000000, 330.000000, 406.000000),
            vec3(314.000000, 330.000000, 456.000000),
            vec3(314.000000, 0.000000, 456.000000),
            vec3(314.000000, 0.000000, 456.000000),
            vec3(314.000000, 330.000000, 456.000000),
            vec3(265.000000, 330.000000, 296.000000),
            vec3(265.000000, 0.000000, 296.000000),
            vec3(265.000000, 0.000000, 296.000000),
            vec3(265.000000, 330.000000, 296.000000),
            vec3(423.000000, 330.000000, 247.000000),
            vec3(423.000000, 0.000000, 247.000000),
            vec3(472.000000, 0.000000, 406.000000),
            vec3(314.000000, 0.000000, 456.000000),
            vec3(265.000000, 0.000000, 296.000000),
            vec3(423.000000, 0.000000, 247.000000) };

        m_scene.addTriangle(cbox_largebox[ 0], cbox_largebox[ 1], cbox_largebox[ 2], scale, 7);
        m_scene.addTriangle(cbox_largebox[ 0], cbox_largebox[ 2], cbox_largebox[ 3], scale, 7);
        m_scene.addTriangle(cbox_largebox[ 4], cbox_largebox[ 5], cbox_largebox[ 6], scale, 7);
        m_scene.addTriangle(cbox_largebox[ 4], cbox_largebox[ 6], cbox_largebox[ 7], scale, 7);
        m_scene.addTriangle(cbox_largebox[ 8], cbox_largebox[ 9], cbox_largebox[10], scale, 7);
        m_scene.addTriangle(cbox_largebox[ 8], cbox_largebox[10], cbox_largebox[11], scale, 7);
        m_scene.addTriangle(cbox_largebox[12], cbox_largebox[13], cbox_largebox[14], scale, 7);
        m_scene.addTriangle(cbox_largebox[12], cbox_largebox[14], cbox_largebox[15], scale, 7);
        m_scene.addTriangle(cbox_largebox[16], cbox_largebox[17], cbox_largebox[18], scale, 7);
        m_scene.addTriangle(cbox_largebox[16], cbox_largebox[18], cbox_largebox[19], scale, 7);
        m_scene.addTriangle(cbox_largebox[20], cbox_largebox[21], cbox_largebox[22], scale, 7);
        m_scene.addTriangle(cbox_largebox[20], cbox_largebox[22], cbox_largebox[23], scale, 7);
#endif

        //////////////////////////////////////////////////////////////////////////
        // Camera
        m_camera = new Camera(vec3(278, 273, -800) * scale, vec3(278, 273, -799) * scale, vec3(0, 1, 0), 512, 512, 39.3077f);

        //////////////////////////////////////////////////////////////////////////
        // Materials
        Material_debug mat;
        // 0) black diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0);
        m_materials.push_back(mat);

        // 1) red diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0.75, 0.25, 0.25);
        //         mat.m_phongReflectance = vec3(0.0f);
        //         mat.m_phongExponent = 90.f;
        m_materials.push_back(mat);

        // 2) green diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(0.25, 0.75, 0.25);
        m_materials.push_back(mat);

        // 3) white diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(1);
        m_materials.push_back(mat);

        // 4) white diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(1);
        m_materials.push_back(mat);

        // 5) white diffuse
        mat.Reset();
        mat.m_diffuseReflectance = vec3(1);
        m_materials.push_back(mat);

        // 6) glass ball
        mat.Reset();
        mat.m_mirrorReflectance = vec3(1.0f);
        mat.m_IOR = 1.33f; // 1.6f;
        m_materials.push_back(mat);

        // 7) mirror ball
        mat.Reset();
        mat.m_mirrorReflectance = vec3(1.0f);
        m_materials.push_back(mat);

        // 8) light, will only emit
        mat.Reset();
        m_materials.push_back(mat);

        // 9) light, will only emit
        mat.Reset();
        m_materials.push_back(mat);


        //////////////////////////////////////////////////////////////////////////
        // Lights
        if (0)
        {
            SphericalLight *l = new SphericalLight(vec3(10, 70, 51.6f) * scale, 6.0f);
            l->m_intensity = vec3(1, 1, 1) / NUM_PI * sq(l->m_radius);
            m_lights.push_back(l);
            m_material2Light.insert(std::make_pair(8, (int)m_lights.size() - 1));
        }

        if (0)
        {
            DirectionalLight *l = new DirectionalLight(vec3(-1.0f, 1.5f, -1.0f));
            l->m_intensity = vec3(0.5f, 0.2f, 0.0f) * 20.0f;
            m_lights.push_back(l);
        }

        if (0)
        {
            PointLight *l = new PointLight(vec3(300, 450, 2.5) * scale);
            l->m_intensity = vec3(8000.0f * (0.25f / NUM_PI));
            m_lights.push_back(l);
        }

        if (0)
        {
            BackgroundLight *l = new BackgroundLight;
            l->m_scale = 1.0f;
            m_lights.push_back(l);
            m_background = l;
        }
    }

    ~Scene_debug()
    {
        for (size_t i = 0; i < m_lights.size(); i++)
        {
            delete m_lights[i];
        }
    }

    int toLightID(int mat_id) const
    {
        std::unordered_map<int, int>::const_iterator it =
            m_material2Light.find(mat_id);

        if (it != m_material2Light.end())
        {
            return it->second;
        }
        else
        {
            return -1;
        }
    }

    HitInfo_debug intersect(const Ray& ray) const
    {
        HitInfo_debug ret = m_scene.intersect(ray);
        return ret;
    }

    bool Occluded(
        const vec3 &_point,
        const vec3 &_dir,
        float _TMax) const
    {
        //         return false;

        Ray ray;
        ray.orig = _point + _dir * NUM_EPS_RAY;
        ray.dir = _dir;
        float tfar = _TMax - NUM_EPS_RAY;

        return m_scene.occluded(ray, tfar);
    }

    const Material_debug& GetMaterial(const int _materialIdx) const
    {
        return m_materials[_materialIdx];
    }

    int GetMaterialCount() const
    {
        return (int)m_materials.size();
    }

    const AbstractLight* GetLightPtr(int _lightIdx) const
    {
        _lightIdx = std::min<int>(_lightIdx, static_cast<int>(m_lights.size()) - 1);
        return m_lights[_lightIdx];
    }

    int GetLightCount() const
    {
        return (int)m_lights.size();
    }

    const BackgroundLight* GetBackground() const
    {
        return m_background;
    }

    void BuildSceneSphere()
    {
        m_scene.finalize(m_sceneSphere);

        std::cout << "scene sphere : ( "
            << m_sceneSphere.m_cenetr.x << ", "
            << m_sceneSphere.m_cenetr.y << ", "
            << m_sceneSphere.m_cenetr.z << " ), "
            << m_sceneSphere.m_radius << std::endl;
    }
};

