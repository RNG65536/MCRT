#include <omp.h>
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "photonmap.h"
#include "scene.h"
#include "geometry.h"
#include "sample.h"
#include "numeric.h"
#include "intersection.h"
#include "ray.h"
#include "material.h"
#include "film.h"
#include "camera.h"
#include "timer.h"
#include "aabb.h"
#include "numeric.h"
#include "mesh.h"

#define INPUT_DIR "../input/"
class DemoPhotonMap
{
public:
    void run();
};

static void loadModel(Scene& scene_bvh);

class CPointLight
{
    float _cos_theta_min; // to sample solid angle uniformly, sample cos_theta uniformly,
                        //similar to uniform sphere and hemisphere sampling
    vec3 _nl;
    OrthonormalFrame _onb;

    vec3 _pos;
    vec3 rad_power;

public:

    CPointLight(const vec3& pos, const vec3& nl = vec3(0, 1, 0), float cos_theta_min = -1)
        : _pos(pos), _nl(normalize(nl)), _cos_theta_min(cos_theta_min), _onb(_nl)
    {
        rad_power = vec3(1, 1, 1) * (800/4) * 100;
    }

    vec3 directIllumination(const vec3& p)
    {
        float l2 = lengthSquared(_pos - p);
        return rad_power / l2 / (4 * NUM_PI);
    }

    Photon emit() const
    {
        float pdf;
        vec3 ss = sampleSolidAngleUniform(randf(), randf(), _cos_theta_min, &pdf);
        vec3 sample_dir = normalize(_onb.toWorld(ss));

        Photon p;
        p.pos[0] = _pos.x;
        p.pos[1] = _pos.y;
        p.pos[2] = _pos.z;
        p.dir[0] = sample_dir.x;
        p.dir[1] = sample_dir.y;
        p.dir[2] = sample_dir.z;
        p.power[0] = rad_power.x / pdf / (4 * NUM_PI);
        p.power[1] = rad_power.y / pdf / (4 * NUM_PI);
        p.power[2] = rad_power.z / pdf / (4 * NUM_PI);

        return p;
    }

    const vec3& position() const
    {
        return _pos;
    }
};

static std::unique_ptr<CPointLight> light;
static PhotonMap pmap(10000000);  // care for overflow if too many photons
const int total_photons = 10000000;
const int spp = 10;
const int win_width = 512;
const int win_height = 512;
AABB specular_obj_bbox;

static vec3 rayTrace(const Ray& primary_ray, const Scene& scene)
{
    HitInfo hit;
    Ray cr(primary_ray);
    vec3 final_radiance(0, 0, 0);
    vec3 throughput(1, 1, 1);
    const int maxDepth = 5;
    bool first_diffuse = 1;

    for (int n = 0; n < maxDepth; n++)
    {
        hit = scene.intersect(cr);
        if (!hit)
        {
            break;
        }

        vec3 p_hit = cr.at(hit.distance());
        vec3 nl(hit.shadingNormal()); // shading normal

        vec3 local_radiance_out(0, 0, 0);

        const TriangleObject *obj = hit.triangleObject();
        const Material& mat = todo_getMaterial(obj->materialID());
        MaterialType mat_type = mat.type();
        vec3 obj_color = mat.color();
        const vec3 objEmission = vec3(mat.emission());

        if (mat_type == DIFF)
        {
            // direct lighting
            Ray shadow_ray(p_hit, normalize(light->position() - p_hit));
            if (!scene.occluded(shadow_ray))
            {
                vec3 Li = light->directIllumination(p_hit);
                float cos_wo;
                vec3 brdf = mat.evaluate(-cr.dir, nl, shadow_ray.dir, &cos_wo);
                vec3 local_radiance_out = Li * (cos_wo * brdf);
                final_radiance += mult(throughput, local_radiance_out);
            }

            // only for caustics, which cannot be sampled otherwise
            {
                vec3 rad;
                pmap.irradiance_estimate(&rad.x, &p_hit.x, &nl.x, 0.3f, 40);
                vec3 brdf = mat.evaluate(vec3(0), nl, vec3(0));
                vec3 local_radiance_out = rad * brdf;
                final_radiance += mult(throughput, local_radiance_out);
            }
        }

        // continue tracing for indirect lighting
        vec3 total_weight;
        vec3 wo = mat.sample(-cr.dir, nl, randf(), randf(), total_weight);
        throughput *= total_weight;
        cr = Ray(p_hit, wo);
    }

    if (final_radiance != final_radiance) final_radiance = vec3(0, 0, 0);
    return final_radiance;
}

static void tracePhoton(const Photon& p, const Scene& scene)
{
    Ray cr = Ray(vec3(p.pos[0], p.pos[1], p.pos[2]), vec3(p.dir[0], p.dir[1], p.dir[2]));
    vec3 throughput(p.power[0], p.power[1], p.power[2]);

    for (int i = 0; i < 10; i++)
    {
        HitInfo hit = scene.intersect(cr);
        if (!hit)
        {
            break;
        }

        vec3 p_hit = cr.at(hit.distance());
        vec3 nl(hit.shadingNormal()); // shading normal

        const TriangleObject *obj = hit.triangleObject();
        const Material& mat = todo_getMaterial(obj->materialID());
        MaterialType mat_type = mat.type();
        vec3 obj_color = mat.color();

        if (mat_type == DIFF)
        {
            if (i > 0) // only for indirect lighting
            {
                pmap.store(&throughput.x, &p_hit.x, &cr.dir.x);
            }
            break; // only for caustics
        }

        vec3 bsdf_weight;
        vec3 wo = mat.sample(-cr.dir, nl, randf(), randf(), bsdf_weight);
        cr = Ray(p_hit, wo);

//         float prob = (throughput.x * 0.3f + throughput.y * 0.6f + throughput.z * 0.1f);
        float prob = 1;// (bsdf_weight.x * 0.3f + bsdf_weight.y * 0.6f + bsdf_weight.z * 0.1f);
        if (randf() < prob)
        {
//             throughput = mult(throughput, obj_color) / prob;
            throughput *= bsdf_weight / prob;
        }
        else
        {
            break;
        }
    }
}

static void lightPass(const Scene& scene)
{
    std::cout << "light pass begin" << std::endl;

    int n_step = total_photons / 100;

    for (int i = 0; i < total_photons; i++)
    {
        if (0 == i % n_step)
        {
            printf("\r%d / %d", i / n_step, 100);
        }

        Photon p = light->emit();
        tracePhoton(p, scene);
    }

    pmap.scale_photon_power(1.0f / total_photons);
    pmap.balance();

    std::cout << "light pass end" << std::endl;
}

static void renderPass(const Scene& scene)
{
    FrameBuffer fb;
    fb.resize(win_width, win_height);

    //////////////////////////////////////////////////////////////////////////
    float invWidth = 1 / float(win_width), invHeight = 1 / float(win_height);

    float zTrans = 7;
    float yTrans = 2.5;
    vec3 lookat(scene.boxCenter());
    vec3 rayorig = lookat + vec3(0, yTrans*0.8f, zTrans*0.5f) - vec3(0, 2.1f, 0); // no x-shift
    float viewdist = (lookat - rayorig).length();
    float _rot = -120 / 180.0f * NUM_PI;
    rayorig = vec3(-viewdist*sin(_rot), yTrans*0.5f, viewdist*cos(_rot)) + lookat;
    {
        printf("lookat: %f, %f, %f\n", lookat.x, lookat.y, lookat.z);
        printf("rayorig: %f, %f, %f\n", rayorig.x, rayorig.y, rayorig.z);
    }

    Camera cam(rayorig, lookat, vec3(0, 1, 0), win_width, win_height, 45.0f);

    Timer tm;

    for (int y = 0; y < win_height; ++y)
    {
        printf("\r render %d / %d  ", y + 1, win_height);
#pragma omp parallel for schedule(dynamic)
        for (int x = 0; x < win_width; ++x)
        {
            vec3 g;
            for (int s = 0; s < spp; ++s)
            {
                float dx = 0;
                float dy = 0;
                if (spp > 1)
                {
                    float r1 = 2 * randf();
                    float r2 = 2 * randf();
                    dx = r1 < 1 ? std::sqrt(r1) - 1 : 1 - std::sqrt(2 - r1);
                    dy = r2 < 1 ? std::sqrt(r2) - 1 : 1 - std::sqrt(2 - r2);
                }

                Ray ray = cam.makeRay((x + 0.5f + dx) * invWidth, (y + 0.5f + dy) * invHeight);
                vec3 v = rayTrace(ray, scene);
                g += vec3(f_max(0.0f, v.x), f_max(0.0f, v.y), f_max(0.0f, v.z));
            }

            fb.accumulatePixel(x, y, g);
        }
    }

    fb.scale(1.0f / spp);

    fb.dumpHDR("test.hdr");

//     fb.tonemapReinhard();
    fb.tonemapGamma(2.2f);

    float duration = tm.getTime();
    printf("took < %f > second\n", duration);

    fb.dumpPPM("test.ppm");
}

void DemoPhotonMap::run()
{
    Scene scene;
    //     envmap.open("white.bin"); // works only with random variations added
    scene.setEnvmap(INPUT_DIR "envmap.bin"); //++++++++++++++++++++++++++
    //     envmap.open("envmap2.bin");
    //     envmap.open("at_the_window2.bin");
    //     envmap.open("sun.bin");
    //     envmap.open("sun2.bin");    envmap.scale(5);
    //     envmap.scale(2.5);
    //     envmap.scale(3);

    loadModel(scene);
    scene.finalize();


    glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
    glm::vec4 v = glm::vec4(
//         10 * -167.492416,
//         10 * 168.008728,
//         10 * -66.623940, 1.0f);
        10 * -167.492416,
        10 * 68.008728,
        10 * -66.623940, 1.0f);
    v = scaling * rotation * v;
    vec3 light_pos = vec3(v.x, v.y, v.z);

//     vec3 d = scene.boundingSphere().m_cenetr - light_pos;
//     float cos_min = std::sqrt(1.0f - sq(scene.boundingSphere().m_radius) / d.lengthSquared());
    vec3 d = specular_obj_bbox.center() - light_pos;
    float cos_min = std::sqrt(1.0f - sq(specular_obj_bbox.diagonal().length() * 0.5f) / d.lengthSquared());
    light.reset(new CPointLight(light_pos, normalize(d), cos_min));

    printf("init complete\n");

    lightPass(scene);
    renderPass(scene);
}

static void loadModel(Scene& scene_bvh)
{
    std::vector<Mesh> obj_files;
    obj_files.push_back(Mesh(INPUT_DIR "0_plane.obj"));
//     obj_files.push_back(Mesh("1_diffuse.obj"));
    obj_files.push_back(Mesh(INPUT_DIR "1_bunny.obj"));
//     obj_files.push_back(Mesh("2_torus.obj"));
    obj_files.push_back(Mesh(INPUT_DIR "2_ring.obj"));

    obj_files[0].flipSide(); // the model's normal points inward
//     obj_files[1].flipSide(); // the model's normal points inward

    int id[3] =
    {
        todo_addMaterial(Material(vec3(0.5f, 0.5f, 0.5f), 0, DIFF)),
        todo_addMaterial(Material(vec3(0.8f, 0.7f, 0.6f), 0, REFR)),
//         todo_addMaterial(Material(vec3(1.0f, 1.0f, 1.0f), 0, REFR)),
        todo_addMaterial(Material(vec3(1.0f, 1.0f, 1.0f), 0, SPEC))
    };

    for (int j = 0; j < obj_files.size(); j++)
    {
        Mesh& obj = obj_files[j];
        {
            glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
            glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
            obj.transform(scaling * rotation);
            obj.genNormals();
        }

        bool specaular = todo_getMaterial(id[j]).type() == SPEC || todo_getMaterial(id[j]).type() == REFR;

        for (int i = 0; i < obj.m_faces.size(); i++){
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
            tri.setMaterialID(id[j]);

            scene_bvh.add(tri);

            if (specaular)
            {
                specular_obj_bbox.enclose(tri.boundingBox());
            }
        }
    }
}

int runTest(int argc, char *argv[])
{
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));
    //     omp_set_num_threads(1); // single thread
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());


    DemoPhotonMap *demo = new DemoPhotonMap;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char **argv)
{
    runTest(argc, argv);
}

