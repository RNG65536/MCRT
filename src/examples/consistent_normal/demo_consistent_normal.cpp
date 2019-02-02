#include <omp.h>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include "material.h"
#include "intersection.h"
#include "ray.h"
#include "scene.h"
#include "mesh.h"
#include "envmap.h"
#include "camera.h"
#include "numeric.h"
#include "film.h"
#include "timer.h"

#define INPUT_DIR "../input/"
class DemoConsistentNormal
{
public:
    void run();
};
static void loadModel(Scene& scene_CUDA)
{
    //     add_CUDA(BoxObject(vec3(0.4,1.01,3.4)+world_shift,vec3(0.6,1.11,3.6)+world_shift),   
    //         vec3(1,1,1),   30,       DIFF);
    //     add_CUDA(BoxObject(vec3(0,0,3)+world_shift,vec3(1,1,4)+world_shift), vec3(1,1,1),   0,       REFR);
    //     add_CUDA(BoxObject(vec3(0,0,3)+world_shift,vec3(1,1,4)+world_shift), vec3(1,1,1),   0,       DIFF);

//     objLoader obj("bunny_1k.obj");
//     obj.set(4, 1, -0.28, 4);
//     objLoader obj("dragon_65k.obj");
//     objLoader obj("..\\dragon_fix.obj");
//     obj.set(.015, 0, 0.6, 4.5);


     Mesh obj(INPUT_DIR "buddha.obj"); //obj.set(1, 0,0.6,4.5);
     {
         // Rotation rx(-90, 1, 0, 0);
         glm::mat4 rx = glm::rotate(glm::mat4(1.0f), glm::radians(270.0f), glm::vec3(1.0f, 0.0f, 0.0f));
         // Rotation ry(120, 0, 1, 0);
         glm::mat4 ry = glm::rotate(glm::mat4(1.0f), glm::radians(120.0f), glm::vec3(0.0f, 1.0f, 0.0f));
         // Scaling sl(1.2);
         glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(1.2f));
         for (int n = 0; n < obj.m_verts.size(); n++)
         {
             glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
             v = scaling * ry * rx * v;
             obj.m_verts[n] = vec3(v.x, v.y, v.z);
             //             obj.vertexArray[n] = sl.onPosition(ry.onPosition(rx.onPosition(obj.vertexArray[n])));
         }
         obj.genNormals();
     }

//     CMesh obj("Armadillo.obj");
//     CMesh obj("arma.obj");
//     {
//         glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//         glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
//         for (int n = 0; n < obj.m_verts.size(); n++) {
//             glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
//             v = scaling * rotation * v;
//             obj.m_verts[n] = vec3(v.x, v.y, v.z);
//         }
//         obj.genNormals();
//     }

//     CMesh obj("unitsph.obj");
//     CMesh obj("sphere.obj");
//     CMesh obj("finesphere.obj");
//     CMesh obj("fineplane.obj");
//     {
//         //         Rotation rx(-90, 1, 0, 0);
//         glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//         glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.9f));
//         //         Rotation ry(300, 0, 1, 0);
//         //         Scaling sl(0.015);
//         for (int n = 0; n < obj.m_verts.size(); n++) {
//             glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
//             v = scaling * rotation * v;
//             obj.m_verts[n] = vec3(v.x, v.y, v.z);
//         }
//         obj.genNormals();
//     }

    int id = todo_addMaterial(Material(vec3(0.9, 0.9, 0.9), 0, SPEC));
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

        tri.setVertexAlpha(vec3(obj.m_vertex_alpha[a],
            obj.m_vertex_alpha[b], obj.m_vertex_alpha[c]));

        scene_CUDA.add(tri);
    }

    printf("load complete");
}

// ref shading
static vec3 rayTrace(const Ray& primary_ray, const Scene& scene)
{
    HitInfo hit;
    Ray cr(primary_ray);
    vec3 radianceIntegral(0, 0, 0);
    vec3 throughput(1, 1, 1);
    const int maxDepth = 10;

    for (int n = 0; n <= maxDepth; n++)
    {
        hit = scene.intersect(cr);
        if (!hit)
        {
//             if (n > 0)
            {
                vec3 r;
                scene.envmap()->sampleRadiance(&cr.dir[0], &r[0]);
                radianceIntegral += mult(throughput, r);
            }
            break;
        }

        vec3 phit = cr.at(hit.distance());
        vec3 nhit(hit.shadingNormal());
//         vec3 nhit_g = obj->getFaceNormal();
//         nhit = hit.getConsistentNormal(-cr.d);

        const Material& mat = todo_getMaterial(hit.triangleObject()->materialID());
        vec3 weight;
        vec3 d = mat.sample(-cr.dir, nhit, randf(), randf(), weight);
        vec3 nc = hit.consistentNormal(-cr.dir, &d); // consistent reflection

        cr = Ray(phit, d);

        const vec3 objEmission = vec3(mat.emission());
        radianceIntegral += mult(throughput, objEmission);
        throughput = mult(throughput, weight);
    }

    if (radianceIntegral != radianceIntegral) radianceIntegral = vec3(0, 0, 0);
    return radianceIntegral;
}

static vec3 rayTrace1(const Ray& primary_ray, const Scene& scene)
{
    Ray cr(primary_ray);
    HitInfo hit = scene.intersect(cr);

    if (!hit)
    {
        vec3 r;
        scene.envmap()->sampleRadiance(&cr.dir[0], &r[0]);
        return r;
    }

    vec3 phit = cr.at(hit.distance());
    vec3 nhit(hit.shadingNormal());
//     nhit = hit.getFaceNormal();
    nhit = hit.consistentNormal(-cr.dir);

//     return nhit;

    if (dot(cr.dir, nhit) >= 0)
    {
        return vec3(1, 0, 0);
    }
    else
    {
        return vec3(-dot(cr.dir, nhit));
    }
}

void DemoConsistentNormal::run()
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

    printf("init complete\n");

    const int win_width = 1024;
    const int win_height = 1024;
    const int spp = 10;

    FrameBuffer fb;
    fb.resize(win_width, win_height);

    //////////////////////////////////////////////////////////////////////////
    float invWidth = 1 / float(win_width), invHeight = 1 / float(win_height);

    float zTrans = 7;
    float yTrans = 2.5;
    vec3 lookat(scene.boxCenter());
    vec3 rayorig = lookat + vec3(0, yTrans*0.8, zTrans*0.5) - vec3(0, 2.1, 0); // no x-shift
    float viewdist = (lookat - rayorig).length();
    float _rot = -120 / 180.0 * NUM_PI;
    rayorig = vec3(-viewdist*sin(_rot), yTrans*0.5, viewdist*cos(_rot)) + lookat;
    {
        printf("lookat: %f, %f, %f\n", lookat.x, lookat.y, lookat.z);
        printf("rayorig: %f, %f, %f\n", rayorig.x, rayorig.y, rayorig.z);
    }

    //     glm::mat4 view = glm::lookAt(glm::vec3(rayorig.x, rayorig.y, rayorig.z), glm::vec3(lookat.x, lookat.y, lookat.z), glm::vec3(0.0f, 1.0f, 0.0f));
    //     glm::mat4 projection = glm::perspective(glm::radians(45.0f), float(win_width) / float(win_height), 0.01f, 1000.0f);
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

                Ray ray = cam.makeRay((x + 0.5 + dx) * invWidth, (y + 0.5 + dy) * invHeight);

                vec3 v = rayTrace(ray, scene);
                g += vec3(f_max(0.0f, v.x), f_max(0.0f, v.y), f_max(0.0f, v.z));
            }

            fb.accumulatePixel(x, y, g);
        }
    }

    fb.scale(1.0 / spp);

    //fb.tonemap_reinhard();
    fb.tonemapGamma(2.2);

    float duration = tm.getTime();
    printf("took < %f > second\n", duration);

    fb.dumpPPM("test.ppm");
}

int runTest(int argc, char *argv[])
{
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));
    //     omp_set_num_threads(1); // single thread
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());


    DemoConsistentNormal *demo = new DemoConsistentNormal;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char **argv)
{
    runTest(argc, argv);
}

