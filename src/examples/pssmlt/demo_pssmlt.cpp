#include <omp.h>
#include <iostream>
#include <memory>
#include <glm/gtc/matrix_transform.hpp>
#include "scene.h"
#include "film.h"
#include "numeric.h"
#include "material.h"
#include "mesh.h"
#include "camera.h"
#include "envmap.h"
#include "ray.h"
#include "geometry.h"
#include "intersection.h"
#include "timer.h"
#include "sample.h"

#define INPUT_DIR "../input/"

class DemoPSSMLT
{
public:
    void run();
};

static void loadSceneArmadillo(Scene& scene, FrameBuffer& film)
{
    film.resize(640, 480);

    //     add_CUDA(BoxObject(vec3(0.4,1.01,3.4)+world_shift,vec3(0.6,1.11,3.6)+world_shift),   
    //         vec3(1,1,1),   30,       DIFF);
    //     add_CUDA(BoxObject(vec3(0,0,3)+world_shift,vec3(1,1,4)+world_shift), vec3(1,1,1),   0,       REFR);
    //     add_CUDA(BoxObject(vec3(0,0,3)+world_shift,vec3(1,1,4)+world_shift), vec3(1,1,1),   0,       DIFF);

    //     objLoader obj("bunny_1k.obj");
    //     obj.set(4, 1,-0.28,4);
    //     objLoader obj("dragon_65k.obj");
    //     objLoader obj("..\\dragon_fix.obj");
    //     obj.set(.015, 0,0.6,4.5);

    glm::mat4 shift = glm::translate(glm::mat4(1.0f), glm::vec3(0.5f, 0.0f, 0.0f));
    scene.loadModel(INPUT_DIR "bunny_1k.obj", shift, Material(vec3(10, 10, 10), 0, LGHT));

    //     objLoader obj("buddha.obj"); //obj.set(1, 0,0.6,4.5);
    //     {
    //         Rotation rx(-90, 1, 0, 0);
    //         Rotation ry(120, 0, 1, 0);
    //         Scaling sl(1.2);
    //         for (int n = 0; n < obj.vertCount; n++) {
    //             obj.vertexArray[n] = sl.onPosition(ry.onPosition(rx.onPosition(obj.vertexArray[n])));
    //         }
    //         obj.genNorm();
    //     }

    glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
    scene.loadModel(INPUT_DIR "Armadillo.obj", scaling * rotation, Material(vec3(1, 1, 1), 0, DIFF));

    // camera
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
    scene.setCamera(rayorig, lookat, film.width(), film.height(), 45.0f);

    printf("load complete");
}

static void loadSceneCornellBox(Scene& scene, FrameBuffer& film)
{
    film.resize(512, 512);

    float scale = 0.1;
    scene.setCamera(vec3(278, 273, -800) * scale, vec3(278, 273, -799) * scale, film.width(), film.height(), 39.3077);

    auto add = [&scene, &scale](const vec3& a, const vec3& b, const vec3& c, const Material& mat) -> void
    {
        vec3 n = normalize(cross(b - a, c - a));
        TriangleObject tri(a * scale, b * scale, c * scale);
        int id = todo_addMaterial(mat);
        tri.setMaterialID(id);
        scene.add(tri);
    };

    vec3 white(1, 1, 1);
    vec3 red(.75, .25, .25);
    vec3 green(.25, .75, .25);
    vec3 box(1, 1, 1);
    float brightness = 10; // 20;

    //////////////////////////////////////////////////////////////////////////
#if 1
    vec3 cbox_luminaire[] = {
        vec3(343, 548.79999 - 0.5, 227),
        vec3(343, 548.79999 - 0.5, 332),
        vec3(213, 548.79999 - 0.5, 332),
        vec3(213, 548.79999 - 0.5, 227) };
    //     vec3 cbox_luminaire[] = {
    //         vec3(343, 548.79999 - 0.5 - 100, 227),
    //         vec3(343, 548.79999 - 0.5 - 100, 332),
    //         vec3(213, 548.79999 - 0.5 - 100, 332),
    //         vec3(213, 548.79999 - 0.5 - 100, 227) };

    vec3 cbox_luminaire2[] = {
        vec3(343, 548.79999 - 0.5, 227),
        vec3(343, 548.79999 - 0.5, 332),
        vec3(213, 548.79999 - 0.5, 332),
        vec3(213, 548.79999 - 0.5, 227),
        (cbox_luminaire[0] + cbox_luminaire[1]) * 0.5f,
        (cbox_luminaire[0] + cbox_luminaire[2]) * 0.5f,
        (cbox_luminaire[1] + cbox_luminaire[2]) * 0.5f,
        (cbox_luminaire[0] + cbox_luminaire[3]) * 0.5f,
        (cbox_luminaire[2] + cbox_luminaire[3]) * 0.5f };

    // smaller triangles
#if 0
    add(cbox_luminaire2[0], cbox_luminaire2[4], cbox_luminaire2[5], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[4], cbox_luminaire2[1], cbox_luminaire2[6], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[5], cbox_luminaire2[6], cbox_luminaire2[2], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[4], cbox_luminaire2[6], cbox_luminaire2[5], Material(vec3(brightness), 0, LGHT));
#else
    add(cbox_luminaire[0], cbox_luminaire[1], cbox_luminaire[2], Material(vec3(brightness), 0, LGHT));
    //     add(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[1], Material(vec3(brightness), 0, LGHT));
#endif

#if 0
    add(cbox_luminaire2[0], cbox_luminaire2[5], cbox_luminaire2[7], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[7], cbox_luminaire2[8], cbox_luminaire2[3], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[5], cbox_luminaire2[2], cbox_luminaire2[8], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[7], cbox_luminaire2[5], cbox_luminaire2[8], Material(vec3(brightness), 0, LGHT));
#else
    // larger triangle
    add(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[3], Material(vec3(brightness), 0, LGHT));
    //     add(cbox_luminaire[0], cbox_luminaire[3], cbox_luminaire[2], Material(vec3(brightness), 0, LGHT));
#endif

#endif
    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_redwall[] = {
        vec3(556, 0, 0),
        vec3(556, 0, 559.20001),
        vec3(556, 548.79999, 559.20001),
        vec3(556, 548.79999, 0) };

    add(cbox_redwall[0], cbox_redwall[1], cbox_redwall[2], Material(red, 0, DIFF));
    add(cbox_redwall[0], cbox_redwall[2], cbox_redwall[3], Material(red, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_greenwall[] = {
        vec3(0, 0, 559.20001),
        vec3(0, 0, 0),
        vec3(0, 548.79999, 0),
        vec3(0, 548.79999, 559.20001) };

    add(cbox_greenwall[0], cbox_greenwall[1], cbox_greenwall[2], Material(green, 0, DIFF));
    add(cbox_greenwall[0], cbox_greenwall[2], cbox_greenwall[3], Material(green, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_floor[] = {
        vec3(556, 0, 0),
        vec3(0, 0, 0),
        vec3(0, 0, 559.20001),
        vec3(556, 0, 559.20001) };

    add(cbox_floor[0], cbox_floor[1], cbox_floor[2], Material(white, 0, DIFF));
    add(cbox_floor[0], cbox_floor[2], cbox_floor[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_back[] = {
        vec3(556, 0, 559.20001),
        vec3(0, 0, 559.20001),
        vec3(0, 548.79999, 559.20001),
        vec3(556, 548.79999, 559.20001) };

    add(cbox_back[0], cbox_back[1], cbox_back[2], Material(white, 0, DIFF));
    add(cbox_back[0], cbox_back[2], cbox_back[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_ceiling[] = {
        vec3(556, 548.79999, 0),
        vec3(556, 548.79999, 559.20001),
        vec3(0, 548.79999, 559.20001),
        vec3(0, 548.79999, 0) };

    add(cbox_ceiling[0], cbox_ceiling[1], cbox_ceiling[2], Material(white, 0, DIFF));
    add(cbox_ceiling[0], cbox_ceiling[2], cbox_ceiling[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    MaterialType smallbox_mat = REFR; //  REFR, DIFF
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

    add(cbox_smallbox[0], cbox_smallbox[1], cbox_smallbox[2], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[0], cbox_smallbox[2], cbox_smallbox[3], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[4], cbox_smallbox[5], cbox_smallbox[6], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[4], cbox_smallbox[6], cbox_smallbox[7], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[8], cbox_smallbox[9], cbox_smallbox[10], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[8], cbox_smallbox[10], cbox_smallbox[11], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[12], cbox_smallbox[13], cbox_smallbox[14], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[12], cbox_smallbox[14], cbox_smallbox[15], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[16], cbox_smallbox[17], cbox_smallbox[18], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[16], cbox_smallbox[18], cbox_smallbox[19], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[20], cbox_smallbox[21], cbox_smallbox[22], Material(box, 0, smallbox_mat));
    add(cbox_smallbox[20], cbox_smallbox[22], cbox_smallbox[23], Material(box, 0, smallbox_mat));

    //////////////////////////////////////////////////////////////////////////
#if 1
    MaterialType tallbox_mat = SPEC; // SPEC, REFR, DIFF
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

    add(cbox_largebox[0], cbox_largebox[1], cbox_largebox[2], Material(box, 0, tallbox_mat));
    add(cbox_largebox[0], cbox_largebox[2], cbox_largebox[3], Material(box, 0, tallbox_mat));
    add(cbox_largebox[4], cbox_largebox[5], cbox_largebox[6], Material(box, 0, tallbox_mat));
    add(cbox_largebox[4], cbox_largebox[6], cbox_largebox[7], Material(box, 0, tallbox_mat));
    add(cbox_largebox[8], cbox_largebox[9], cbox_largebox[10], Material(box, 0, tallbox_mat));
    add(cbox_largebox[8], cbox_largebox[10], cbox_largebox[11], Material(box, 0, tallbox_mat));
    add(cbox_largebox[12], cbox_largebox[13], cbox_largebox[14], Material(box, 0, tallbox_mat));
    add(cbox_largebox[12], cbox_largebox[14], cbox_largebox[15], Material(box, 0, tallbox_mat));
    add(cbox_largebox[16], cbox_largebox[17], cbox_largebox[18], Material(box, 0, tallbox_mat));
    add(cbox_largebox[16], cbox_largebox[18], cbox_largebox[19], Material(box, 0, tallbox_mat));
    add(cbox_largebox[20], cbox_largebox[21], cbox_largebox[22], Material(box, 0, tallbox_mat));
    add(cbox_largebox[20], cbox_largebox[22], cbox_largebox[23], Material(box, 0, tallbox_mat));

#else //cylinder
    //////////////////////////////////////////////////////////////////////////
    //     vec3 upper_center(400,200,300);
    //     vec3 lower_center(400,  0,300);
    vec3 upper_center(400, 200, 300);
    vec3 lower_center(400, 100, 300);
    for (int n = 0; n < 50; n++){
        vec3 a(lower_center + vec3(sin(n / 50.0 * 2 * M_PI), 0, -cos(n / 50.0 * 2 * M_PI)) * 100);
        vec3 b(upper_center + vec3(sin(n / 50.0 * 2 * M_PI), 0, -cos(n / 50.0 * 2 * M_PI)) * 100);
        vec3 c(upper_center + vec3(sin((n + 1) / 50.0 * 2 * M_PI), 0, -cos((n + 1) / 50.0 * 2 * M_PI)) * 100);
        vec3 d(lower_center + vec3(sin((n + 1) / 50.0 * 2 * M_PI), 0, -cos((n + 1) / 50.0 * 2 * M_PI)) * 100);
        add_triangle(TriangleObject(a, b, c), box, 0, REFR);
        add_triangle(TriangleObject(a, c, d), box, 0, REFR);
        add_triangle(TriangleObject(upper_center, c, b), box, 0, REFR);
        add_triangle(TriangleObject(lower_center, a, d), box, 0, REFR);
    }
    //     for(int n=0; n<50; n++){
    //         vec3 a(vec3(  0+300,200,300)+vec3(0,sin( n   /50.0*2*M_PI),-cos( n   /50.0*2*M_PI))*100);
    //         vec3 b(vec3(200+300,200,300)+vec3(0,sin( n   /50.0*2*M_PI),-cos( n   /50.0*2*M_PI))*100);
    //         vec3 c(vec3(200+300,200,300)+vec3(0,sin((n+1)/50.0*2*M_PI),-cos((n+1)/50.0*2*M_PI))*100);
    //         vec3 d(vec3(  0+300,200,300)+vec3(0,sin((n+1)/50.0*2*M_PI),-cos((n+1)/50.0*2*M_PI))*100);
    //         add_triangle(TriangleObject(a,b,c),box,0,REFR);
    //         add_triangle(TriangleObject(a,c,d),box,0,REFR);
    //     }
#endif
}

static void loadCaustics(Scene& scene, FrameBuffer& film)
{
    film.resize(512, 512);

    float scale = 0.1;
    scene.setCamera(vec3(278, 273, -800) * scale, vec3(278, 273, -799) * scale, film.width(), film.height(), 39.3077);
//     scene.setCamera(vec3(278, 273, -100) * scale, vec3(278, 273, 0) * scale, film.width(), film.height(), 90);

    auto add = [&scene, &scale](const vec3& a, const vec3& b, const vec3& c, const Material& mat) -> void
    {
        vec3 n = normalize(cross(b - a, c - a));
        TriangleObject tri(a * scale, b * scale, c * scale);
        int id = todo_addMaterial(mat);
        tri.setMaterialID(id);
        scene.add(tri);
    };

    vec3 white(1, 1, 1);
    vec3 red(.75, .25, .25);
    vec3 green(.25, .75, .25);
    vec3 box(1, 1, 1);
//     float brightness = 500;
    float brightness = 100;

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_luminaire[] = {
        vec3(343 - 30 * 1, 548.79999 - 0.5, 227 + 30 * 1),
        vec3(343 - 30 * 1, 548.79999 - 0.5, 332 - 30 * 1),
        vec3(213 + 30 * 1, 548.79999 - 0.5, 332 - 30 * 1),
        vec3(213 + 30 * 1, 548.79999 - 0.5, 227 + 30 * 1) };
//     vec3 cbox_luminaire[] = {
//         vec3(343, 548.79999 - 0.5, 227),
//         vec3(343, 548.79999 - 0.5, 332),
//         vec3(213, 548.79999 - 0.5, 332),
//         vec3(213, 548.79999 - 0.5, 227) };
    add(cbox_luminaire[0], cbox_luminaire[1], cbox_luminaire[2], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[3], Material(vec3(brightness), 0, LGHT));

//     vec3 cbox_luminaire[] = {
//         vec3(343, 548.79999 - 0.5 - 100, 227),
//         vec3(343, 548.79999 - 0.5 - 100, 332),
//         vec3(213, 548.79999 - 0.5 - 100, 332),
//         vec3(213, 548.79999 - 0.5 - 100, 227) };
//     add(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[1], Material(vec3(brightness), 0, LGHT));
//     add(cbox_luminaire[0], cbox_luminaire[3], cbox_luminaire[2], Material(vec3(brightness), 0, LGHT));

    //////////////////////////////////////////////////////////////////////////
#if 0
    vec3 cbox_redwall[] = {
        vec3(556, 0, 0),
        vec3(556, 0, 559.20001),
        vec3(556, 548.79999, 559.20001),
        vec3(556, 548.79999, 0) };

    add(cbox_redwall[0], cbox_redwall[1], cbox_redwall[2], Material(red, 0, DIFF));
    add(cbox_redwall[0], cbox_redwall[2], cbox_redwall[3], Material(red, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_greenwall[] = {
        vec3(0, 0, 559.20001),
        vec3(0, 0, 0),
        vec3(0, 548.79999, 0),
        vec3(0, 548.79999, 559.20001) };

    add(cbox_greenwall[0], cbox_greenwall[1], cbox_greenwall[2], Material(green, 0, DIFF));
    add(cbox_greenwall[0], cbox_greenwall[2], cbox_greenwall[3], Material(green, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_floor[] = {
        vec3(556, 0, 0),
        vec3(0, 0, 0),
        vec3(0, 0, 559.20001),
        vec3(556, 0, 559.20001) };

    add(cbox_floor[0], cbox_floor[1], cbox_floor[2], Material(white, 0, DIFF));
    add(cbox_floor[0], cbox_floor[2], cbox_floor[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_back[] = {
        vec3(556, 0, 559.20001),
        vec3(0, 0, 559.20001),
        vec3(0, 548.79999, 559.20001),
        vec3(556, 548.79999, 559.20001) };

    add(cbox_back[0], cbox_back[1], cbox_back[2], Material(white, 0, DIFF));
    add(cbox_back[0], cbox_back[2], cbox_back[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_ceiling[] = {
        vec3(556, 548.79999, 0),
        vec3(556, 548.79999, 559.20001),
        vec3(0, 548.79999, 559.20001),
        vec3(0, 548.79999, 0) };

    add(cbox_ceiling[0], cbox_ceiling[1], cbox_ceiling[2], Material(white, 0, DIFF));
    add(cbox_ceiling[0], cbox_ceiling[2], cbox_ceiling[3], Material(white, 0, DIFF));
#endif

    //////////////////////////////////////////////////////////////////////////
    glm::mat4 t(1.0f); // right multiplied, just like OpenGL API does
//     t = glm::translate(t, glm::vec3(200.000000, 200.000000, 300.000000) * scale);
//     t = glm::scale(t, glm::vec3(6.0f));
    t = glm::translate(t, glm::vec3(200.000000, 240.000000, 300.000000) * scale);
    t = glm::scale(t, glm::vec3(10.0f));
    scene.loadModel(INPUT_DIR "unitsph.obj", t, Material(vec3(1), 0, REFR));
}

static float MIS_weight(float pdfA, float pdfB)
{
    return (pdfA) / (pdfA + pdfB);
}

class EnvmapSampler
{
    const Scene *pscene;

public:
    EnvmapSampler(const Scene* scene)
    {
        pscene = scene;
    }
    vec3 sampleEnvmapUsingLightPDF(const vec3& wi, const vec3& pos, const vec3& nl, float& pdf, vec3& dir,
        const Material& bsdf, std::shared_ptr<Sampler> sampler)
    {
        vec3 irr(0);
        pdf = 1; // using pdf = 0 causes problem due to NANs
        dir = vec3(0, 0, 0);
        DirLight dl = pscene->envmap()->sample(sampler->next(), sampler->next());
        vec3 lightDirection(dl.dir[0], dl.dir[1], dl.dir[2]);
        bool occluded = pscene->occluded(Ray(pos, lightDirection));
        if (!occluded)
        {
            float g = f_max(0, dot(lightDirection, nl)); // cosine factor for irradiance
            if (g > 0)
            {
                //float weightEnvmap = pscene->envmap()->getPdf(dl.dir) / (pscene->envmap()->getPdf(dl.dir) + g / M_PI);
                pdf = pscene->envmap()->getPdf(dl.dir);
                dir = lightDirection;
                irr += (g * bsdf.evaluate(wi, nl, lightDirection) / dl.pdf) * vec3(dl.rad[0], dl.rad[1], dl.rad[2]); // note the 1/M_PI in the BRDF!!!!!
            }
        }
        if (irr != irr)
        {
            irr = vec3(0);
        }
        return irr;
    }
    vec3 sampleEnvmapUsingSurfacePDF(const vec3& wi, const vec3& pos, const vec3& nl, float& pdf, vec3& dir,
        const Material& bsdf, std::shared_ptr<Sampler> sampler)
    {
        // NOW JUST INLINES ideal diffuse BSDF
        vec3 irr(0);
        pdf = 1; // using pdf = 0 causes problem due to NANs
        dir = vec3(0, 0, 0);

//         float r1 = 2 * NUM_PI * randf(), r2 = randf(), r2s = std::sqrt(r2);
//         vec3 w = nl, u = normalize(cross(vec3(.3, .4, .5), w)), v = cross(w, u);
//         vec3 wo = (cos(r1)*r2s*u + sin(r1)*r2s*v + std::sqrt(1 - r2)*w).norm();
        vec3 weight;
        vec3 wo = bsdf.sample(wi, nl, sampler->next(), sampler->next(), weight);

        bool occluded = pscene->occluded(Ray(pos, wo));
        if (!occluded)
        {
            float g = f_max(0, dot(wo, nl)); // cosine factor for irradiance
            if (g > 0)
            {
                //float weightBRDF = (g / M_PI) / (pscene->envmap()->getPdf(&d[0]) + g / M_PI);
                // pdf = g / M_PI;
                pdf = bsdf.pdfW(wi, nl, wo);
                dir = wo;
                vec3 rad;
                pscene->envmap()->sampleRadiance(&wo[0], &rad[0]);
                // BSDF("1/pi") * cos(theta) / pdf("cos(theta)/pi") = 1
                irr += rad;
                // irr += (g * bsdf.eval(wi, d, nl) / pdf) * rad;
            }
        }
        if (irr != irr)
        {
            irr = vec3(0);
        }
        return irr;
    }
};

enum EventTypes { DIFFUSE, SPECULAR };

static vec3 rayTrace(const Ray& primary_ray, const Scene& scene, Sampler& sampler)
{
    const bool passive = false;
    const bool medium = true;

    //////////////////////////////////////////////////////////////////////////
    HitInfo hit;
    Ray cr(primary_ray);
    vec3 final_radiance(0, 0, 0);
    vec3 throughput(1, 1, 1);
    const int maxDepth = 5;
    MaterialType last_mat;
    std::vector<EventTypes> path;

    float sigma_s = 1 * 0.01;
    float sigma_a = 0 * 0.01;
//     float sigma_s = 1 * 0.02;
//     float sigma_a = 0 * 0.02;
//     float sigma_s = 1 * 0.002;
//     float sigma_a = 0 * 0.002;
//     float sigma_s = 1 * 0.0001;
//     float sigma_a = 0 * 0.0001;
    float sigma_t = sigma_a + sigma_s;
    float albedo = sigma_s / sigma_t;
//     std::cout << albedo << std::endl;

    bool last_event_specular = false;

    int n;
    for (n = 0; n <= maxDepth; n++)
    {
        if (!medium)
        {
            hit = scene.intersect(cr);
            if (!hit)
            {
                if (scene.envmap())
                {
                    vec3 r;
                    scene.envmap()->sampleRadiance(&cr.dir[0], &r[0]);
                    final_radiance += throughput * r;
                }
                break;
            }
        }
        else
        {
            hit = scene.intersect(cr);
            float dist = -std::log(f_max(sampler.next(), FLT_MIN)) / sigma_t;

            bool to_scatter = false;
            if (!hit)
            {
                to_scatter = true;
            }
            else
            {
                const TriangleObject *obj = hit.triangleObject();
                const Material& mat = todo_getMaterial(obj->materialID());

                bool inside = (dot(cr.dir, hit.faceNormal()) > 0) && mat.type() != LGHT;

                if (!inside && dist < hit.distance())
                {
                    to_scatter = true;
                }
            }

            if (to_scatter)
            {
                //             vec3 dir = sampleHG(0.1f, sampler.next(), sampler.next());
                vec3 dir = sampleSphereUniform(sampler.next(), sampler.next());
                OrthonormalFrame onb(cr.dir);
                dir = normalize(onb.toWorld(dir));
                vec3 pos = cr.at(dist);
                cr = Ray(pos, dir);

                if (!passive)
                {
                    int light_id;
                    float pdf_light;
                    vec3 sample_pos, sample_normal;
                    scene.sampleLight(pdf_light, light_id, sample_pos, sample_normal,
                        sampler.next(), sampler.next(), sampler.next()); // pdfA
                    //                 scene.sampleLightArea(pos, pdf_light, light_id, &sample_pos, &sample_normal); // pdfW

                    vec3 dir_to_light = sample_pos - pos;
                    float r2 = dir_to_light.lengthSquared();
                    dir_to_light = normalize(dir_to_light);
                    float cos_on_light = -dot(dir_to_light, sample_normal);
                    if (pdf_light > 0 && cos_on_light > 0)
                    {
                        // attenuated
                        vec3 rad_light = scene.lightRadiance(pos, dir_to_light, light_id) * std::exp(-sigma_t * sqrt(r2));

                        pdf_light = __pdfAtoW(pdf_light, r2, cos_on_light);

                        // albedo = sigma_s / sigma_t, but note that in case of volume emission, sigma_s should be excluded
                        // 1/pi is isotropic phase function
                        vec3 rad = throughput * albedo * M_1_4PI * rad_light / pdf_light;
                        final_radiance += rad;
                    }
                }

                throughput *= albedo;
                last_event_specular = false;

                continue;
            }
        }

        const TriangleObject *obj = hit.triangleObject();
        const Material& mat = todo_getMaterial(obj->materialID());
        vec3 phit = cr.at(hit.distance());
//         vec3 nhit(hit.getShadingNormal());
        vec3 nhit(hit.faceNormal());
        vec3 nhit_g = obj->normal();

//         bool inside = false;
//         vec3 nhit0(nhit);
//         if (dot(cr.d, nhit) > 0)
//         {
//             nhit = -nhit;
//             inside = true;
//         }

        if (mat.type() == LGHT)
        {
            if (passive)
            {
                if (dot(cr.dir, nhit) < 0)
                {
                    final_radiance += throughput * mat.color();
                }
            }
            else
            {
                if ((n == 0 || last_event_specular) && dot(cr.dir, nhit) < 0)
                {
                    final_radiance += throughput * mat.color();
                }
            }
            break;
        }

        if (mat.type() == REFR || mat.type() == SPEC)
        {
            last_event_specular = true;
        }
        else
        {
            last_event_specular = false;
        }

        // NEE for diffuse surface here
        if (!passive && !last_event_specular)
        {
            int light_id;
            float pdf_light;
            vec3 sample_pos, sample_normal;
            scene.sampleLight(pdf_light, light_id, sample_pos, sample_normal,
                sampler.next(), sampler.next(), sampler.next()); // pdfA
//             vec3 dir_to_light = scene.sampleLightArea(phit, pdf_light, light_id, &sample_pos, &sample_normal); // pdfW

            vec3 dir_to_light = sample_pos - phit;
            float r2 = dir_to_light.lengthSquared();
            dir_to_light = normalize(dir_to_light);
            float cos_on_light = -dot(dir_to_light, sample_normal);

            if (pdf_light > 0 && cos_on_light > 0)
            {
                pdf_light = __pdfAtoW(pdf_light, r2, cos_on_light);

                // attenuated
                vec3 rad_light = scene.lightRadiance(phit, dir_to_light, light_id);
                if (medium)
                {
                    rad_light *= std::exp(-sigma_t * sqrt(r2));
                }

                // albedo = sigma_s / sigma_t, but note that in case of volume emission, sigma_s should be excluded
                // 1/pi is isotropic phase function
                vec3 weight = mat.evaluate(-cr.dir, nhit, dir_to_light);
                float cos_theta = f_max(0.0f, dot(nhit, dir_to_light));
                vec3 rad = throughput * weight * rad_light * (cos_theta / pdf_light);
                final_radiance += rad;
            }
        }

        MaterialType mat_type = mat.type();
        const vec3 objEmission = vec3(mat.emission());
        last_mat = mat_type;

        final_radiance += throughput * objEmission;

        vec3 weight;
        vec3 d = mat.sample(-cr.dir, nhit, sampler.next(), sampler.next(), weight);
        throughput *= weight;
        cr = Ray(phit, d);
    }

//     std::cout << str(throughput) << std::endl;
//     std::cout << n << std::endl;

    if (final_radiance != final_radiance) final_radiance = vec3(0, 0, 0);
    return final_radiance;
}

int sscount = 0;
PathSample samplePath(const Scene& scene, Sampler& sampler)
{
    Camera *cam = scene.camera();
    int width = cam->filmWidth();
    int height = cam->filmHeight();

    // sample camera
    float x = sampler.next();
    float y = sampler.next();

    Ray ray = cam->makeRay(x, y);

    // sample path space
    vec3 v = rayTrace(ray, scene, sampler);

    PathSample s;
    s.x = x * width;
    s.y = y * height;
    s.F = v;
    s.weight = width * height;
    return s;
}

void renderPT(FrameBuffer &fb, const int spp, const Scene& scene)
{
    std::vector<std::shared_ptr<StandardSampler>> sampler(omp_get_max_threads());
    for (int i = 0; i < sampler.size(); i++)
    {
        sampler[i] = std::make_shared<StandardSampler>();
    }
    const int width = fb.width();
    const int height = fb.height();

    for (int y = 0; y < height; ++y)
    {
        printf("\r render %d / %d  ", y + 1, height);
#pragma omp parallel for schedule(dynamic)
        for (int x = 0; x < width; ++x)
        {
            vec3 g;
            for (int s = 0; s < spp; ++s)
            {
                PathSample ss = samplePath(scene,
                    *sampler[omp_get_thread_num()]);
                fb.accumulatePixel(ss.x, ss.y, ss.F);
            }
        }
    }

    std::cout << sscount << " / " << sscount * 1.0f / spp / width << std::endl;

    fb.scale(1.0 / spp);
}

void renderMLT(FrameBuffer &fb, const int spp, const Scene& scene)
{
    const int width = fb.width();
    const int height = fb.height();

    int num_mutations = spp * width * height;

    int num_threads = omp_get_max_threads();
    std::vector<FrameBuffer> temp_fb(num_threads);
    for (auto& f : temp_fb)
    {
        f.resize(width, height);
    }

    std::vector<std::shared_ptr<MetropolisSampler>> sampler(num_threads);
    for (int i = 0; i < sampler.size(); i++)
    {
        sampler[i] = std::make_shared<MetropolisSampler>();
    }

#pragma omp parallel for
    for (int n = 0; n < num_threads; n++)
    {
        int tid = omp_get_thread_num();
        auto& buffer = temp_fb[tid];
        MetropolisSampler& tsampler = *sampler[tid];

        int num_seeds = width * height;
        std::vector<PathSample> seed_paths(num_seeds);

        float sum_I = 0;
        for (int i = 0; i < num_seeds; i++)
        {
            tsampler.reset(true);
            PathSample s = samplePath(scene, tsampler);
            tsampler.accept();

            sum_I += luminance(s.F);
            seed_paths[i] = s;
        }

        // select seed
        float rnd = tsampler.randf() * sum_I;
        int selected_path = 0;
        float cdf = 0;
        for (int i = 0; i < num_seeds; i++)
        {
            cdf += luminance(seed_paths[i].F);
            if (cdf >= rnd)
            {
                selected_path = i;
                break;
            }
        }

        const float b = sum_I / num_seeds;
        const float p_large = 0.2; // 0.6; // 0.5;
        const int M = num_mutations / num_threads;
        int num_accept = 0, num_reject = 0;

        PathSample current = seed_paths[selected_path];

        int M_step = M / 100;
        int progress = 0;
        for (int i = 0; i < M; i++)
        {
            //         if (0 == i % (100000))
            //         std::cout << i << " / " << M << " [ " << num_accept << ", " << num_reject << " ] " << std::endl;
            if ((i + 1) % M_step == 0)
            {
                progress++;
                std::cout << progress << "% " << "Accept: " << num_accept << " Reject: " << num_reject << " Acc rate: " << 100.0 * num_accept / (num_accept + num_reject) << "%" << std::endl;
            }

            bool large_step = tsampler.randf() < p_large;
            tsampler.reset(large_step);
            PathSample proposal = samplePath(scene, tsampler);

            float a = std::min(1.0f, luminance(proposal.F) / luminance(current.F));
            float weight_proposal = (a + float(large_step)) / (luminance(proposal.F) / b + p_large) / M;
            float weight_current = (1.0 - a) / (luminance(current.F) / b + p_large) / M;
            //         float weight_proposal = (a) / (luminance(proposal.F) / b) / M;
            //         float weight_current = (1.0 - a) / (luminance(current.F) / b) / M;
            //         float weight_proposal = (0) / (luminance(proposal.F) / b) / M;
            //         float weight_current = (1.0 - 0) / (luminance(current.F) / b) / M;
            if (weight_proposal != weight_proposal)
            {
                weight_proposal = 0;
            }
            if (weight_current != weight_current)
            {
                weight_current = 0;
            }

            buffer.accumulatePixel(current.x, current.y, current.F * (current.weight * weight_current / num_threads));
            buffer.accumulatePixel(proposal.x, proposal.y, proposal.F * (proposal.weight * weight_proposal / num_threads));

            if (tsampler.randf() < a)
            {
                current = proposal;
                num_accept++;
                tsampler.accept();
            }
            else
            {
                num_reject++;
                tsampler.reject();
            }
        }
    }

    for (auto& f : temp_fb)
    {
        fb.accumulateBuffer(f);
    }
}

void DemoPSSMLT::run()
{
    Scene scene;
    //     envmap.open("white.bin"); // works only with random variations added
//     scene.setEnvmap("envmap.bin"); //++++++++++++++++++++++++++
    //     envmap.open("envmap2.bin");
    //     envmap.open("at_the_window2.bin");
    //     envmap.open("sun.bin");
    //     envmap.open("sun2.bin");    envmap.scale(5);
    //     envmap.scale(2.5);
    //     envmap.scale(3);

    FrameBuffer fb;

//     loadSceneArmadillo(scene, fb);
    //loadSceneCornellBox(scene, fb);
    loadCaustics(scene, fb);
    scene.finalize();

    printf("init complete\n");

    const int spp = 400;// 00;// 2000;

    //////////////////////////////////////////////////////////////////////////

    Timer tm;

//     renderPT(fb, spp, scene);
    renderMLT(fb, spp, scene);

    fb.dumpHDR("test.hdr");

    //fb.tonemap_reinhard();
    fb.tonemapGamma(2.2);

    float duration = tm.getTime();
    printf("took < %f > second\n", duration);

    fb.dumpPPM("test_mlt.ppm");
}

int runTest(int argc, char *argv[])
{
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));
    //     omp_set_num_threads(1); // single thread
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());


    DemoPSSMLT *demo = new DemoPSSMLT;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char **argv)
{
    runTest(argc, argv);
}

