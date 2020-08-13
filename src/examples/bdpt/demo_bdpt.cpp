#include <omp.h>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "bdpt_impl.h"
#include "mesh.h"
#include "scene.h"
#include "timer.h"
#include "sample.h"
#include "logger.h"

// TODO : involve envmap lighting

#define INPUT_DIR "../input/"

class DemoBDPT
{
public:
    void run();
};

///
static void loadModel0(Scene& scene);

#define PixelWidth 512
#define PixelHeight 512
ReferenceBDPT* bdpt = nullptr;

//////////////////////////////////////////////////////////////////////////

static void addSphere(Scene&       scene,
                      const vec3&  position,
                      float        radius,
                      const vec3&  color,
                      MaterialType mat_type)
{
    Mesh obj("unitsph.obj");
    //     CMesh obj("sphere.obj");
    {
        glm::mat4 translation = glm::translate(
            glm::mat4(1.0f), glm::vec3(position.x, position.y, position.z));
        glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(radius));

        for (int n = 0; n < obj.m_verts.size(); n++)
        {
            glm::vec4 v = glm::vec4(
                obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
            v = translation * scaling * v;
            obj.m_verts[n] = vec3(v.x, v.y, v.z);
        }
        obj.genNormals();
    }

    int id = todo_addMaterial(Material(color, 0, mat_type));

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

        scene.add(tri);
    }
}

static void addBunny(Scene& scene)
{
    Mesh obj("bunny_1k.obj");
    {
        //             vec3 translation = vec3(2, 0, 7);
        vec3 translation = vec3(2, 4, 7);  // +++++++++++++++
        //             vec3 translation = vec3(2, 4, 2);
        float scaling = 5;

        for (int n = 0; n < obj.m_verts.size(); n++)
        {
            vec3 v = vec3(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z);
            obj.m_verts[n] = vec3(v.x, v.y, v.z) * scaling + translation;
        }
        obj.genNormals();
    }

    vec3 emission = vec3(10);
    vec3 color = vec3(10);

    int id = todo_addMaterial(Material(color, 0, LGHT));

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

        scene.add(tri);
    }
}

static void loadModel0(Scene& scene, FrameBuffer& film)
{
    film.resize(PixelWidth, PixelHeight);

    float scale = 0.1;
    scene.setCamera(vec3(50.0, 40.8, 220.0) * scale,
                    vec3(50.0, 40.8, 0.0) * scale,
                    film.width(),
                    film.height(),
                    40.0f);

    // light
    addSphere(scene,
              vec3(10, 70, 51.6) * scale,
              6.0 * scale,
              vec3(100., 100., 100.),
              LGHT);
//     addSphere(scene, vec3(10, 70, 51.6) * scale, 10.0 * scale, vec3(1.0,
//     1.0, 1.0), LGHT);

    // walls
#if 0
    addSphere(scene, vec3(1e5 + 1, 40.8, 81.6), 1e5, vec3(0.75, 0.25, 0.25), DIFF);
    addSphere(scene, vec3(-1e5 + 99, 40.8, 81.6), 1e5, vec3(0.25, 0.25, 0.75), DIFF);
    addSphere(scene, vec3(50, 40.8, 1e5), 1e5, vec3(0.75, 0.65, 0.75), DIFF);
    addSphere(scene, vec3(50, 40.8, -1e5 + 350), 1e5, vec3(0.50, 0.50, 0.50), DIFF);
    addSphere(scene, vec3(50, 1e5, 81.6), 1e5, vec3(0.65, 0.75, 0.75), DIFF);
    addSphere(scene, vec3(50, -1e5 + 81.6, 81.6), 1e5, vec3(0.75, 0.75, 0.65), DIFF);
#else
    //     addSphere(scene, vec3(1e5 + 1 - 2e5, 40.8, 81.6) * scale, 1e5 *
    //     scale, vec3(0.75, 0.25, 0.25), DIFF); addSphere(scene, vec3(-1e5 + 99
    //     + 2e5, 40.8, 81.6) * scale, 1e5 * scale, vec3(0.25, 0.25, 0.75),
    //     DIFF); addSphere(scene, vec3(50, 40.8, 1e5 - 2e5) * scale, 1e5 *
    //     scale, vec3(0.75, 0.65, 0.75), DIFF); addSphere(scene, vec3(50, 40.8,
    //     -1e5 + 350 + 2e5) * scale, 1e5 * scale, vec3(0.50, 0.50, 0.50),
    //     DIFF);
    addSphere(scene,
              vec3(50, 1e5 - 2e5, 81.6) * scale,
              1e5 * scale,
              vec3(0.65, 0.75, 0.75),
              DIFF);  // floor
//     addSphere(scene, vec3(50, 1e5 - 2e5, 81.6) * scale, 1e5 * scale,
//     vec3(1.0, 1.0, 1.0), DIFF); // floor addSphere(scene, vec3(50, -1e5
//     + 81.6 + 2e5, 81.6) * scale, 1e5 * scale, vec3(0.75, 0.75, 0.65), DIFF);
#endif

    // objects
    addSphere(scene,
              vec3(50, 20, 50) * scale,
              20 * scale,
              vec3(0.25, 0.75, 0.25),
              REFR);
    addSphere(scene,
              vec3(19, 16.5, 25) * scale,
              16.5 * scale,
              vec3(0.99, 0.99, 0.99),
              SPEC);
    //     addSphere(scene, vec3(50, 20, 50) * scale, 20 * scale, vec3(1.0,
    //     1.0, 1.0), DIFF); addSphere(scene, vec3(19, 16.5, 25) * scale, 16.5 *
    //     scale, vec3(1.0, 1.0, 1.0), SPEC); addSphere(scene, vec3(77, 16.5,
    //     78) * scale, 16.5 * scale, vec3(0.99, 0.99, 0.99), REFR);

    addBunny(scene);
}

static void loadModel1(Scene& scene, FrameBuffer& film)
{
    film.resize(PixelWidth, PixelHeight);

    scene.setCamera(vec3(3.137931, 1.571292, -1.968698),
                    vec3(0.105605, 0.321292, -0.217984),
                    film.width(),
                    film.height(),
                    40.0f);

    Mesh obj("Armadillo.obj");
    {
        //         Rotation rx(-90, 1, 0, 0);
        glm::mat4 rotation = glm::rotate(
            glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
        //         Rotation ry(300, 0, 1, 0);
        //         Scaling sl(0.015);
        for (int n = 0; n < obj.m_verts.size(); n++)
        {
            glm::vec4 v = glm::vec4(
                obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
            v = scaling * rotation * v;
            obj.m_verts[n] = vec3(v.x, v.y, v.z);
        }
        obj.genNormals();
    }

    int id = todo_addMaterial(Material(vec3(0.7, 0.6, 0.5), 0, DIFF));

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

        scene.add(tri);
    }

    // light
    addSphere(scene, vec3(1, 0, 0), 0.2, vec3(100., 100., 100.), LGHT);
}

void loadScenel2(Scene& scene, FrameBuffer& film)
{
    film.resize(PixelWidth, PixelHeight);

    float scale = 0.1;
    scene.setCamera(vec3(278, 273, -800) * scale,
                    vec3(278, 273, -799) * scale,
                    film.width(),
                    film.height(),
                    39.3077);

    auto add = [&scene, &scale](const vec3&     a,
                                const vec3&     b,
                                const vec3&     c,
                                const Material& mat) -> void {
        vec3           n = normalize(cross(b - a, c - a));
        TriangleObject tri(a * scale, b * scale, c * scale);
        int            id = todo_addMaterial(mat);
        tri.setMaterialID(id);
        scene.add(tri);
    };

    vec3  white(1, 1, 1);
    vec3  red(.75, .25, .25);
    vec3  green(.25, .75, .25);
    vec3  box(1, 1, 1);
    float brightness = 10;  // 20;

    //////////////////////////////////////////////////////////////////////////
#if 1
    vec3 cbox_luminaire[] = {vec3(343, 548.79999 - 0.5, 227),
                             vec3(343, 548.79999 - 0.5, 332),
                             vec3(213, 548.79999 - 0.5, 332),
                             vec3(213, 548.79999 - 0.5, 227)};
    //     vec3 cbox_luminaire[] = {
    //         vec3(343, 548.79999 - 0.5 - 100, 227),
    //         vec3(343, 548.79999 - 0.5 - 100, 332),
    //         vec3(213, 548.79999 - 0.5 - 100, 332),
    //         vec3(213, 548.79999 - 0.5 - 100, 227) };

    vec3 cbox_luminaire2[] = {vec3(343, 548.79999 - 0.5, 227),
                              vec3(343, 548.79999 - 0.5, 332),
                              vec3(213, 548.79999 - 0.5, 332),
                              vec3(213, 548.79999 - 0.5, 227),
                              (cbox_luminaire[0] + cbox_luminaire[1]) * 0.5f,
                              (cbox_luminaire[0] + cbox_luminaire[2]) * 0.5f,
                              (cbox_luminaire[1] + cbox_luminaire[2]) * 0.5f,
                              (cbox_luminaire[0] + cbox_luminaire[3]) * 0.5f,
                              (cbox_luminaire[2] + cbox_luminaire[3]) * 0.5f};

    // smaller triangles
#if 0
    add(cbox_luminaire2[0], cbox_luminaire2[4], cbox_luminaire2[5], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[4], cbox_luminaire2[1], cbox_luminaire2[6], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[5], cbox_luminaire2[6], cbox_luminaire2[2], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[4], cbox_luminaire2[6], cbox_luminaire2[5], Material(vec3(brightness), 0, LGHT));
#else
    add(cbox_luminaire[0],
        cbox_luminaire[1],
        cbox_luminaire[2],
        Material(vec3(brightness), 0, LGHT));
//     add(cbox_luminaire[0], cbox_luminaire[2], cbox_luminaire[1],
//     Material(vec3(brightness), 0, LGHT));
#endif

#if 0
    add(cbox_luminaire2[0], cbox_luminaire2[5], cbox_luminaire2[7], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[7], cbox_luminaire2[8], cbox_luminaire2[3], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[5], cbox_luminaire2[2], cbox_luminaire2[8], Material(vec3(brightness), 0, LGHT));
    add(cbox_luminaire2[7], cbox_luminaire2[5], cbox_luminaire2[8], Material(vec3(brightness), 0, LGHT));
#else
    // larger triangle
    add(cbox_luminaire[0],
        cbox_luminaire[2],
        cbox_luminaire[3],
        Material(vec3(brightness), 0, LGHT));
//     add(cbox_luminaire[0], cbox_luminaire[3], cbox_luminaire[2],
//     Material(vec3(brightness), 0, LGHT));
#endif

#endif
    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_redwall[] = {vec3(556, 0, 0),
                           vec3(556, 0, 559.20001),
                           vec3(556, 548.79999, 559.20001),
                           vec3(556, 548.79999, 0)};

    add(cbox_redwall[0],
        cbox_redwall[1],
        cbox_redwall[2],
        Material(red, 0, DIFF));
    add(cbox_redwall[0],
        cbox_redwall[2],
        cbox_redwall[3],
        Material(red, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_greenwall[] = {vec3(0, 0, 559.20001),
                             vec3(0, 0, 0),
                             vec3(0, 548.79999, 0),
                             vec3(0, 548.79999, 559.20001)};

    add(cbox_greenwall[0],
        cbox_greenwall[1],
        cbox_greenwall[2],
        Material(green, 0, DIFF));
    add(cbox_greenwall[0],
        cbox_greenwall[2],
        cbox_greenwall[3],
        Material(green, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_floor[] = {vec3(556, 0, 0),
                         vec3(0, 0, 0),
                         vec3(0, 0, 559.20001),
                         vec3(556, 0, 559.20001)};

    add(cbox_floor[0], cbox_floor[1], cbox_floor[2], Material(white, 0, DIFF));
    add(cbox_floor[0], cbox_floor[2], cbox_floor[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_back[] = {vec3(556, 0, 559.20001),
                        vec3(0, 0, 559.20001),
                        vec3(0, 548.79999, 559.20001),
                        vec3(556, 548.79999, 559.20001)};

    add(cbox_back[0], cbox_back[1], cbox_back[2], Material(white, 0, DIFF));
    add(cbox_back[0], cbox_back[2], cbox_back[3], Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    vec3 cbox_ceiling[] = {vec3(556, 548.79999, 0),
                           vec3(556, 548.79999, 559.20001),
                           vec3(0, 548.79999, 559.20001),
                           vec3(0, 548.79999, 0)};

    add(cbox_ceiling[0],
        cbox_ceiling[1],
        cbox_ceiling[2],
        Material(white, 0, DIFF));
    add(cbox_ceiling[0],
        cbox_ceiling[2],
        cbox_ceiling[3],
        Material(white, 0, DIFF));

    //////////////////////////////////////////////////////////////////////////
    MaterialType smallbox_mat = REFR;  //  REFR, DIFF
    vec3         cbox_smallbox[] = {vec3(130.000000, 165.000000, 65.000000),
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
                            vec3(130.000000, 0.000000, 65.000000)};

    add(cbox_smallbox[0],
        cbox_smallbox[1],
        cbox_smallbox[2],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[0],
        cbox_smallbox[2],
        cbox_smallbox[3],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[4],
        cbox_smallbox[5],
        cbox_smallbox[6],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[4],
        cbox_smallbox[6],
        cbox_smallbox[7],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[8],
        cbox_smallbox[9],
        cbox_smallbox[10],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[8],
        cbox_smallbox[10],
        cbox_smallbox[11],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[12],
        cbox_smallbox[13],
        cbox_smallbox[14],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[12],
        cbox_smallbox[14],
        cbox_smallbox[15],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[16],
        cbox_smallbox[17],
        cbox_smallbox[18],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[16],
        cbox_smallbox[18],
        cbox_smallbox[19],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[20],
        cbox_smallbox[21],
        cbox_smallbox[22],
        Material(box, 0, smallbox_mat));
    add(cbox_smallbox[20],
        cbox_smallbox[22],
        cbox_smallbox[23],
        Material(box, 0, smallbox_mat));

    //////////////////////////////////////////////////////////////////////////
#if 1
    MaterialType tallbox_mat = SPEC;  // SPEC, REFR, DIFF
    vec3         cbox_largebox[] = {vec3(423.000000, 330.000000, 247.000000),
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
                            vec3(423.000000, 0.000000, 247.000000)};

    add(cbox_largebox[0],
        cbox_largebox[1],
        cbox_largebox[2],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[0],
        cbox_largebox[2],
        cbox_largebox[3],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[4],
        cbox_largebox[5],
        cbox_largebox[6],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[4],
        cbox_largebox[6],
        cbox_largebox[7],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[8],
        cbox_largebox[9],
        cbox_largebox[10],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[8],
        cbox_largebox[10],
        cbox_largebox[11],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[12],
        cbox_largebox[13],
        cbox_largebox[14],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[12],
        cbox_largebox[14],
        cbox_largebox[15],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[16],
        cbox_largebox[17],
        cbox_largebox[18],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[16],
        cbox_largebox[18],
        cbox_largebox[19],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[20],
        cbox_largebox[21],
        cbox_largebox[22],
        Material(box, 0, tallbox_mat));
    add(cbox_largebox[20],
        cbox_largebox[22],
        cbox_largebox[23],
        Material(box, 0, tallbox_mat));

#else  // cylinder
    //////////////////////////////////////////////////////////////////////////
    //     vec3 upper_center(400,200,300);
    //     vec3 lower_center(400,  0,300);
    vec3 upper_center(400, 200, 300);
    vec3 lower_center(400, 100, 300);
    for (int n = 0; n < 50; n++)
    {
        vec3 a(lower_center +
               vec3(sin(n / 50.0 * 2 * M_PI), 0, -cos(n / 50.0 * 2 * M_PI)) *
                   100);
        vec3 b(upper_center +
               vec3(sin(n / 50.0 * 2 * M_PI), 0, -cos(n / 50.0 * 2 * M_PI)) *
                   100);
        vec3 c(upper_center + vec3(sin((n + 1) / 50.0 * 2 * M_PI),
                                   0,
                                   -cos((n + 1) / 50.0 * 2 * M_PI)) *
                                  100);
        vec3 d(lower_center + vec3(sin((n + 1) / 50.0 * 2 * M_PI),
                                   0,
                                   -cos((n + 1) / 50.0 * 2 * M_PI)) *
                                  100);
        add_triangle(TriangleObject(a, b, c), box, 0, REFR);
        add_triangle(TriangleObject(a, c, d), box, 0, REFR);
        add_triangle(TriangleObject(upper_center, c, b), box, 0, REFR);
        add_triangle(TriangleObject(lower_center, a, d), box, 0, REFR);
    }
        //     for(int n=0; n<50; n++){
        //         vec3 a(vec3(  0+300,200,300)+vec3(0,sin( n
        //         /50.0*2*M_PI),-cos( n   /50.0*2*M_PI))*100); vec3
        //         b(vec3(200+300,200,300)+vec3(0,sin( n   /50.0*2*M_PI),-cos( n
        //         /50.0*2*M_PI))*100); vec3
        //         c(vec3(200+300,200,300)+vec3(0,sin((n+1)/50.0*2*M_PI),-cos((n+1)/50.0*2*M_PI))*100);
        //         vec3 d(vec3(
        //         0+300,200,300)+vec3(0,sin((n+1)/50.0*2*M_PI),-cos((n+1)/50.0*2*M_PI))*100);
        //         add_triangle(TriangleObject(a,b,c),box,0,REFR);
        //         add_triangle(TriangleObject(a,c,d),box,0,REFR);
        //     }
#endif
}

// ref shading
static vec3 rayTrace(const Ray& primary_ray, const Scene& scene)
{
    HitInfo   hit;
    Ray       cr(primary_ray);
    vec3      final_radiance(0, 0, 0);
    vec3      throughput(1, 1, 1);
    const int maxDepth = 10;

    for (int n = 0; n <= maxDepth; n++)
    {
        hit = scene.intersect(cr);
        if (!hit)
        {
            //             if (n > 0) {
            //                 vec3 r;
            //                 scene.envmap()->sampleRadiance(&cr.d[0], &r[0]);
            //                 final_radiance += mult(throughput, r);
            //             }
            break;
        }
        //         if (dot(cr.d, n_hit) > 0)
        //         {
        //             n_hit = -n_hit;
        //         }

        const Material& mat =
            todo_getMaterial(hit.triangleObject()->materialID());

        // now only light can emit radiance
        if (mat.type() == LGHT)
        {
            if (dot(cr.dir, hit.shadingNormal()) < 0)
            {
                vec3 obj_color = mat.color();
                //             vec3 Lo = obj_color / NUM_PI;
                vec3 Lo = obj_color;
                final_radiance += throughput * Lo;
            }
            break;
        }

        vec3 p_hit = cr.at(hit.distance());
        vec3 n_hit = hit.shadingNormal();
        // vec3 nhit_g = hit.getObject()->getFaceNormal();

        // continue tracing
        vec3 total_weight;
        vec3 wo = mat.sample(-cr.dir, n_hit, randf(), randf(), total_weight);
        throughput *= total_weight;
        cr = Ray(p_hit, wo);
    }

    if (final_radiance != final_radiance)
    {
        final_radiance = vec3(0, 0, 0);
    }
    return final_radiance;
}

static void render_pt(const Scene& scene, FrameBuffer& film)
{
    const int spp = 1000;  // 10 * 50 * 20;// *5;// *20;// *100;

    //////////////////////////////////////////////////////////////////////////
    int   width = film.width();
    int   height = film.height();
    float invWidth = 1 / float(width), invHeight = 1 / float(height);

    Timer tm;

    for (int y = 0; y < height; ++y)
    {
        Logger::info() << Logger::carriage_return << " render " << y+1 << " / " << height << "     ";
#pragma omp parallel for schedule(dynamic)
        for (int x = 0; x < width; ++x)
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

                Ray ray = scene.camera()->makeRay((x + 0.5 + dx) * invWidth,
                                                  (y + 0.5 + dy) * invHeight);

                vec3 v = rayTrace(ray, scene);
                g += vec3(f_max(0.0f, v.x), f_max(0.0f, v.y), f_max(0.0f, v.z));
            }

            film.accumulatePixel(x, y, g);
        }
    }

    film.scale(1.0 / spp);

    film.dumpHDR("test_pt_ref.hdr");

    // film.tonemap_reinhard();
    film.tonemapGamma(2.2);

    float duration = tm.getTime();
    Logger::info() << "took < " << duration << " > second" << std::endl;

    film.dumpPPM("test_pt_ref.ppm");
}

const int G_SPP = 100;

static void render_bdpt(const Scene& scene, FrameBuffer& film)
{
    bdpt = new ReferenceBDPT(*scene.camera(), film, scene);

    Timer tm;

    // integration
    unsigned int total_samples =
        film.width() * film.height() * G_SPP;  // 100;// *100;// *1000; // 100000;
    int step = total_samples / 100;

    std::vector<StandardSampler> samplers(omp_get_max_threads());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < total_samples; i++)
    {
        bdpt->render(MaxEvents, samplers[omp_get_thread_num()]);

        if (0 == i % step)
        {
            Logger::info() << Logger::carriage_return << " " << i / step + 1 << " / 100  ";
        }
    }

    float brightness_scale =
        float(film.width()) * float(film.height()) / float(total_samples);

    if (0)
    {
        bdpt->dumpPyramid(brightness_scale);
    }

    film.scale(brightness_scale);

    // film.tonemapReinhard();
    film.tonemapGamma(2.2);

    float duration = tm.getTime();
    Logger::info() << "took < " << duration << " > second" << std::endl;

    film.dumpPPM("test_bdpt.ppm");
}
//////////////////////////////////////////////////////////////////////////

void DemoBDPT::run()
{
    Scene scene;
    // scene.setEnvmap(INPUT_DIR "envmap.bin");

    FrameBuffer film;

    //     loadModel1(scene, film);
    //     loadModel0(scene, film);
    loadScenel2(scene, film);

    scene.finalize();

    //     for (auto& t : scene.m_triangles)
    //     {
    //         if (t.getMaterial().getType() == LGHT)
    //         {
    //             cout << t.getMaterial().getType() << endl;
    //         }
    //     }

    Logger::info() << "init complete" << std::endl;

    render_bdpt(scene, film);
    //     render_pt(scene, film);
}

int runTest(int argc, char* argv[])
{
    Logger::info() << "#procs: " << omp_get_num_procs() << std::endl;
    Logger::info() << "#threads: " << omp_get_max_threads() << std::endl;
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));
    //     omp_set_num_threads(1); // single thread
    Logger::info() << "#procs: " << omp_get_num_procs() << std::endl;
    Logger::info() << "#threads: " << omp_get_max_threads() << std::endl;

    DemoBDPT* demo = new DemoBDPT;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char** argv)
{
    runTest(argc, argv);
}
