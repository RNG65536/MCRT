#include <omp.h>
#include <iostream>
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include "scene.h"
#include "film.h"
#include "material.h"
#include "scene.h"
#include "mesh.h"
#include "geometry.h"
#include "envmap.h"
#include "numeric.h"
#include "ray.h"
#include "camera.h"
#include "intersection.h"
#include "timer.h"

class DemoSimple
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
    loadModel(scene, "bunny_1k.obj", shift, Material(vec3(10, 10, 10), 0, LGHT));

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
    loadModel(scene, "Armadillo.obj", scaling * rotation, Material(vec3(1, 1, 1), 0, DIFF));

    // camera
    float zTrans = 7;
    float yTrans = 2.5;
    vec3 lookat(scene.getCenter());
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
        const Material& bsdf) {
        vec3 irr(0);
        pdf = 1; // using pdf = 0 causes problem due to NANs
        dir = vec3(0, 0, 0);
        DirLight dl = pscene->getEnvmap()->sample(randf(), randf());
        vec3 lightDirection(dl.dir[0], dl.dir[1], dl.dir[2]);
        bool occluded = pscene->occluded(Ray(pos, lightDirection));
        if (!occluded) {
            float g = f_max(0, dot(lightDirection, nl)); // cosine factor for irradiance
            if (g > 0) {
                //float weightEnvmap = pscene->getEnvmap()->getPdf(dl.dir) / (pscene->getEnvmap()->getPdf(dl.dir) + g / M_PI);
                pdf = pscene->getEnvmap()->getPdf(dl.dir);
                dir = lightDirection;
                irr += (g * bsdf.evaluate(wi, nl, lightDirection) / dl.pdf) * vec3(dl.rad[0], dl.rad[1], dl.rad[2]); // note the 1/M_PI in the BRDF!!!!!
            }
        }
        if (irr != irr) {
            irr = vec3(0);
        }
        return irr;
    }
    vec3 sampleEnvmapUsingSurfacePDF(const vec3& wi, const vec3& pos, const vec3& nl, float& pdf, vec3& dir,
        const Material& bsdf) {
        // NOW JUST INLINES ideal diffuse BSDF
        vec3 irr(0);
        pdf = 1; // using pdf = 0 causes problem due to NANs
        dir = vec3(0, 0, 0);
        float r1 = 2 * NUM_PI * randf(), r2 = randf(), r2s = std::sqrt(r2);
        vec3 w = nl, u = normalize(cross(vec3(.3, .4, .5), w)), v = cross(w, u);
        vec3 wo = (cos(r1)*r2s*u + sin(r1)*r2s*v + std::sqrt(1 - r2)*w).norm();
        bool occluded = pscene->occluded(Ray(pos, wo));
        if (!occluded) {
            float g = f_max(0, dot(wo, nl)); // cosine factor for irradiance
            if (g > 0) {
                //float weightBRDF = (g / M_PI) / (pscene->getEnvmap()->getPdf(&d[0]) + g / M_PI);
                // pdf = g / M_PI;
                pdf = bsdf.pdfW(wi, nl, wo);
                dir = wo;
                vec3 rad;
                pscene->getEnvmap()->sampleRadiance(&wo[0], &rad[0]);
                // BSDF("1/pi") * cos(theta) / pdf("cos(theta)/pi") = 1
                irr += rad;
                // irr += (g * bsdf.eval(wi, d, nl) / pdf) * rad;
            }
        }
        if (irr != irr) {
            irr = vec3(0);
        }
        return irr;
    }
};

static vec3 rayTrace(const Ray& primary_ray, const Scene& scene)
{
    bool passive = 0;

    HitInfo hit;
    Ray cr(primary_ray);
    vec3 final_radiance(0, 0, 0);
    vec3 throughput(1, 1, 1);
    const int maxDepth = 10;
    MaterialType last_mat;

    for (int n = 0; n <= maxDepth; n++)
    {
        hit = scene.intersect(cr);
        if (!hit)
        {
            if ((passive || n == 0) && scene.getEnvmap())
            {
                vec3 r;
                scene.getEnvmap()->sampleRadiance(&cr.d[0], &r[0]);
                final_radiance += mult(throughput, r);
            }
            break;
        }

        const TriangleObject *obj = hit.getObject();
        const Material& mat = todo_getMaterial(obj->materialID());
        vec3 phit = cr.at(hit.getDistance());
        vec3 nhit(hit.getShadingNormal());
//         bool inside = false;
//         vec3 nhit0(nhit);
//         if (dot(cr.d, nhit) > 0)
//         {
//             nhit = -nhit;
//             inside = true;
//         }

        if (mat.type() == LGHT)
        {
            if ((n == 0 || passive || last_mat != DIFF) && dot(cr.d, nhit) < 0)
            {
                final_radiance += throughput * mat.color();
            }
            break;
        }

        vec3 nhit_g = obj->getFaceNormal();
        vec3 local_radiance_out(0, 0, 0);
        MaterialType mat_type = mat.type();
        vec3 objColor = mat.color();
        const vec3 objEmission = vec3(mat.emission());
        last_mat = mat_type;

        if (!passive && mat_type == DIFF)
        {
            if (scene.getEnvmap())
            {
                EnvmapSampler rs(&scene);
                float envmapPDF, bsdfPDF;
                vec3 envmapDir, bsdfDir;
                vec3 rad_envmap = rs.sampleEnvmapUsingLightPDF(-cr.d, phit, nhit, envmapPDF, envmapDir, mat);
                vec3 rad_bsdf = rs.sampleEnvmapUsingSurfacePDF(-cr.d, phit, nhit, bsdfPDF, bsdfDir, mat);

                float weightEnvmap = MIS_weight(envmapPDF, mat.pdfW(-cr.d, nhit, envmapDir));
                float weightBSDF = MIS_weight(bsdfPDF, scene.getEnvmap()->getPdf(&bsdfDir[0]));

                vec3 rad = rad_envmap * weightEnvmap + rad_bsdf * weightBSDF;
                //             vec3 rad = rad_envmap;
                //             vec3 rad = rad_bsdf;

                local_radiance_out += mult(objColor, rad);
            }

            if (scene.getLightArea() > 0)
            {
                // BRDF sampling
                vec3 weight;
                float pdfW_BRDF;
                vec3 brdf_sampled_dir = mat.sample(-cr.d, nhit, randf(), randf(), weight, &pdfW_BRDF);

                float pdfW_light_passive = scene.lightPdfW(phit, nhit, brdf_sampled_dir);
                float mis_weight_brdf = MIS_weight(pdfW_BRDF, pdfW_light_passive);
                vec3 rad_brdf_dir = scene.lightRadiance(phit, /*nhit,*/ brdf_sampled_dir);

                if (pdfW_light_passive == 0 && rad_brdf_dir.length() > 0) // these are the cases that cause fireflies
                {
                    std::cout << "{ " << mis_weight_brdf << " }" << std::endl;
                }

                local_radiance_out += weight * (rad_brdf_dir) * (mis_weight_brdf);

                //////////////////////////////////////////////////////////////////////////

                //light sampling
                int light_id;
                float pdfW_light;
                vec3 light_sampled_dir = scene.sampleLightArea(phit, /*nhit,*/ pdfW_light, light_id);

                if (pdfW_light > 0)
                {
                    // light weighting
                    float pdfW_BRDF_passive = mat.pdfW(-cr.d, nhit, light_sampled_dir);
                    float mis_weight_light = MIS_weight(pdfW_light, pdfW_BRDF_passive);

                    vec3 rad_light_dir = scene.lightRadiance(phit, /*nhit,*/ light_sampled_dir, light_id);
                    vec3 weight = mat.evaluate(-cr.d, nhit, light_sampled_dir);
                    float cos_theta_light = std::max(dot(nhit, light_sampled_dir), 0.0f);
                    local_radiance_out += weight * (rad_light_dir) * (cos_theta_light / pdfW_light * mis_weight_light);
                }
            }
        }

        local_radiance_out += objEmission;
        final_radiance += mult(throughput, local_radiance_out);

        vec3 weight;
        vec3 d = mat.sample(-cr.d, nhit, randf(), randf(), weight);
        throughput *= weight;
        cr = Ray(phit, d);

    }

    if (final_radiance != final_radiance) final_radiance = vec3(0, 0, 0);
    return final_radiance;
}

void DemoSimple::run()
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
    loadSceneCornellBox(scene, fb);
    scene.finalize();

    printf("init complete\n");

    const int spp = 100;

    //////////////////////////////////////////////////////////////////////////
    CCamera *cam = scene.getCamera();

    WTimer tm;

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

                float dx = 0;
                float dy = 0;
                if (spp > 1)
                {
                    float r1 = 2 * randf();
                    float r2 = 2 * randf();
                    dx = r1 < 1 ? std::sqrt(r1) - 1 : 1 - std::sqrt(2 - r1);
                    dy = r2 < 1 ? std::sqrt(r2) - 1 : 1 - std::sqrt(2 - r2);
                }

                Ray ray = cam->makeRay((x + 0.5 + dx) / width, (y + 0.5 + dy) / height);

                vec3 v = rayTrace(ray, scene);
                g += vec3(f_max(0.0f, v.x), f_max(0.0f, v.y), f_max(0.0f, v.z));
            }

            fb.accumulatePixel(x, y, g);
        }
    }

    fb.scale(1.0 / spp);
    fb.dumpHDR("test.hdr");

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


    DemoSimple *demo = new DemoSimple;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char **argv)
{
    runTest(argc, argv);
}

