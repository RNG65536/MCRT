#include <omp.h>
#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include "numeric.h"
#include "octree.h"
#include "mesh_sample.h"
#include "bssrdf.h"
#include "ray.h"
#include "scene.h"
#include "mesh.h"
#include "material.h"
#include "consoledebug.h"
#include "envmap.h"
#include "intersection.h"
#include "sample.h"
#include "geometry.h"
#include "volumetric.h"
#include "film.h"
#include "camera.h"
#include "timer.h"

class DemoBSSRDF
{
public:
    void run();
};
// TODO : investigate whether the samples are chosen on the normalized surface or the original surface, and find out a suitable sampling efficiency

#define DIRPOLE 1

#define REBUILD 1
#define REBAKE 1

#define HIDE_BACKGROUND 1
#define G_TOTAL_SAMPLE  3e6  //1e6 // 
#define G_IRRADIANCE_SAMPLE 64 // 256

//#define G_MIN_DIST 1 / f_max(f_max(sigmap_t.x, sigmap_t.y), sigmap_t.z) / std::sqrt(20.)
#define G_MIN_DIST 1 / f_max(f_max(sigmap_t.x, sigmap_t.y), sigmap_t.z) / std::sqrt(80.)
#define MEDIUM_SCALE 2// 2//3//  1.25 // 1.5 //0.5 //6 //    (5) //2    //+++++++++++ 1.5 for 30
//+++++++1.25 for comparison with MITSUBA
// #define M_1_4PIPI 0.0253302959105844


#define INPUT_DIR "../input/"

static vec3 radiance_sss(const Ray& primary_ray, const Scene& scene, int irr_samples, bool include_sss);

static void loadModel(Scene& scene)
{
    std::cout << "Begin loading models" << std::endl;

#if 1
    //     add_CUDA(BoxObject(vec3(0.4,1.01,3.4)+world_shift,vec3(0.6,1.11,3.6)+world_shift),   
    //         vec3(1,1,1),   30,       DIFF);
    //     add_CUDA(BoxObject(vec3(0,0,3)+world_shift,vec3(1,1,4)+world_shift), vec3(1,1,1),   0,       REFR);
    //     add_CUDA(BoxObject(vec3(0,0,3)+world_shift,vec3(1,1,4)+world_shift), vec3(1,1,1),   0,       DIFF);

    //     objLoader obj("bunny_1k.obj");
    //     obj.set(4, 1,-0.28,4);
    //     objLoader obj("dragon_65k.obj");
    //     objLoader obj("..\\dragon_fix.obj");
    //     obj.set(.015, 0,0.6,4.5);

//////////////////////////////////////////////////////////////////////////

     Mesh obj(INPUT_DIR "buddha.obj");  // obj.set(1, 0,0.6,4.5);
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
//     {
//         //         Rotation rx(-90, 1, 0, 0);
//         glm::mat4 rotation = glm::rotate(glm::mat4(1.0f), glm::radians(300.0f), glm::vec3(0.0f, 1.0f, 0.0f));
//         glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
//         //         Rotation ry(300, 0, 1, 0);
//         //         Scaling sl(0.015);
//         for (int n = 0; n < obj.m_verts.size(); n++) {
//             glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
//             v = scaling * rotation * v;
//             obj.m_verts[n] = vec3(v.x, v.y, v.z);
//         }
//         obj.genNormals();
//     }

//     CMesh obj("dragon3.obj"); //obj.set(1, 0,0.6,4.5);
//     {
//         for (int n = 0; n < obj.m_verts.size(); n++)
//         {
//             glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
//             obj.m_verts[n] = vec3(v.x, v.y, v.z);
//         }
//         obj.genNormals();
//     }

    //////////////////////////////////////////////////////////////////////////

    int id = todo_addMaterial(Material(vec3(1, 1, 1), 0, SSSS));

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

        scene.add(tri);
    }
#endif

    //////////////////////////////////////////////////////////////////////////
#if 0
    {
        CMesh obj("bssrdf_floor.obj");
        {
            for (int n = 0; n < obj.m_verts.size(); n++)
            {
                glm::vec4 v = glm::vec4(obj.m_verts[n].x, obj.m_verts[n].y, obj.m_verts[n].z, 1.0f);
                obj.m_verts[n] = vec3(v.x, v.y, v.z);
            }
            obj.genNormals();
        }

        //     int id = todo_addMaterial(Material(vec3(1, 1, 1) * 0.4, 0, DIFF));
        //     int id = todo_addMaterial(Material(vec3(1, 1, 1), 0, DIFF));
        int id = todo_addMaterial(Material(vec3(1, 1, 1) * 0.5, 0, DIFF));

        for (int i = 0; i < obj.m_faces.size(); i++)
        {
            int a = obj.m_faces[i].x;
            int b = obj.m_faces[i].y;
            int c = obj.m_faces[i].z;

            TriangleObject tri(
                vec3(obj.m_verts[a].x, (obj.m_verts[a].y), obj.m_verts[a].z),
                vec3(obj.m_verts[b].x, (obj.m_verts[b].y), obj.m_verts[b].z),
                vec3(obj.m_verts[c].x, (obj.m_verts[c].y), obj.m_verts[c].z)
                );
            tri.setMaterialID(id);

            scene.add(tri);
        }
    }
#endif

#if 0
    {
        CMesh obj("dragonwall.obj");

        int id = todo_addMaterial(Material(vec3(1, 1, 1) * 0.9, 0, SPEC));

        for (int i = 0; i < obj.m_faces.size(); i++)
        {
            int a = obj.m_faces[i].x;
            int b = obj.m_faces[i].y;
            int c = obj.m_faces[i].z;

            TriangleObject tri(
                vec3(obj.m_verts[a].x, (obj.m_verts[a].y), obj.m_verts[a].z),
                vec3(obj.m_verts[b].x, (obj.m_verts[b].y), obj.m_verts[b].z),
                vec3(obj.m_verts[c].x, (obj.m_verts[c].y), obj.m_verts[c].z));
            tri.setMaterialID(id);

            scene.add(tri);
        }
    }
    {
        CMesh obj("dragonball.obj");
        {
            for (int n = 0; n < obj.m_verts.size(); n++)
            {
                obj.m_verts[n] = obj.m_verts[n] + vec3(-1, 1.5, 0);
            }
            obj.genNormals();
        }

        int id = todo_addMaterial(Material(vec3(1, 1, 1) * 0.9, 0, SPEC));

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

            scene.add(tri);
        }
    }
#endif

    std::cout << "Finished loading models" << std::endl;
}

Octree *octree;

// extern SceneBVH scene_CUDA;
// extern EnvmapLoader envmap;

//////////////////////////////////////////////////////////////////////////
float fdr(float x)
{
    return -1.44 / sq(x) + 0.71 / x + 0.668 + 0.0636*x;
}

float C1(const float n)
{
    float r;
    if (n > 1.0) {
        r = -9.23372f + n * (22.2272f + n * (-20.9292f + n * (10.2291f + n * (-2.54396f + 0.254913f * n))));
    }
    else {
        r = 0.919317f + n * (-3.4793f + n * (6.75335f + n *  (-7.80989f + n *(4.98554f - 1.36881f * n))));
    }
    return r / 2.0f;
}

float C2(const float n)
{
    float r = -1641.1f + n * (1213.67f + n * (-568.556f + n * (164.798f + n * (-27.0181f + 1.91826f * n))));
    r += (((135.926f / n) - 656.175f) / n + 1376.53f) / n;
    return r / 3.0f;
}

const float g = 0;
const vec3 sigma_a = vec3(0.013, 0.07, 0.145) * 20 * MEDIUM_SCALE;
const vec3 sigma_s = vec3(1.09, 1.59, 1.79) * 20 * MEDIUM_SCALE;
// const vec3 sigma_a = vec3(0.0011, 0.0024, 0.014) * 20 * MEDIUM_SCALE;
// const vec3 sigma_s = vec3(2.55, 3.21, 3.77)   * 20 * MEDIUM_SCALE;
const vec3 sigma_t = sigma_a + sigma_s;
// const vec3 sigma_a = vec3(0.0014, 0.0025, 0.0142) * 50;
// const vec3 sigma_s = vec3(0.70, 1.22, 1.90)       * 50;
const vec3 sigmap_s = sigma_s * (1 - g);
const vec3 sigmap_t = sigmap_s + sigma_a;
const vec3 sigma_tr = sqrt(3 * mult(sigma_a, sigmap_t));
const vec3 alphap = div(sigmap_s, sigmap_t); // reduced albedo
const vec3 D1 = div(vec3(1, 1, 1), (3 * sigmap_t));
const vec3 eta(1.49, 1.49, 1.49);
const float ETA = 1.49;
const vec3 F_dr(fdr(eta.x), fdr(eta.y), fdr(eta.z));
const vec3 F_dt = vec3(1, 1, 1) - F_dr;
const vec3 A = div((1 + F_dr), (1 - F_dr));
const vec3 z_r = div(vec3(1, 1, 1), sigmap_t);
const vec3 z_v = z_r + 4 * mult(A, D1);

const vec3 de = 2.131 * div(D1, sqrt(alphap));
const float Cp_norm = 1.0 / (1.0 - 2.0 * C1(1.0 / ETA));
const float Cp = (1.0 - 2.0 * C1(ETA)) / 4.0;
const float Ce = (1.0 - 3.0 * C2(ETA)) / 2.0;
const float AA = (1.0 - Ce) / (2.0 * Cp);

void printBSSRDFParams()
{
    std::cout << "=======================================" << std::endl;
    std::cout << "g : " << str(g) << std::endl;
    std::cout << "sigma_s : " << str(sigma_s) << std::endl;
    std::cout << "sigma_a : " << str(sigma_a) << std::endl;
    std::cout << "sigma_t : " << str(sigma_t) << std::endl;
    std::cout << "sigma_s_prime : " << str(sigmap_s) << std::endl;
    std::cout << "sigma_t_prime : " << str(sigmap_t) << std::endl;
    std::cout << "sigma_tr : " << str(sigma_tr) << std::endl;
    std::cout << "alpha_prime : " << str(alphap) << std::endl;
    std::cout << "D : " << str(D1) << std::endl;
    std::cout << "eta : " << str(eta) << std::endl;
    std::cout << "F_dr : " << str(F_dr) << std::endl;
    std::cout << "F_dt : " << str(F_dt) << std::endl;
    std::cout << "A : " << str(A) << std::endl;
    std::cout << "z_r : " << str(z_r) << std::endl;
    std::cout << "z_v : " << str(z_v) << std::endl;
    std::cout << "de : " << str(de) << std::endl;
    std::cout << "Cp_norm : " << str(Cp_norm) << std::endl;
    std::cout << "Cp : " << str(Cp) << std::endl;
    std::cout << "Ce : " << str(Ce) << std::endl;
    std::cout << "AA : " << str(AA) << std::endl;
    std::cout << "=======================================" << std::endl;
};

float getEta() {
    return ETA;
}

float fresnel(const float cos_theta_i, const float eta /* eta_t / eta_i */)
{
    const float sin_theta_t_sqr = 1.0 / (eta * eta) * (1.0 - cos_theta_i * cos_theta_i);
    if (sin_theta_t_sqr >= 1.0) return 1.0;
    const float cos_theta_t = std::sqrt(1.0 - sin_theta_t_sqr);
    const float r_s = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t); // perpendicular
    const float r_p = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t); // parallel
    return (r_s * r_s + r_p * r_p) * 0.5;
}

// real get_Fr(vec3 w_i, vec3 w_t, vec3 nl, real eta)
// {
//     real nc=1, nt=eta;
//     real nnt=nc/nt, ddn=dot(w_i,nl);
//     real cos2t=1-nnt*nnt*(1-ddn*ddn);
// 
//     real cos_theta_t = std::sqrt(cos2t);
//     real cos_theta_i = ddn;
//     real a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-cos_theta_i;
//     //Fresnel
//     real R_p = (nt*cos_theta_i - nc*cos_theta_t)/(nt*cos_theta_i + nc*cos_theta_t);
//     real R_o = (nc*cos_theta_i - nt*cos_theta_t)/(nc*cos_theta_i + nt*cos_theta_t);
//     real Re = 0.5*(R_p*R_p+R_o*R_o);
//     return Re;
// }

struct lightingSample
{
    vec3 lighting;
    vec3 position;
    vec3 normal;
    vec3 meanLightDir; // radiance-averaged dominating light direction
};
std::vector<lightingSample> perFaceLighting;

float sampleIrradiance(vec3& meanDir, const vec3& pl, const vec3& nhit, int channel, const Scene *scene)
{
    float irr = (0);
    meanDir = vec3(0);
    {
        DirLight dl = scene->envmap()->sample(randf(), randf());
        vec3 lightDirection(dl.dir[0], dl.dir[1], dl.dir[2]);
        HitInfo hit = scene->intersect(Ray(pl, lightDirection));
        if (!hit)
        {
            float g = f_max(0, dot(lightDirection, nhit)); // cosine factor for irradiance
            if (g > 0)
            {
                float weightEnvmap = scene->envmap()->getPdf(dl.dir) / (scene->envmap()->getPdf(dl.dir) + g / NUM_PI);
                irr += weightEnvmap * g / NUM_PI * dl.rad[channel] / dl.pdf; // note the 1/M_PI in the BRDF!!!!!
                meanDir += weightEnvmap * lightDirection;
            }
        }
    }
    {
        float r1 = 2 * NUM_PI*randf(), r2 = randf(), r2s = std::sqrt(r2);
        vec3 w = nhit, u = normalize(cross(vec3(.3, .4, .5), w)), v = cross(w, u);
        vec3 d = (cos(r1)*r2s*u + sin(r1)*r2s*v + std::sqrt(1 - r2)*w).normalize();
        HitInfo hit = scene->intersect(Ray(pl, d));
        if (!hit)
        {
            float g = f_max(0, dot(d, nhit)); // cosine factor for irradiance
            if (g > 0)
            {
                float weightBRDF = (g / NUM_PI) / (scene->envmap()->getPdf(&d[0]) + g / NUM_PI);
                irr += weightBRDF * scene->envmap()->sampleRadiance(&d[0], channel);
                meanDir += weightBRDF * d;
            }
        }
    }
    return irr;
}

vec3 sampleTransmittedIrradiance(vec3& meanDir, const vec3& pl, const vec3& nhit, float eta, const Scene *scene)
{
    vec3 irr(0);
    meanDir = vec3(0);

    // irradiance estimator is Li * cos_theta / pdfW

    // sample by light
    if (1)
    {
        DirLight dl = scene->envmap()->sample(randf(), randf());
        vec3 lightDirection(dl.dir[0], dl.dir[1], dl.dir[2]);
        //         Hit_CUDA hit = scene->getSceneBVH()->intersect(ray3(pl, lightDirection));
        //         bool occluded = hit.t != __INFINITY;
        bool occluded = scene->occluded(Ray(pl, lightDirection));
        if (!occluded) {
            float g = f_max(0, dot(lightDirection, nhit)); // cosine factor for irradiance
            if (g > 0) {
                float weightEnvmap = scene->envmap()->getPdf(dl.dir) / (scene->envmap()->getPdf(dl.dir) + g / NUM_PI);
                irr += (weightEnvmap * g / dl.pdf) * vec3(dl.rad[0], dl.rad[1], dl.rad[2]); // note the 1/M_PI in the BRDF!!!!! which gets cancelled
                meanDir += weightEnvmap * lightDirection;
            }
        }
    }

    // sample by lambertian BRDF
    if (1)
    {
        float r1 = 2 * NUM_PI*randf(), r2 = randf(), r2s = std::sqrt(r2);
        vec3 w = nhit, u = normalize(cross(vec3(.3, .4, .5), w)), v = cross(w, u);
        vec3 d = (cos(r1)*r2s*u + sin(r1)*r2s*v + std::sqrt(1 - r2)*w).normalize();
        //         Hit_CUDA hit = scene->getSceneBVH()->intersect(ray3(pl, d));
        //         bool occluded = hit.t != __INFINITY;
        bool occluded = scene->occluded(Ray(pl, d));
        if (!occluded) {
            float g = f_max(0, dot(d, nhit)); // cosine factor for irradiance
            if (g > 0) {
                float weightBRDF = (g / NUM_PI) / (scene->envmap()->getPdf(&d[0]) + g / NUM_PI);
                vec3 rad;
                scene->envmap()->sampleRadiance(&d[0], &rad[0]);
                irr += (weightBRDF * NUM_PI) * rad; // pi is from irradiance estimator; cos_wo cancelled with pdf
                meanDir += weightBRDF * d;
            }
        }
    }

    if (meanDir.x != meanDir.x || meanDir.y != meanDir.y || meanDir.z != meanDir.z) {
        meanDir = vec3(0);
    }

    if (irr.x != irr.x || irr.y != irr.y || irr.z != irr.z) {
        irr = vec3(0);
    }

    return irr;
}

vec3 sampleIrradianceFinal(vec3& meanDir, const vec3& pl, const vec3& nhit, float eta, const Scene *scene)
{
    vec3 irr(0);
    meanDir = vec3(0);

    // irradiance estimator is Li * cos_theta / pdfW
    // sample by weighted cosine
    float pdfW;
    vec3 s = sampleHemisphereCosineW(randf(), randf(), &pdfW);
    OrthonormalFrame frame(nhit);
    vec3 d = frame.toWorld(s);

    float cos_theta = f_max(0, dot(d, nhit)); // cosine factor for irradiance
    if (cos_theta > 0)
    {
        vec3 rad = radiance_sss(Ray(pl, d), *scene, 1, false);
        irr += (cos_theta / pdfW) * rad;
        meanDir = d;
    }
    meanDir = normalize(meanDir);

    if (meanDir.x != meanDir.x ||
        meanDir.y != meanDir.y ||
        meanDir.z != meanDir.z)
    {
        meanDir = vec3(0);
    }

    if (irr.x != irr.x || irr.y != irr.y || irr.z != irr.z) {
        irr = vec3(0);
    }

    return irr;
}

static void initScene(Scene &scene)
{
    loadModel(scene);

    //     add_CUDA(vec3(1,1,5), 0.5, DIFF);
    //     add_CUDA(vec3(1.3,1.4,5), 0.3, DIFF);

    scene.finalize();

//     center = scene->getSceneBVH()->get_center();

    //////////////////////////////////////////////////////////////////////////
    //     real minDist = 0.2 * 0.5 * 0.5 * 0.5 *  1 / f_max(f_max(reduced_sigma_t.x, reduced_sigma_t.y), reduced_sigma_t.z); //for 40000 samples
    //     real minDist = 16 * 2 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 *  1 / f_max(f_max(reduced_sigma_t.x, reduced_sigma_t.y), reduced_sigma_t.z); //for 40000 samples
    //     real minDist = 10 * 2 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 * 0.5 *  1 / f_max(f_max(reduced_sigma_t.x, reduced_sigma_t.y), reduced_sigma_t.z); //for 40000 samples
    float minDist = G_MIN_DIST;
    int gSkipCount = 0;

    printf("F_dr: %f, %f, %f\n", F_dr[0], F_dr[1], F_dr[2]);
    printf("F_dt: %f, %f, %f\n", F_dt[0], F_dt[1], F_dt[2]);

    printf("sigma_a: %f, %f, %f\n", sigma_a[0], sigma_a[1], sigma_a[2]);
    printf("sigma_s: %f, %f, %f\n", sigma_s[0], sigma_s[1], sigma_s[2]);
    printf("sigma_t: %f, %f, %f\n", sigma_t[0], sigma_t[1], sigma_t[2]);
    printf("sample radius: %f\n", minDist);

    perFaceLighting.clear();

    float totalArea = 0;
    float lightingFactor;

    if (REBUILD) //rebuild samples
    {
        printf("total face count: %d\n", scene.numTrangles());

        std::vector<int> bssrdf_to_global_idx;
        std::vector<float> triAreas;

        for (int i = 0; i < scene.numTrangles(); i++)
        {
            if (todo_getMaterial(scene.triangle(i).materialID()).type() == SSSS)
            {
                triAreas.push_back(scene.triangleArea(i));
                bssrdf_to_global_idx.push_back(i);
                totalArea += triAreas[i];
            }
        }
        printf("bssrdf mesh surface area: %f\n", totalArea);
        printf("bssrdf face count: %d\n", bssrdf_to_global_idx.size());

        DiscreteSampler dsample(triAreas);

        int sampleTotalCount = G_TOTAL_SAMPLE;
        std::vector<Position> p_in;
        while (p_in.size() < sampleTotalCount)
        {
            int n = p_in.size();
            if (n % 1000 == 0)
            {
                printf("\r" "presampling %f %%", n / (sampleTotalCount - 1.0f) * 100.0f);
            }

            int currentIdx = dsample.sample();
            currentIdx = bssrdf_to_global_idx[currentIdx];

            auto& currentTri = scene.triangle(currentIdx);

            float u, v;
            vec3 sampleLocation = currentTri.samplePoint(u, v);
            sampleLocation += currentTri.normal(u, v) * 1e-4; // a little offset

            Position temp(sampleLocation.x, sampleLocation.y, sampleLocation.z);
            temp.payload.triangle_idx = currentIdx;
            temp.payload.u = u;
            temp.payload.v = v;
            p_in.push_back(temp);
        }

        std::vector<Position> p_out;
        poissonSample(p_in, p_out, minDist);

        printf("Sampling efficiency: %f%%\n", 100.0f * p_out.size() / p_in.size());

        if (0) check_min_distance(p_out); // takes too much time for large sample set

        for (int n = 0; n < p_out.size(); n++) {
            TriangleObject &sampleTriangle = scene.triangle(p_out[n].payload.triangle_idx);

            lightingSample sample;
            sample.position = vec3(p_out[n].x, p_out[n].y, p_out[n].z);
            sample.normal = sampleTriangle.normal(p_out[n].payload.u, p_out[n].payload.v);
            //             sample.normal = sampleTriangle.getFaceNormal();

            if (sample.normal.lengthSquared() < 0.001*0.001) {
                printf("\t ----- degenerate normal: %f, %f, %f\n", sample.normal.x, sample.normal.y, sample.normal.z);
            }
            sample.normal = sample.normal.normalize();

            perFaceLighting.push_back(sample);
        }

        {
            FILE *fp = fopen("sample_cache.bin", "wb");
            int totalBytes = perFaceLighting.size() * sizeof(lightingSample);
            fwrite(&totalBytes, 1, sizeof(int), fp);
            fwrite(&perFaceLighting[0], 1, totalBytes, fp);
            fclose(fp);

            lightingFactor = totalArea / perFaceLighting.size();
            printf("lightingFactor: %f = Area(%f) / #sample(%d)\n", lightingFactor, totalArea, perFaceLighting.size());
        }
    }
    else
    {
        FILE *fp = fopen("sample_cache.bin", "rb");
        int totalBytes;
        fread(&totalBytes, 1, sizeof(int), fp);

        int nSamps = totalBytes / sizeof(lightingSample);
        perFaceLighting.resize(nSamps);

        fread(&perFaceLighting[0], 1, totalBytes, fp);
        fclose(fp);
    }

    printf("\n");
    printf("# samples on mesh: %d\n", perFaceLighting.size());
    printf("mean free path: %f, %f, %f\n", 1 / sigmap_t.x, 1 / sigmap_t.y, 1 / sigmap_t.z);
    printf("diffuse mfp: %f, %f, %f\n", 1 / sigma_tr.x, 1 / sigma_tr.y, 1 / sigma_tr.z);
    printf("skip count: %d\n", gSkipCount);

    printf("\n");

    if (REBAKE)
    {
        for (int n = 0; n < perFaceLighting.size(); n++)
        {
            printf("\rlight baking %f%%", float(n) / (perFaceLighting.size() - 1) * 100.0f);
            lightingSample& sample = perFaceLighting[n];

            sample.lighting = vec3(0, 0, 0);
            sample.meanLightDir = vec3(0, 0, 0);

            int tnum = omp_get_max_threads();
            std::vector<lightingSample> ompSamples(tnum, sample);

            int nn = G_IRRADIANCE_SAMPLE;
            float oneOverNumIrrSamples = 1.0 / nn;

            vec3 meanLightDir(0, 0, 0);

#if 0 // old method
#pragma omp parallel for
            for (int i = 0; i < nn; i++)
            {
                int tid = omp_get_thread_num();
                DirLight dl = scene->envmap()->generateLight();
                vec3 lightDirection(dl.dir[0], dl.dir[1], dl.dir[2]);
                Hit_CUDA hit = scene->getSceneBVH()->intersect(ray3(sample.position, lightDirection));
                if (hit.t == __INFINITY) {
                    float g = f_max(0, dot(lightDirection, sample.normal)); // cosine factor for irradiance
                    real f_t = 1 - fresnel(g, ETA); // Fresnel reflectance

                    ompSamples[tid].lighting[0] += g * dl.rad[0] / dl.pdf      * f_t;
                    ompSamples[tid].lighting[1] += g * dl.rad[1] / dl.pdf      * f_t;
                    ompSamples[tid].lighting[2] += g * dl.rad[2] / dl.pdf      * f_t;

                    float w = EnvmapLoader::luminance(dl.rad[0], dl.rad[1], dl.rad[2]);
                    meanLightDir = meanLightDir + lightDirection * w;
                }
            }
#else // MIS
//             struct DebugRecord
//             {
//                 vec3 dir;
//                 float w;
//             };
//             std::vector<DebugRecord> rec(tnum);

#pragma omp parallel for
            for (int i = 0; i < nn; i++) {
                int tid = omp_get_thread_num();
                vec3 meanLightDirThread(0);
//                 vec3 rad = sampleIrradianceFinal(meanLightDirThread, sample.position, sample.normal, ETA, &scene);
                vec3 rad = sampleTransmittedIrradiance(meanLightDirThread, sample.position, sample.normal, ETA, &scene);
                ompSamples[tid].lighting += rad;
                float w = EnvmapLoader::luminance(rad[0], rad[1], rad[2]);
                meanLightDir += meanLightDirThread * w;

//                 rec[tid].dir = meanLightDirThread;
//                 rec[tid].w = w;
            }

//             vec3 debug_meanLightDir(0, 0, 0);
//             for (auto& r : rec)
//             {
//                 debug_meanLightDir += r.dir * r.w;
//             }
//             if (length(debug_meanLightDir) == 0)
//             {
//                 cout << str(debug_meanLightDir) << std::endl;
//             }
#endif

            for (int n = 0; n < tnum; n++) {
                if (ompSamples[n].lighting[0] != ompSamples[n].lighting[0] ||
                    ompSamples[n].lighting[1] != ompSamples[n].lighting[1] ||
                    ompSamples[n].lighting[2] != ompSamples[n].lighting[2])
                {
                    printf("radiance sample discarded!\n");
                    continue;
                }

                sample.lighting += ompSamples[n].lighting;
            }

            // the inv_pi from the BSSRDF and the pi from irradiance are canceled out
            // FIX: the pi from irradiance is not universal, unless a pdf like cos(theta)/pi or 1/2pi is used,
            // thus, the inv_pi from the bssrdf must be explicitly considered !!!!!!!!
            sample.lighting[0] *= oneOverNumIrrSamples; // * F_dt[0]
            sample.lighting[1] *= oneOverNumIrrSamples; // * F_dt[1]
            sample.lighting[2] *= oneOverNumIrrSamples; // * F_dt[2]
                    // F_dt moved to bssrdf evaluation
                    // in order not to disrupt the light cache

            float l = length(meanLightDir);
            if (l != 0)
            {
                sample.meanLightDir = meanLightDir / l;
            }
//             sample.meanLightDir = meanLightDir.norm();

            {
                //                 printf("li: %f, %f, %f\n", sample.lighting.x, sample.lighting.y, sample.lighting.z);
                //                 printf("nl: %f, %f, %f\n", sample.normal.x, sample.normal.y, sample.normal.z);
                //                 printf("di: %f, %f, %f\n", sample.meanLightDir.x, sample.meanLightDir.y, sample.meanLightDir.z);
            }
        }

        {
            std::cout << std::endl;
            std::cout << ">> before fix : " << perFaceLighting.size() << std::endl;
            std::vector<lightingSample> temp;
            for (auto& r : perFaceLighting)
            {
                if (length(r.meanLightDir) > 0)
                {
                    temp.push_back(r);
                }
            }
            perFaceLighting.swap(temp);
            std::cout << "after fix : " << perFaceLighting.size() << std::endl;
        }

        {
            FILE *fp = fopen("light_cache.bin", "wb");
            int totalBytes = perFaceLighting.size() * sizeof(lightingSample);
            fwrite(&totalBytes, 1, sizeof(int), fp);
            fwrite(&perFaceLighting[0], 1, totalBytes, fp);
            fclose(fp);

            lightingFactor = totalArea / perFaceLighting.size();
            printf("lightingFactor: %f = Area(%f) / #sample(%d)\n", lightingFactor, totalArea, perFaceLighting.size());

            // parameters baked with only the bssrdf object
            fp = fopen("parameter_cache.bin", "wb");
            fwrite(&totalArea, 1, sizeof(float), fp);
            fwrite(&lightingFactor, 1, sizeof(float), fp);
            fclose(fp);
        }
    }
    else
    {
        FILE *fp = fopen("light_cache.bin", "rb");
        int totalBytes;
        fread(&totalBytes, 1, sizeof(int), fp);

        int nSamps = totalBytes / sizeof(lightingSample);
        perFaceLighting.resize(nSamps);

        fread(&perFaceLighting[0], 1, totalBytes, fp);
        fclose(fp);

//         lightingFactor = totalArea / perFaceLighting.size();
//         printf("lightingFactor: %f\n", lightingFactor);

        fp = fopen("parameter_cache.bin", "rb");
        fread(&totalArea, 1, sizeof(float), fp);
        fread(&lightingFactor, 1, sizeof(float), fp);
        fclose(fp);
        printf("lightingFactor: %f\n", lightingFactor);
    }

    printf("\n");
    printf("# samples on mesh: %d\t\t\t\t\n", perFaceLighting.size());
    printf("mean free path: %f, %f, %f\n", 1 / sigmap_t.x, 1 / sigmap_t.y, 1 / sigmap_t.z);
    printf("diffuse mfp: %f, %f, %f\n", 1 / sigma_tr.x, 1 / sigma_tr.y, 1 / sigma_tr.z);
    printf("skip count: %d\n", gSkipCount);

    {
        vec3 posMin(1e10, 1e10, 1e10);
        vec3 posMax(-1e10, -1e10, -1e10);
        for (auto it = perFaceLighting.begin(); it != perFaceLighting.end(); ++it) {
            lightingSample& rad = *it;
            posMin = f_min(posMin, rad.position);
            posMax = f_max(posMax, rad.position);
        }
        posMin -= vec3(.1, .1, .1);  // some margin
        posMax += vec3(.1, .1, .1);  // some margin

        printf("Octree bounding box:\n");
        printf("\tmin: %f, %f, %f\n", posMin.x, posMin.y, posMin.z);
        printf("\tmax: %f, %f, %f\n", posMax.x, posMax.y, posMax.z);

        octree = new Octree(
            (posMin + posMax)*.5f,
            (posMax - posMin)*.5f);
        for (auto it = perFaceLighting.begin(); it != perFaceLighting.end(); ++it){
            lightingSample& rad = *it;

            if (rad.position != rad.position)
            {
                std::cout << "1!!!!!!!!!!!!!" << str(rad.position) << std::endl;
            }
            if (rad.normal != rad.normal)
            {
                std::cout << "2!!!!!!!!!!!!!" << str(rad.normal) << std::endl;
            }
            if (rad.meanLightDir != rad.meanLightDir)
            {
                std::cout << "3!!!!!!!!!!!!!" << str(rad.meanLightDir) << std::endl;
            }
            if (rad.lighting != rad.lighting)
            {
                std::cout << "4!!!!!!!!!!!!!" << str(rad.lighting) << std::endl;
            }

//             rad.lighting *= F_dt;
            octree->insert(new OctreeSample(rad.position, rad.normal, rad.meanLightDir, rad.lighting, lightingFactor));
            //             octree->insert(new OctreeSample(rad.position, rad.normal, rad.meanLightDir, mult(rad.lighting, F_dt), lightingFactor));
        }
        octree->finalize();
    }

    printf("init complete\n");
}

float Rd(float r, int channel)
{
    float d_r = std::sqrt(sq(r) + sq(z_r[channel]));
    float d_v = std::sqrt(sq(r) + sq(z_v[channel]));
    float d_r_3 = d_r * d_r * d_r;
    float d_v_3 = d_v * d_v * d_v;

    float c1 = z_r[channel] * (sigma_tr[channel] * d_r + 1);
    float c2 = z_v[channel] * (sigma_tr[channel] * d_v + 1);

    float R = alphap[channel] / (4 * NUM_PI) * (

        c1 * expf(-sigma_tr[channel] * d_r) / (d_r_3)+
        c2 * expf(-sigma_tr[channel] * d_v) / (d_v_3));

    return R;
}

#if DIRPOLE
float Sp_d(const vec3 x, const vec3 w, const float r /* radius */, const vec3 n /* normal */, const int j /* channel */)
{
    // evaluate the profile
    const float s_tr_r = sigma_tr[j] * r;
    const float s_tr_r_one = 1.0 + s_tr_r;
    const float x_dot_w = dot(x, w);
    const float r_sqr = r * r;

    const float t0 = Cp_norm * M_1_4PIPI * expf(-s_tr_r) / (r * r_sqr);
    const float t1 = r_sqr / D1[j] + 3.0 * s_tr_r_one * x_dot_w;
    const float t2 = 3.0 * D1[j] * s_tr_r_one * dot(w, n);
    const float t3 = (s_tr_r_one + 3.0 * D1[j] * (3.0 * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * dot(x, n);

    return t0 * (Cp * t1 - Ce * (t2 - t3));
}

float bssrdf(const vec3 xi, const vec3 ni, const vec3 wi, /* light */ const vec3 xo, const vec3 no, const vec3 wo, /* eye */ const int j /* channel */)
{
    // distance
    const vec3 xoxi = xo - xi;
    float r = xoxi.length();

    //     r = f_max(r, Rmin);
    //     r = std::sqrt(sq(r) + sq(Rmin));

    // modified normal
    //     const vec3 ni_s = cross(xoxi.norm(), cross(ni, xoxi).norm());
    float r_sqr = dot(xoxi, xoxi);
    const vec3 ni_s = r_sqr < 1.0e-12 ? ni : (ni*r_sqr - xoxi*dot(xoxi, ni)).normalize();

    // directions of ray sources
    const float nnt = 1.0 / ETA, ddn = -dot(wi, ni);
    const vec3 wr = (wi * -nnt - ni * (ddn * nnt + std::sqrt(1.0 - nnt * nnt * (1.0 - ddn * ddn)))).normalize();  // refraction
    //     if (0 == wr.length())
    //     {
    //         printf("wr: %f, %f, %f,\n", wr.x, wr.y, wr.z);
    //         printf("wi: %f, %f, %f,\n", wi.x, wi.y, wi.z);
    //         printf("ni: %f, %f, %f,\n", ni.x, ni.y, ni.z);
    //         printf("nnt: %f,\n", nnt);
    //         printf("ddn: %f,\n", ddn);
    //     }
    const vec3 wv = wr - ni_s * (2.0 * dot(wr, ni_s)); // reflection

    // distance to real sources
    const float cos_beta = -std::sqrt((r * r - dot(xoxi, wr) * dot(xoxi, wr)) / (r * r + de[j] * de[j]));
    float dr;
    const float mu0 = -dot(no, wr);
    if (mu0 > 0.0) {
        dr = std::sqrt((D1[j] * mu0) * ((D1[j] * mu0) - de[j] * cos_beta * 2.0) + r * r);
    }
    else {
        dr = std::sqrt(1.0 / (3.0 * sigma_t[j] * 3.0 * sigma_t[j]) + r * r);
    }

    // distance to virtual source
    const vec3 xoxv = xo - (xi + ni_s * (2.0 * AA * de[j]));
    float dv = xoxv.length();

    //     dr = f_max(dr, Rmin);
    //     dv = f_max(dv, Rmin);

    //     printf("debug: %f, %f, %f,\n", xoxi.x, xoxi.y, xoxi.z);
    //     printf("debug: %f, %f, %f,\n", wr.x, wr.y, wr.z);
    //     printf("debug: %f,\n", dr);
    //     printf("debug: %f, %f, %f,\n", no.x, no.y, no.z);

    // BSSRDF
    const float result = Sp_d(xoxi, wr, dr, no, j) - Sp_d(xoxv, wv, dv, no, j);

    // clamping to zero
    return (result < 0.0) ? 0.0 : result;
}

#else

float Sp_d(const float zr, const float zv, const float r /* dist */, const int j /* channel */)
{
    // evaluate the profile
    const float dr = std::sqrt(zr*zr + r*r);
    const float dv = std::sqrt(zv*zv + r*r);

    const float r_s_tr_r = sigma_tr[j] * dr;
    const float r_s_tr_r_one = 1.0 + r_s_tr_r;
    const float v_s_tr_r = sigma_tr[j] * dv;
    const float v_s_tr_r_one = 1.0 + v_s_tr_r;

    return alphap[j] * M_1_4PIPI * (
        zr * r_s_tr_r_one * expf(-r_s_tr_r) / (dr * dr * dr) +
        zv * v_s_tr_r_one * expf(-v_s_tr_r) / (dv * dv * dv)
        );
}

float bssrdf( const vec3 xi, const vec3 ni, const vec3 wi, /* light */ const vec3 xo, const vec3 no, const vec3 wo, /* eye */ const int j /* channel */ )
{
    // distance
    const float r = (xo - xi).length();

    // distance to real sources
    float zr = 3 * D1[j];
    float zv = zr + 4 * A[j] * D1[j];

    // BSSRDF
    const float result = Sp_d(zr, zv, r, j);

    // clamping to zero
    //     return ((result < 0.0) ? 0.0 : result) / M_PI; // inv_pi is moved out of the loop
    return ((result < 0.0) ? 0.0 : result);
}
#endif

#if DIRPOLE
DirectionalDipoleBSSRDF sss(sigma_s, sigma_a, g, eta);
#else
StandardDipoleBSSRDF sss(sigma_s, sigma_a, g, eta);
#endif

vec3 Bssrdf(const vec3& xi, const vec3& ni, const vec3& wi, /* light */ const vec3& xo, const vec3& no, const vec3& wo) // eye
{
#if 0
    // old code
//     return vec3(
//         bssrdf(xi, ni, wi, xo, no, wo, 0),
//         bssrdf(xi, ni, wi, xo, no, wo, 1),
//         bssrdf(xi, ni, wi, xo, no, wo, 2)) * F_dt * M_PI; // there's no where to put this pi, but it helps to modulate the result nicely when F_dt is used; can replace with Ft and pi is not needed
#else
    // new code from shader
    return sss.evaluate(xi, ni, wi, xo, no, wo);
#endif
}

vec3 Constant(const vec3& xi, const vec3& ni, const vec3& wi, /* light */ const vec3& xo, const vec3& no, const vec3& wo) // eye
{
    return vec3(1);
}
//////////////////////////////////////////////////////////////////////////

static float MIS_weight(float pdfA, float pdfB) {
    return (pdfA) / (pdfA + pdfB);
}

class BSDF
{
public:
    float getPDF(const vec3& nl, const vec3& dir) const {
        return f_max(0, dot(nl, dir)) / NUM_PI;
    }
    float eval(const vec3& wi, const vec3& wo, const vec3& nl) const {
        return 1 / NUM_PI;
    }
};

class EnvmapSampler {
    const Scene *pscene;

public:
    EnvmapSampler(const Scene* scene)
    {
        pscene = scene;
    }
    vec3 sampleEnvmapUsingLightPDF(const vec3& wi, const vec3& pos, const vec3& nl, float& pdf, vec3& dir,
        const BSDF& bsdf) {
        vec3 irr(0);
        pdf = 1; // using pdf = 0 causes problem due to NANs
        dir = vec3(0, 0, 0);
        DirLight dl = pscene->envmap()->sample(randf(), randf());
        vec3 lightDirection(dl.dir[0], dl.dir[1], dl.dir[2]);
        bool occluded = pscene->occluded(Ray(pos, lightDirection));
        if (!occluded) {
            float g = f_max(0, dot(lightDirection, nl)); // cosine factor for irradiance
            if (g > 0) {
                //float weightEnvmap = pscene->envmap()->getPdf(dl.dir) / (pscene->envmap()->getPdf(dl.dir) + g / M_PI);
                pdf = pscene->envmap()->getPdf(dl.dir);
                dir = lightDirection;
                irr += (g * bsdf.eval(wi, lightDirection, nl) / dl.pdf) * vec3(dl.rad[0], dl.rad[1], dl.rad[2]); // note the 1/M_PI in the BRDF!!!!!
            }
        }
        if (irr != irr) {
            irr = vec3(0);
        }
        return irr;
    }
    vec3 sampleEnvmapUsingSurfacePDF(const vec3& wi, const vec3& pos, const vec3& nl, float& pdf, vec3& dir,
        const BSDF& bsdf) {
        // NOW JUST INLINES ideal diffuse BSDF
        vec3 irr(0);
        pdf = 1; // using pdf = 0 causes problem due to NANs
        dir = vec3(0, 0, 0);
        float r1 = 2 * NUM_PI * randf(), r2 = randf(), r2s = std::sqrt(r2);
        vec3 w = nl, u = normalize(cross(vec3(.3, .4, .5), w)), v = cross(w, u);
        vec3 d = (cos(r1)*r2s*u + sin(r1)*r2s*v + std::sqrt(1 - r2)*w).normalize();
        bool occluded = pscene->occluded(Ray(pos, d));
        if (!occluded) {
            float g = f_max(0, dot(d, nl)); // cosine factor for irradiance
            if (g > 0) {
                //float weightBRDF = (g / M_PI) / (pscene->envmap()->getPdf(&d[0]) + g / M_PI);
                // pdf = g / M_PI;
                pdf = bsdf.getPDF(nl, d);
                dir = d;
                vec3 rad;
                pscene->envmap()->sampleRadiance(&d[0], &rad[0]);
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

static vec3 radiance_sss(const Ray& primary_ray, const Scene& scene, int irr_samples = 64, bool include_sss = true) //Monte Carlo path tracing
{
    HitInfo hit;
    Ray cr(primary_ray);//current ray
    vec3 radianceIntegral(0, 0, 0);
    vec3 colorState(1, 1, 1);
    const int maxDepth = 4; // 3;//2;
    //     _CONST_ int maxDepth = maxDepths[channel];
    //     _CONST_ int maxDepth = default_medium.maxDepths[channel];

    //     real iors[3] = {1.5, 1.5, 1.5};
    //     real iors[3] = {1e5, 1e5, 1e5};
    const VolumetricProps &volProps = default_medium;

    MaterialType last_mat;

    for (int n = 0; n < maxDepth; n++)
    {

        hit = scene.intersect(cr);
        if (!hit)
        { //background color

#if HIDE_BACKGROUND
            if (n > 0 && last_mat != DIFF)
#else
            if (n == 0 || last_mat != DIFF)
#endif
            {
                vec3 r;
                scene.envmap()->sampleRadiance(&cr.dir[0], &r[0]);
                radianceIntegral += colorState * r;
            }

                //             radianceIntegral += colorState * ambient_intensity(cr.d)[channel]; 
            //             radianceIntegral += colorState * ambient_intensity(cr.d)[channel] * (cr.d.y>0.2); 
            break;
        }

        vec3 phit = cr.at(hit.distance());//calculate the hit point, the last term is the surface fudge factor
        vec3 nhit(hit.shadingNormal());
        bool inside = false;
        vec3 nhit0(nhit);
        const TriangleObject *obj = hit.triangleObject();
        vec3 nhit_g = obj->normal();

        // in case of flipped triangles
        if (dot(cr.dir, nhit) > 0)
        {
            nhit = -nhit;
            nhit_g = -nhit_g;
            inside = true;
        }

//         return vec3(fabs(dot(cr.d, nhit)));

//         if (dot(cr.d, nhit) > 0)
//         {
//             break;
//         }

//         nhit = nhit_g;

        // normal fix (not useful)
//         if (dot(nhit, nhit_g) <= 0) nhit = nhit_g;

        vec3 emissionTotal(0, 0, 0);
        const Material& mat = todo_getMaterial(obj->materialID());
        vec3 objColor = mat.color();
        const vec3 objEmission = vec3(1, 1, 1) * mat.emission();

        switch (mat.type())
        {
        default:
        case DIFF:
        {
            vec3 weight;
            vec3 d = mat.sample(-cr.dir, nhit, randf(), randf(), weight);

            vec3 rad(0);
            for (int i = 0; i < irr_samples; i++)
            {
                EnvmapSampler rs(&scene);
                float envmapPDF, bsdfPDF;
                vec3 envmapDir, bsdfDir;
                vec3 rad_envmap = rs.sampleEnvmapUsingLightPDF(cr.dir, phit, nhit, envmapPDF, envmapDir, BSDF());
                vec3 rad_bsdf = rs.sampleEnvmapUsingSurfacePDF(cr.dir, phit, nhit, bsdfPDF, bsdfDir, BSDF());

                float weightEnvmap = MIS_weight(envmapPDF, BSDF().getPDF(nhit, envmapDir));
                float weightBSDF = MIS_weight(bsdfPDF, scene.envmap()->getPdf(&bsdfDir[0]));

                rad += rad_envmap * weightEnvmap + rad_bsdf * weightBSDF;
            }
            radianceIntegral += colorState * mat.color() * rad / irr_samples;

            colorState *= weight;
            cr = Ray(phit, d);
        }
        break;

        case SPEC:
        {
            vec3 weight;
            vec3 d = mat.sample(-cr.dir, nhit, randf(), randf(), weight);
            float f_r = 1;// fresnel(f_max(-dot(cr.d, nhit), 0.0), getEta()); // Fresnel reflectance
            colorState *= weight * f_r;
            cr = Ray(phit, d);
        }
        break;

#if 1
        case SSSS:
        {
            //             real r1=2*M_PI*randf(), r2=randf(), r2s=std::sqrt(r2);//r2s = sin(phi), where phi is reflection angle, cos(phi) as tangent component, sin(phi) as normal component of the reflected ray
            //             vec3 w=nhit, u=(vec3(.3,.4,.5)^w).norm(), v=w^u;//w,u,v are all unit vectors, constructs orthonormal frame
            //             vec3 d = (cos(r1)*r2s*u + sin(r1)*r2s*v + std::sqrt(1-r2)*w).norm();//uniformly distributed sampling ray on the reflection hemisphere with uniformly sampled reflection angle
//             emissionTotal += objEmission;
            //             colorState = mult(colorState, objColor); // TODO: check if use this with dipole !!!!!!!

            float f_r = fresnel(f_max(-dot(cr.dir, nhit), 0.0), getEta()); // Fresnel reflectance
            f_r = clampf(f_r, 0, 1);

            // diffuse shading
            //             radianceIntegral += mult(colorState, vec3(f_max(0, -dot(cr.d, nhit)))); //+++++++++++

            float threshold = n == 0 ? 0.2 : 0.4;
//             threshold = sq(0.25);

//             if (n == 0)
//             if (n == 0 || last_mat == BSSRDF)
            {
                //             real f_rr = fresnel(-dot(cr.d,nhit_g), ETA); // Fresnel reflectance
                //             f_rr = clampf(f_rr, 0, 1);
                //             const real T21 = 1.0 - f_min(0.1, fresnel(f_max(0, dot(-cr.d, nhit)), ETA));
                //                 radianceIntegral += octree->fastsumRadiance(pl, nhit, -cr.d)  * (1 - f_r); // correct but strange looking? (may need deeper paths, because the primary intersection may have lower transmittance)
                //             radianceIntegral = octree->fastsumRadiance(pl, nhit, -cr.d)  * (1 - f_rr);
                //             radianceIntegral = octree->fastsumRadiance(pl, nhit, -cr.d);
                //             radianceIntegral = vec3(f_max(0, dot(-cr.d, nhit)));
                //             radianceIntegral = (1 - f_r);
                //             radianceIntegral = vec3(f_max(0, dot(-cr.d, nhit_g)));
                //             radianceIntegral = vec3(f_max(0, dot(-cr.d, nhit))) * (1 - f_r);

                // -----------------------------------------------------------
                // for old bssrdf code !!!!!!!!!!!!!!!!!!!
                // there's no where to put this pi, but it helps to modulate the result nicely when F_dt is used
                //vec3 octree_rad = octree->fastsumRadiance(phit, nhit_g, -cr.d, theta)  * (1 - f_r) * NUM_PI; // todo : figure out where to put this pi
                
                // for new shader bssrdf code !!!!!!!!!!!!!!!!!!!
                if (include_sss)
                {
                    vec3 octree_rad = octree->fastsumRadiance(phit, nhit_g, -cr.dir, threshold, Bssrdf)  * (1 - f_r);
                    radianceIntegral += colorState * octree_rad;
                }
            }

            vec3 reflDir((cr.dir - 2 * dot(cr.dir, nhit) * nhit).normalize());     // Ideal dielectric REFRACTION

    // debug
//             nhit = hit.getConsistentNormal(-cr.d, &reflDir); // consistent reflection
//             return radianceIntegral += colorState * f_max(-dot(cr.d, nhit), 0.0);

            //             if (0) // fresnel reflection // disable this code block and use passive evaluation instead
            //             {
            // #if 0
            //                 // the trick is to always use lighting for fresnel reflection without occlusion,
            //                 // such that the dark edges of dipole evaluation are covered ( is it used in mitsuba)
            //                 // update: clearly this trick is not used by mitsuba, so deep tracing (>= 2) is a must,
            //                 // which means at least the dipole contribution from the primary hit and the next (modulated by f_r) must be considered
            //                 radianceIntegral += vec3(
            //                     envmap.sampleRadiance(&reflDir[0], 0),
            //                     envmap.sampleRadiance(&reflDir[0], 1),
            //                     envmap.sampleRadiance(&reflDir[0], 2)) * f_r;
            // #else
            //                 Hit_CUDA hit = scene->getSceneBVH()->intersect(ray3(pl, reflDir));
            //                 if (hit.t == infinity) {
            // //                     real nc=1, nt=ETA;//water
            // //                     real nnt=nc/nt, ddn=dot(cr.d,nhit);
            // //                     real cos2t=1-nnt*nnt*(1-ddn*ddn);
            // //                     vec3 tdir = (cr.d*nnt - nhit*(ddn*nnt+std::sqrt(cos2t))).norm();
            // //                     real a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1 + ddn;//R0 is the reflectance at normal incidence
            // //                     real Re=R0+(1-R0)*c*c*c*c*c;//Schlick's approximation for Fresnel reflectance
            // 
            //                     radianceIntegral += vec3(
            //                         envmap.sampleRadiance(&reflDir[0], 0),
            //                         envmap.sampleRadiance(&reflDir[0], 1),
            //                         envmap.sampleRadiance(&reflDir[0], 2)) * f_r;
            //                 }
            // #endif
            //             }

            // clearly at least depth = 2 is required for the dipole rendering to get rid of edge artifacts,
            // the dipole contribution from primary hit is modulated by fresnel transmittance, which is vanishing
            // at the gazing angle, but the reflection is high, thus it requires tracing deeper to evaluate
            // dipole contribution that manifest through this fresnel reflection. higher order dipole contribution
            // from even deeper trace is not as important though
            colorState *= f_r;
            cr = Ray(phit, reflDir);
        }break;
#endif

        }

        last_mat = mat.type();

//         radianceIntegral += mult(colorState, emissionTotal);
    }

    if (radianceIntegral != radianceIntegral)
    {
        radianceIntegral = vec3(1, 0, 0);
    }

    return radianceIntegral;
}

static float renderEdge(const Ray& primary_ray, const Scene& scene)
{
    Ray cr(primary_ray);
    float ret = 0;

    for (int i = 0; i < 10; i++)
    {
        HitInfo hit = scene.intersect(cr);
        if (!hit)
        {
            break;
        }

        vec3 phit = cr.at(hit.distance());
        vec3 nhit(hit.shadingNormal());
        const TriangleObject *obj = hit.triangleObject();
        vec3 nhit_g = obj->normal();
        auto& mat = todo_getMaterial(obj->materialID());

        if (!mat.isDelta())
        {
            break;
        }

        if (/*mat.type() == SSSS &&*/ -dot(cr.dir, nhit) < 0.1)
        {
            ret += 1;
        }

        vec3 weight;
        vec3 d = mat.sample(-cr.dir, nhit, randf(), randf(), weight);
        cr = Ray(phit, d);
    }

    return ret > 0;
}

static vec3 rayTraceRef(const Ray& primary_ray, const Scene& scene)
{
    HitInfo hit;
    Ray cr(primary_ray);
    vec3 radianceIntegral(0, 0, 0);
    vec3 throughput(1, 1, 1);
    const int maxDepth = 10;

    for (int n = 0; n <= maxDepth; n++) {
        hit = scene.intersect(cr);
        if (!hit) {
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
        const TriangleObject *obj = hit.triangleObject();
        vec3 nhit_g = obj->normal();

        auto& mat = todo_getMaterial(obj->materialID());
        vec3 weight;
        vec3 d = mat.sample(-cr.dir, nhit, randf(), randf(), weight);

        cr = Ray(phit, d);

        throughput = mult(throughput, weight);
    }

    if (radianceIntegral != radianceIntegral) radianceIntegral = vec3(0, 0, 0);
    return radianceIntegral;
}

void DemoBSSRDF::run()
{
    Scene scene;
    //     envmap.open("white.bin"); // works only with random variations added
    scene.setEnvmap(INPUT_DIR "envmap.bin"); //++++++++++++++++++++++++++
    //     envmap.open("envmap2.bin");
    //     envmap.open("at_the_window2.bin");
    //     envmap.open("sun.bin");
    //     envmap.open("sun2.bin");    envmap.scale(5);
    //     envmap.scale(2.5);
//     scene.setEnvmap("envmap.pfm", 512, 256);
    scene.envmap()->scale(1);

    initScene(scene);

    const int win_width = 256 * 4;
    const int win_height = 256 * 4;
    const int spp = 1;
//     const int win_width = 64;
//     const int win_height = 64;
//     const int spp = 64;

    FrameBuffer fb(win_width, win_height);
    FrameBuffer edge_map(win_width, win_height); // for edge anti-aliasing

    //////////////////////////////////////////////////////////////////////////
    float invWidth = 1 / float(win_width), invHeight = 1 / float(win_height);

    vec3 lookat;
    vec3 rayorig;
    {
        float zTrans = 7;
        float yTrans = 2.5;
        lookat = (scene.boxCenter());
        rayorig = lookat + vec3(0, yTrans*0.8, zTrans*0.5) - vec3(0, 2.1, 0); // no x-shift
        float viewdist = (lookat - rayorig).length();
        float _rot = -120 / 180.0 * NUM_PI;
        rayorig = vec3(-viewdist*sin(_rot), yTrans*0.5, viewdist*cos(_rot)) + lookat;

        lookat = vec3(0);
        rayorig = vec3(3.032326, 1.250000, -1.750714);
        {
            printf("lookat: %f, %f, %f\n", lookat.x, lookat.y, lookat.z);
            printf("rayorig: %f, %f, %f\n", rayorig.x, rayorig.y, rayorig.z);
        }
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
            Ray ray = cam.makeRay((x + 0.5) * invWidth, (y + 0.5) * invHeight);
            vec3 v = (vec3)renderEdge(ray, scene);
            edge_map.accumulatePixel(x, y, v);
        }
    }
    edge_map.tonemapGamma(2.2);
    edge_map.dumpPPM("edge.ppm");

    printBSSRDFParams();

#if 1
    for (int y = 0; y < win_height; ++y)
    {
        printf("\r render %d / %d  ", y + 1, win_height);
#pragma omp parallel for schedule(dynamic)
        for (int x = 0; x < win_width; ++x)
        {
            vec3 g;
            int samps = spp;
            if (edge_map.pixel(x, y).x > 0)
            {
                samps *= 16;
            }

            for (int s = 0; s < samps; ++s)
            {

                float dx = 0;
                float dy = 0;
                if (samps > 1)
                {
                    float r1 = 2 * randf();
                    float r2 = 2 * randf();
                    dx = r1 < 1 ? std::sqrt(r1) - 1 : 1 - std::sqrt(2 - r1);
                    dy = r2 < 1 ? std::sqrt(r2) - 1 : 1 - std::sqrt(2 - r2);
                }

                Ray ray = cam.makeRay((x + 0.5 + dx) * invWidth, (y + 0.5 + dy) * invHeight);

                vec3 v = radiance_sss(ray, scene);
                g += vec3(f_max(0.0f, v.x), f_max(0.0f, v.y), f_max(0.0f, v.z));
            }

            fb.accumulatePixel(x, y, g / float(samps));
        }
    }
#endif

    fb.dumpHDR("test_sss.hdr");

    //fb.tonemapReinhard();
    fb.tonemapGamma(2.2);

    float duration = tm.getTime();
    printf("took < %f > second\n", duration);

    fb.dumpPPM("test_sss.ppm");
}

int runTest(int argc, char *argv[])
{
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());
    omp_set_num_threads(std::max(1, omp_get_max_threads() - 1));
    //     omp_set_num_threads(1); // single thread
    printf("#procs: %d\n", omp_get_num_procs());
    printf("#threads: %d\n", omp_get_max_threads());


    DemoBSSRDF *demo = new DemoBSSRDF;

    demo->run();
    delete demo;

    return 0;
}

int main(int argc, char **argv)
{
    runTest(argc, argv);
}

