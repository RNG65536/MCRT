#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>
#include "constants.h"
#include "ray.h"
#include "scenebvh.h"

#define USE_EMBREE WITH_EMBREE

#if USE_EMBREE
#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>

static void errCB(const RTCError code, const char* str)
{
    const char* error = "";
    switch (code)
    {
        case RTC_UNKNOWN_ERROR:
        {
            error = "RTC_UNKNOWN_ERROR";
            break;
        }
        case RTC_INVALID_ARGUMENT:
        {
            error = "RTC_INVALID_ARGUMENT";
            break;
        }
        case RTC_INVALID_OPERATION:
        {
            error = "RTC_INVALID_OPERATION";
            break;
        }
        case RTC_OUT_OF_MEMORY:
        {
            error = "RTC_OUT_OF_MEMORY";
            break;
        }
        case RTC_UNSUPPORTED_CPU:
        {
            error = "RTC_UNSUPPORTED_CPU";
            break;
        }
        default:
        {
            error = "Invalid error code";
            break;
        }
    }
    printf("Embree error : %s\n", error);
}

HitInfo SceneBVH::intersect(const Ray& ray) const
{
    RTCRay rtc_ray;
    rtc_ray.org[0] = (float)(ray.orig[0]);
    rtc_ray.org[1] = (float)(ray.orig[1]);
    rtc_ray.org[2] = (float)(ray.orig[2]);
    rtc_ray.dir[0] = (float)(ray.dir[0]);
    rtc_ray.dir[1] = (float)(ray.dir[1]);
    rtc_ray.dir[2] = (float)(ray.dir[2]);
    rtc_ray.tnear = (float)(NUM_EPS_RAY);
    rtc_ray.tfar = (float)(NUM_INFINITY);
    rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
    rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;
    rtc_ray.instID = RTC_INVALID_GEOMETRY_ID;
    rtc_ray.mask = 0xFFFFFFFF;
    rtc_ray.time = 0;

    // intersection query
    assert(nullptr != m_scene);
    rtcIntersect(static_cast<RTCScene>(m_scene), rtc_ray);
    if ((unsigned int)(rtc_ray.geomID) == RTC_INVALID_GEOMETRY_ID)
    {
        return HitInfo();  // no hit
    }

    vec3 sn = m_triangles_ref[rtc_ray.primID].normal(rtc_ray.u, rtc_ray.v);

    HitInfo isect = HitInfo(rtc_ray.tfar,
                            rtc_ray.u,
                            rtc_ray.v,
                            sn,
                            &m_triangles_ref[rtc_ray.primID]);
    isect.setPrimitiveID(rtc_ray.primID);
    isect.setMaterialID(m_triangles_ref[rtc_ray.primID].materialID());

    return isect;
}

bool SceneBVH::occluded(const Ray& ray, float tfar) const
{
    RTCRay rtc_ray;
    rtc_ray.org[0] = (float)(ray.orig[0]);
    rtc_ray.org[1] = (float)(ray.orig[1]);
    rtc_ray.org[2] = (float)(ray.orig[2]);
    rtc_ray.dir[0] = (float)(ray.dir[0]);
    rtc_ray.dir[1] = (float)(ray.dir[1]);
    rtc_ray.dir[2] = (float)(ray.dir[2]);
    rtc_ray.tnear = (float)(NUM_EPS_RAY);
    rtc_ray.tfar = tfar - NUM_EPS_RAY;
    rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
    rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;
    rtc_ray.instID = RTC_INVALID_GEOMETRY_ID;
    rtc_ray.mask = 0xFFFFFFFF;
    rtc_ray.time = 0;

    // intersection query
    rtcOccluded(static_cast<RTCScene>(m_scene), rtc_ray);
    if ((unsigned int)(rtc_ray.geomID) == RTC_INVALID_GEOMETRY_ID)
    {
        return false;
    }
    else
    {
        return true;
    }
}

SceneBVH::~SceneBVH()
{
    if (m_scene)
    {
        rtcDeleteScene(static_cast<RTCScene>(m_scene));
        m_scene = nullptr;
    }
    rtcDeleteDevice(static_cast<RTCDevice>(m_device));

    std::cout << "EMBREE FREE" << std::endl;
}

SceneBVH::SceneBVH(const std::vector<TriangleObject>& triangles)
    : m_triangles_ref(triangles)
{
    auto rtc_device = rtcNewDevice(nullptr);
    // create scene
    auto rtc_scene = rtcDeviceNewScene(
        rtc_device, RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT1);
    m_device = rtc_device;
    m_scene = rtc_scene;

    rtcDeviceSetErrorFunction(rtc_device, errCB);

    std::cout << "EMBREE INIT" << std::endl;

    // add meshes to the scene
    {
        size_t num_triangles = m_triangles_ref.size();

        unsigned int geomID = rtcNewTriangleMesh(
            rtc_scene, RTC_GEOMETRY_STATIC, num_triangles, num_triangles * 3);

        // copy vertices & faces
        auto* mapped_positions = reinterpret_cast<float*>(
            rtcMapBuffer(rtc_scene, geomID, RTC_VERTEX_BUFFER));
        auto* mapped_faces = reinterpret_cast<int*>(
            rtcMapBuffer(rtc_scene, geomID, RTC_INDEX_BUFFER));

        for (int j = 0; j < num_triangles; j++)
        {
            // positions
            const vec3& p1 = m_triangles_ref[j].a;
            const vec3& p2 = m_triangles_ref[j].b;
            const vec3& p3 = m_triangles_ref[j].c;

            // store into mapped buffers
            int mi1 = 3 * (int)(j);
            int mi2 = 3 * (int)(j) + 1;
            int mi3 = 3 * (int)(j) + 2;
            mapped_faces[mi1] = mi1;
            mapped_faces[mi2] = mi2;
            mapped_faces[mi3] = mi3;
            for (int k = 0; k < 3; k++)
            {
                mapped_positions[4 * mi1 + k] = (float)(p1[k]);
                mapped_positions[4 * mi2 + k] = (float)(p2[k]);
                mapped_positions[4 * mi3 + k] = (float)(p3[k]);
            }
        }

        rtcUnmapBuffer(rtc_scene, geomID, RTC_VERTEX_BUFFER);
        rtcUnmapBuffer(rtc_scene, geomID, RTC_INDEX_BUFFER);
    }

    rtcCommit(rtc_scene);
}
#else
#include <nanort/nanort.h>

struct NanoRTContext
{
    nanort::BVHAccel<float>   g_accel;
    std::vector<float>        g_verts;
    std::vector<unsigned int> g_faces;
};

HitInfo SceneBVH::intersect(const Ray& ray) const
{
    NanoRTContext* ctx = static_cast<NanoRTContext*>(m_ctx);

    nanort::Ray<float> nano_ray;
    nano_ray.org[0] = (float)(ray.orig[0]);
    nano_ray.org[1] = (float)(ray.orig[1]);
    nano_ray.org[2] = (float)(ray.orig[2]);
    nano_ray.dir[0] = (float)(ray.dir[0]);
    nano_ray.dir[1] = (float)(ray.dir[1]);
    nano_ray.dir[2] = (float)(ray.dir[2]);
    nano_ray.min_t = (float)(NUM_EPS_RAY);
    nano_ray.max_t = (float)(NUM_INFINITY);

    nanort::TriangleIntersector<float> triangleIntersector(
        ctx->g_verts.data(), ctx->g_faces.data(), sizeof(float) * 3);
    nanort::TriangleIntersection<float> result;
    bool hit = ctx->g_accel.Traverse(nano_ray, triangleIntersector, &result);

    if (!hit)
    {
        return HitInfo();
    }

    // shading normal
    vec3 sn = m_triangles_ref[result.prim_id].normal(result.u, result.v);

    HitInfo isect = HitInfo(
        result.t, result.u, result.v, sn, &m_triangles_ref[result.prim_id]);
    isect.setPrimitiveID(result.prim_id);
    isect.setMaterialID(m_triangles_ref[result.prim_id].materialID());

    return isect;
}

bool SceneBVH::occluded(const Ray& ray, float tfar) const
{
    NanoRTContext* ctx = static_cast<NanoRTContext*>(m_ctx);

    nanort::Ray<float> nano_ray;
    nano_ray.org[0] = (float)(ray.orig[0]);
    nano_ray.org[1] = (float)(ray.orig[1]);
    nano_ray.org[2] = (float)(ray.orig[2]);
    nano_ray.dir[0] = (float)(ray.dir[0]);
    nano_ray.dir[1] = (float)(ray.dir[1]);
    nano_ray.dir[2] = (float)(ray.dir[2]);
    nano_ray.min_t = (float)(NUM_EPS_RAY);
    nano_ray.max_t = tfar - NUM_EPS_RAY;

    // Intersection query
    nanort::TriangleIntersector<float> triangle_intersector(
        ctx->g_verts.data(), ctx->g_faces.data(), sizeof(float) * 3);
    nanort::TriangleIntersection<float> isect;
    if (!ctx->g_accel.Traverse(nano_ray, triangle_intersector, &isect))
    {
        return false;
    }
    else
    {
        return true;
    }
}

SceneBVH::~SceneBVH()
{
    delete m_ctx;
    std::cout << "NANORT FREE" << std::endl;
}

SceneBVH::SceneBVH(const std::vector<TriangleObject>& triangles)
    : m_triangles_ref(triangles)
{
    NanoRTContext* ctx = new NanoRTContext();
    assert(ctx);
    m_ctx = ctx;

    nanort::BVHBuildOptions<float> build_options;
    build_options.cache_bbox = false;

    int num_triangles = static_cast<int>(triangles.size());

    // convert mesh
    ctx->g_verts.clear();
    ctx->g_faces.clear();
    for (int j = 0; j < num_triangles; j++)
    {
        // Positions
        const vec3& p1 = m_triangles_ref[j].a;
        const vec3& p2 = m_triangles_ref[j].b;
        const vec3& p3 = m_triangles_ref[j].c;

        ctx->g_verts.push_back(p1.x);
        ctx->g_verts.push_back(p1.y);
        ctx->g_verts.push_back(p1.z);

        ctx->g_verts.push_back(p2.x);
        ctx->g_verts.push_back(p2.y);
        ctx->g_verts.push_back(p2.z);

        ctx->g_verts.push_back(p3.x);
        ctx->g_verts.push_back(p3.y);
        ctx->g_verts.push_back(p3.z);

        ctx->g_faces.push_back(j * 3);
        ctx->g_faces.push_back(j * 3 + 1);
        ctx->g_faces.push_back(j * 3 + 2);
    }

    nanort::TriangleMesh<float> tri_mesh(
        ctx->g_verts.data(), ctx->g_faces.data(), sizeof(float) * 3);
    nanort::TriangleSAHPred<float> tri_pred(
        ctx->g_verts.data(), ctx->g_faces.data(), sizeof(float) * 3);

    bool ret =
        ctx->g_accel.Build(num_triangles, tri_mesh, tri_pred, build_options);
    assert(ret);

    std::cout << "NANORT INIT" << std::endl;
}
#endif