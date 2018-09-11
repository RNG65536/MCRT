#include <cassert>
#include <iostream>
#include <string>
#include "constants.h"
#include "ray.h"
#include "scenebvh.h"

#define USE_EMBREE WITH_EMBREE

#if USE_EMBREE
#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>

void errCB(const RTCError code, const char* str)
{
    std::string error = "";
    switch (code)
    {
    case RTC_UNKNOWN_ERROR:     { error = "RTC_UNKNOWN_ERROR";      break; }
    case RTC_INVALID_ARGUMENT:  { error = "RTC_INVALID_ARGUMENT";   break; }
    case RTC_INVALID_OPERATION: { error = "RTC_INVALID_OPERATION";  break; }
    case RTC_OUT_OF_MEMORY:     { error = "RTC_OUT_OF_MEMORY";      break; }
    case RTC_UNSUPPORTED_CPU:   { error = "RTC_UNSUPPORTED_CPU";    break; }
    default:                    { error = "Invalid error code";	    break; }
    }
    std::cout << ("Embree error : " + error) << std::endl;
}

HitInfo SceneBVH::intersect(const Ray& ray) const
{
    //         if (minT > maxT)
    //         {
    //             return false;
    //         }

    // Create RTCRay
    RTCRay rtcRay;
    rtcRay.org[0] = (float)(ray.o[0]);
    rtcRay.org[1] = (float)(ray.o[1]);
    rtcRay.org[2] = (float)(ray.o[2]);
    rtcRay.dir[0] = (float)(ray.d[0]);
    rtcRay.dir[1] = (float)(ray.d[1]);
    rtcRay.dir[2] = (float)(ray.d[2]);
    rtcRay.tnear = (float)(NUM_EPS_RAY);
    rtcRay.tfar = (float)(NUM_INFINITY);
    rtcRay.geomID = RTC_INVALID_GEOMETRY_ID;
    rtcRay.primID = RTC_INVALID_GEOMETRY_ID;
    rtcRay.instID = RTC_INVALID_GEOMETRY_ID;
    rtcRay.mask = 0xFFFFFFFF;
    rtcRay.time = 0;

    // Intersection query
    assert(nullptr != RtcScene);
    rtcIntersect(static_cast<RTCScene>(RtcScene), rtcRay);
    if ((unsigned int)(rtcRay.geomID) == RTC_INVALID_GEOMETRY_ID)
    {
        //             cout << "no hit" << endl;
        return HitInfo();
    }

    //         cout << rtcRay.tfar << endl;


    // Fill in the intersection structure
    //         Hit_CUDA isect = Hit_CUDA(
    //             rtcRay.tfar,
    //             vec3(-rtcRay.Ng[0], -rtcRay.Ng[1], -rtcRay.Ng[2]).normalize(),
    //             &triangleObjectPool[rtcRay.primID]);

    vec3 sn = triangleObjectPool[rtcRay.primID].normal(rtcRay.u, rtcRay.v);
    //         vec3 gn = vec3(-rtcRay.Ng[0], -rtcRay.Ng[1], -rtcRay.Ng[2]).normalize();
    //         if (dot(gn, sn) <= 0) sn = gn;
//     vec3 sn = triangleObjectPool[rtcRay.primID].getFaceNormal();
//     if (dot(ray.d, sn) >= 0)
//     {
//         sn = -sn;
//     }

    HitInfo isect = HitInfo(
        rtcRay.tfar,
        rtcRay.u, rtcRay.v,
        sn,
        &triangleObjectPool[rtcRay.primID]);
    isect.setPrimitiveID(rtcRay.primID);
    isect.setMaterialID(triangleObjectPool[rtcRay.primID].materialID());

    return isect;
}

bool SceneBVH::occluded(const Ray& ray, float tfar) const
{
    RTCRay rtcRay;
    rtcRay.org[0] = (float)(ray.o[0]);
    rtcRay.org[1] = (float)(ray.o[1]);
    rtcRay.org[2] = (float)(ray.o[2]);
    rtcRay.dir[0] = (float)(ray.d[0]);
    rtcRay.dir[1] = (float)(ray.d[1]);
    rtcRay.dir[2] = (float)(ray.d[2]);
    rtcRay.tnear = (float)(NUM_EPS_RAY * 1);
    rtcRay.tfar = tfar - NUM_EPS_RAY * 1; // /*(float)(NUM_INFINITY)*/;
    rtcRay.geomID = RTC_INVALID_GEOMETRY_ID;
    rtcRay.primID = RTC_INVALID_GEOMETRY_ID;
    rtcRay.instID = RTC_INVALID_GEOMETRY_ID;
    rtcRay.mask = 0xFFFFFFFF;
    rtcRay.time = 0;

    // Intersection query
    rtcOccluded(static_cast<RTCScene>(RtcScene), rtcRay);
    if ((unsigned int)(rtcRay.geomID) == RTC_INVALID_GEOMETRY_ID)
    {
        return false;
    }
    else {
        return true;
    }
}

SceneBVH::~SceneBVH()
{
    if (RtcScene)
    {
        rtcDeleteScene(static_cast<RTCScene>(RtcScene));
        RtcScene = nullptr;
    }
    rtcDeleteDevice(static_cast<RTCDevice>(device));
    std::cout << "EMBREE FREE" << std::endl;
}

SceneBVH::SceneBVH(const std::vector<TriangleObject>& triangles) : triangleObjectPool(triangles)
{
    device = rtcNewDevice(nullptr);
    rtcDeviceSetErrorFunction(static_cast<RTCDevice>(device), errCB);
    std::cout << "EMBREE INIT" << std::endl;

    // Create scene
    RtcScene = rtcDeviceNewScene(static_cast<RTCDevice>(device), RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT1);

    // Add meshes to the scene
    {
        size_t numFaces = triangleObjectPool.size();

        unsigned int geomID = rtcNewTriangleMesh(static_cast<RTCScene>(RtcScene), RTC_GEOMETRY_STATIC, numFaces, numFaces * 3);

        // Copy vertices & faces
        auto* mappedPositions = reinterpret_cast<float*>(rtcMapBuffer(static_cast<RTCScene>(RtcScene), geomID, RTC_VERTEX_BUFFER));
        auto* mappedFaces = reinterpret_cast<int*>(rtcMapBuffer(static_cast<RTCScene>(RtcScene), geomID, RTC_INDEX_BUFFER));

        for (int j = 0; j < numFaces; j++)
        {
            // Positions
            const vec3& p1 = triangleObjectPool[j].a;
            const vec3& p2 = triangleObjectPool[j].b;
            const vec3& p3 = triangleObjectPool[j].c;

            // Store into mapped buffers
            int mi1 = 3 * (int)(j);
            int mi2 = 3 * (int)(j)+1;
            int mi3 = 3 * (int)(j)+2;
            mappedFaces[mi1] = mi1;
            mappedFaces[mi2] = mi2;
            mappedFaces[mi3] = mi3;
            for (int k = 0; k < 3; k++)
            {
                mappedPositions[4 * mi1 + k] = (float)(p1[k]);
                mappedPositions[4 * mi2 + k] = (float)(p2[k]);
                mappedPositions[4 * mi3 + k] = (float)(p3[k]);
            }
        }

        rtcUnmapBuffer(static_cast<RTCScene>(RtcScene), geomID, RTC_VERTEX_BUFFER);
        rtcUnmapBuffer(static_cast<RTCScene>(RtcScene), geomID, RTC_INDEX_BUFFER);
    }

    rtcCommit(static_cast<RTCScene>(RtcScene));
}
#else
#include <nanort/nanort.h>

nanort::BVHAccel<float> *g_accel = nullptr;
std::vector<float> g_verts;
std::vector<unsigned int> g_faces;

HitInfo SceneBVH::intersect(const Ray& ray) const
{
    nanort::Ray<float> nanoRay;
    nanoRay.org[0] = (float)(ray.o[0]);
    nanoRay.org[1] = (float)(ray.o[1]);
    nanoRay.org[2] = (float)(ray.o[2]);
    nanoRay.dir[0] = (float)(ray.d[0]);
    nanoRay.dir[1] = (float)(ray.d[1]);
    nanoRay.dir[2] = (float)(ray.d[2]);
    nanoRay.min_t = (float)(NUM_EPS_RAY);
    nanoRay.max_t = (float)(NUM_INFINITY);

    nanort::TriangleIntersector<float> triangleIntersector(g_verts.data(), g_faces.data(), sizeof(float) * 3);
    nanort::TriangleIntersection<float> result;
    bool hit = g_accel->Traverse(nanoRay, triangleIntersector, &result);

    if (!hit)
    {
        return HitInfo();
    }

    vec3 sn = triangleObjectPool[result.prim_id].normal(result.u, result.v);

    HitInfo isect = HitInfo(
        result.t,
        result.u, result.v,
        sn,
        &triangleObjectPool[result.prim_id]);
    isect.setPrimitiveID(result.prim_id);
    isect.setMaterialID(triangleObjectPool[result.prim_id].materialID());

    return isect;
}

bool SceneBVH::occluded(const Ray& ray, float tfar) const
{
    nanort::Ray<float> nanoRay;
    nanoRay.org[0] = (float)(ray.o[0]);
    nanoRay.org[1] = (float)(ray.o[1]);
    nanoRay.org[2] = (float)(ray.o[2]);
    nanoRay.dir[0] = (float)(ray.d[0]);
    nanoRay.dir[1] = (float)(ray.d[1]);
    nanoRay.dir[2] = (float)(ray.d[2]);
    nanoRay.min_t = (float)(NUM_EPS_RAY * 1);
    nanoRay.max_t = tfar - NUM_EPS_RAY * 1;

    // Intersection query
    nanort::TriangleIntersector<float> triangleIntersector(g_verts.data(), g_faces.data(), sizeof(float) * 3);
    nanort::TriangleIntersection<float> isect;
    if (!g_accel->Traverse(nanoRay, triangleIntersector, &isect))
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
    delete g_accel;
    std::cout << "NANORT FREE" << std::endl;
}

SceneBVH::SceneBVH(const std::vector<TriangleObject>& triangles) : triangleObjectPool(triangles)
{
    nanort::BVHBuildOptions<float> build_options;
    build_options.cache_bbox = false;

    size_t numTriangles = triangles.size();

    // convert mesh
    g_verts.clear();
    g_faces.clear();
    for (int j = 0; j < numTriangles; j++)
    {
        // Positions
        const vec3& p1 = triangleObjectPool[j].a;
        const vec3& p2 = triangleObjectPool[j].b;
        const vec3& p3 = triangleObjectPool[j].c;

        g_verts.push_back(p1.x);
        g_verts.push_back(p1.y);
        g_verts.push_back(p1.z);

        g_verts.push_back(p2.x);
        g_verts.push_back(p2.y);
        g_verts.push_back(p2.z);

        g_verts.push_back(p3.x);
        g_verts.push_back(p3.y);
        g_verts.push_back(p3.z);

        g_faces.push_back(j * 3);
        g_faces.push_back(j * 3 + 1);
        g_faces.push_back(j * 3 + 2);
    }

    nanort::TriangleMesh<float> triMesh(g_verts.data(), g_faces.data(), sizeof(float) * 3);
    nanort::TriangleSAHPred<float> triPred(g_verts.data(), g_faces.data(), sizeof(float) * 3);

    g_accel = new nanort::BVHAccel<float>();
    bool ret = g_accel->Build(numTriangles, triMesh, triPred, build_options);
    assert(ret);

    std::cout << "NANORT INIT" << std::endl;

}
#endif