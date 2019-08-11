//#include <omp.h>
//#include <algorithm>
//#include <glm/gtc/matrix_transform.hpp>
//#include <iostream>
//#include "camera.h"
//#include "envmap.h"
//#include "film.h"
//#include "geometry.h"
//#include "intersection.h"
//#include "material.h"
//#include "mesh.h"
#include "numeric.h"
#include "ray.h"
#include "triangle.h"
//#include "scene.h"
//#include "timer.h"

#define INPUT_DIR "../input/"

class DemoSimple
{
public:
    void run();
};

void runBVHTest();

void DemoSimple::run()
{
    runBVHTest();
}

int runTest(int argc, char* argv[])
{
    DemoSimple* demo = new DemoSimple;
    demo->run();
    delete demo;

    return 0;
}
int main(int argc, char** argv)
{
    runTest(argc, argv);

    return 0;
}

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <string>

#include "acceleration/bvh.h"
#include "acceleration/sahbvh.h"
#include "acceleration/splitbvh.h"
#include "boundedtriangle.h"
#include "debugintersection.h"
#include "logger.h"

#include "debug_mesh.h"
bool checkEqual(float f1, float f2)
{
    return fabs(f1 - f2) < 1e-6f;
}

class IntersectionTriangle : public Intersection
{
public:
    IntersectionTriangle(const Triangle& triangle) : m_triangle(triangle)
    {
    }

    float intersectWith(const Ray& ray) const override
    {
        HitInfo hit_info;
        Ray     r(ray);
        m_triangle.intersect(r, hit_info);
        return hit_info.m_t;
    }

private:
    const Triangle& m_triangle;
};

class IntersectionBoundedTriangle : public Intersection
{
public:
    IntersectionBoundedTriangle(const BoundedTriangle& triangle)
        : m_triangle(triangle)
    {
    }

    float intersectWith(const Ray& ray) const override
    {
        HitInfo hit_info;
        Ray     r(ray);
        m_triangle.intersect(r, hit_info);
        return hit_info.m_t;
    }

private:
    const BoundedTriangle& m_triangle;
};

class IntersectionMesh : public Intersection
{
public:
    IntersectionMesh(const DebugMesh& mesh) : m_mesh(mesh)
    {
    }

    float intersectWith(const Ray& ray) const override
    {
        HitInfo hit_info;
        Ray     r(ray);
        m_mesh.intersect(r, hit_info);
        return hit_info.m_t;
    }

private:
    const DebugMesh& m_mesh;
};

class IntersectionNaiveBVH : public Intersection
{
public:
    IntersectionNaiveBVH(const BVH& bvh) : m_bvh(bvh)
    {
    }

    float intersectWith(const Ray& ray) const override
    {
        HitInfo hit_info;
        Ray     r(ray);
        m_bvh.intersect(r, hit_info);
        return hit_info.m_t;
    }

private:
    const BVH& m_bvh;
};

class IntersectionSAHBVH : public Intersection
{
public:
    IntersectionSAHBVH(const SAHBVH& bvh) : m_bvh(bvh)
    {
    }

    float intersectWith(const Ray& ray) const override
    {
        HitInfo hit_info;
        Ray     r(ray);
        m_bvh.intersect(r, hit_info);
        return hit_info.m_t;
    }

private:
    const SAHBVH& m_bvh;
};

class IntersectionSplitBVH : public Intersection
{
public:
    IntersectionSplitBVH(const SplitBVH& bvh) : m_bvh(bvh)
    {
    }

    float intersectWith(const Ray& ray) const override
    {
        HitInfo hit_info;
        Ray     r(ray);
        m_bvh.intersect(r, hit_info);
        return hit_info.m_t;
    }

private:
    const SplitBVH& m_bvh;
};

//////////////////////////////////////////

class Benchmark
{
public:
    Benchmark(const Intersection& reference,
              const Intersection& subject,
              const std::string&  name)
        : m_reference(reference), m_subject(subject), m_name(name)
    {
    }

    void run()
    {
        // benchmark build performance
        // exit(0);

        // while (0)
        int test_total = 1000000;
        int test_step = test_total / 100;
        for (int i = 0; i < test_total; i++)
        {
            float3 origin(randf() * 2.0f - 1.0f,
                          randf() * 2.0f - 1.0f,
                          randf() * 2.0f - 1.0f);
            float3 direction = normalize(
                float3(randf() - 0.5f, randf() - 0.5f, randf() - 0.5f));
            Ray ray(origin, direction);

            float ta = m_reference.intersectWith(ray);
            float tb = m_subject.intersectWith(ray);

            bool no_hit = std::isinf(ta) && std::isinf(tb);
            if (no_hit)
            {
            }
            else if (checkEqual(ta, tb))
            {
                // Logger::info()
                //    << "passed: ta = " << ta << ", tb = " << tb << std::endl;
            }
            else
            {
                Logger::info()
                    << "failed: ta = " << ta << ", tb = " << tb << std::endl;
                Logger::info() << "ray = " << ray.toString() << std::endl;
                Logger::info() << std::endl;

                // exit(1);
            }

            if (0 == i % test_step)
            {
                Logger::info()
                    << "checking " << i / test_step << "%" << std::endl;
            }
        }

        int                 n_buffer = 1000000;
        std::vector<float3> ray_origins(n_buffer);
        std::vector<float3> ray_directions(n_buffer);
        for (int i = 0; i < n_buffer; i++)
        {
            float3 origin(randfFast() * 2.0f - 1.0f,
                          randfFast() * 2.0f - 1.0f,
                          randfFast() * 2.0f - 1.0f);
            float3 direction = normalize(float3(
                randfFast() - 0.5f, randfFast() - 0.5f, randfFast() - 0.5f));
            ray_origins[i] = origin;
            ray_directions[i] = direction;
        }

        while (1)
        {
            auto begin = std::chrono::steady_clock::now();
#pragma omp parallel for
            for (int i = 0; i < m_n_rays; i++)
            {
                Ray ray(ray_origins[i % n_buffer],
                        ray_directions[i % n_buffer]);

                m_subject.intersectWith(ray);
            }
            auto end = std::chrono::steady_clock::now();
            auto duration =
                std::chrono::duration_cast<std::chrono::microseconds>(end -
                                                                      begin);
            auto  duration_f = float(duration.count());
            float mrays_per_sec = m_n_rays / duration_f;
            Logger::info() << m_name << " : " << mrays_per_sec
                           << " Mrays / s, with " << m_n_rays << " rays"
                           << std::endl;

            // update ray count to roughly specify interval
            float interval = 1.0f;
            m_n_rays = int(mrays_per_sec * 1e6f * interval);
        }
    }

private:
    const Intersection& m_reference;
    const Intersection& m_subject;
    std::string         m_name;
    int                 m_n_rays = 10000;
};

////////////////////////////////////////////////////////

void benchmarkLinearTraversal(const DebugMesh& mesh, const Intersection& ref)
{
    // median split bvh
    IntersectionMesh isect3(mesh);
    Benchmark(ref, isect3, "linear").run();
}

void benchmarkNaiveBVH(const DebugMesh& mesh, const Intersection& ref)
{
    // median split bvh
    BVH                  bvh(mesh.triangles());
    IntersectionNaiveBVH isect4(bvh);
    Benchmark(ref, isect4, "naive bvh").run();
}

void benchmarkSAHBVH(const DebugMesh& mesh, const Intersection& ref)
{
    // SAH bvh
    SAHBVH             sahbvh(mesh.triangles());
    IntersectionSAHBVH isect5(sahbvh);
    Benchmark(ref, isect5, "sah bvh").run();
}

void benchmarkSplitBVH(const DebugMesh& mesh, const Intersection& ref)
{
    // split bvh
    SplitBVH             splitbvh(mesh.triangles());
    IntersectionSplitBVH isect6(splitbvh);
    Benchmark(ref, isect6, "split bvh").run();
}

void runBVHTest()
{
    // naive scan
    DebugMesh        mesh;
    IntersectionMesh isect3(mesh);

    Logger::info() << "mesh AABB = " << mesh.boundingBox().toString()
                   << std::endl;
    Logger::info() << "#######################################################"
                   << std::endl;

    // benchmarkLinearTraversal(mesh, isect3);
    // benchmarkNaiveBVH(mesh, isect3);
    // benchmarkSAHBVH(mesh, isect3);
    benchmarkSplitBVH(mesh, isect3);
}
