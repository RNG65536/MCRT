#pragma once

// NVIDIA split BVH
// Spatial Splits in Bounding Volume Hierarchies

#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "AABB.h"
#include "BoundedTriangle.h"
#include "Logger.h"
#include "constants.h"

namespace debug
{
#define DEBUG_BINS 0
#define DEBUG_INFO 0

#define ENABLE_MEDIAN_OBJECT_SPLIT 1
#define ENABLE_SAH_OBJECT_SPLIT 1
#define ENABLE_SAH_SPATIAL_SPLIT 1

constexpr int k_x_axis = 0;
constexpr int k_y_axis = 1;
constexpr int k_z_axis = 2;
constexpr int k_sbvh_leaf_node_size = 4;
constexpr int k_sbvh_num_binning_buckets =
    16;  // 256;
constexpr float k_edge_split_threshold = 1e-3f;

struct SplitBVHNode
{
    std::shared_ptr<AABB> bounding_box;

    uint32_t    debug_level;
    std::string debug_tag;

    // TOOD : avoid this
    virtual bool intersect(Ray& ray, HitInfo& hit_info) const = 0;
};

struct SplitBVHLeafNode : public SplitBVHNode
{
    std::shared_ptr<BoundedTriangle> primitive[k_sbvh_leaf_node_size];

    bool intersect(Ray& ray, HitInfo& hit_info) const override;
};

struct SplitBVHBranchNode : public SplitBVHNode
{
    std::shared_ptr<SplitBVHNode> left_child;
    std::shared_ptr<SplitBVHNode> right_child;
    int                           split_axis;

    // TODO : remove redundant logics
    bool intersect(Ray& ray, HitInfo& hit_info) const override;
};

//////////////////////////////////////////////
// splitting techniques
//////////////////////////////////////////////

struct ClippedBoundCenterComparator
{
    ClippedBoundCenterComparator(int split_axis);

    bool operator()(const std::shared_ptr<BoundedTriangle>& a,
                    const std::shared_ptr<BoundedTriangle>& b);

    int m_split_axis;
};

struct MedianSplitBVHSplitter
{
    MedianSplitBVHSplitter(
        const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
        const AABB&                                           node_bound);

    float cost() const;

    void sort(int&                                            axis,
              std::vector<std::shared_ptr<BoundedTriangle> >& left_branch,
              std::vector<std::shared_ptr<BoundedTriangle> >& right_branch);

    float m_cost = NUM_INFINITY;
    int   m_split_axis = -1;
    int   m_mid = -1;

    std::vector<std::shared_ptr<BoundedTriangle> > m_triangles;
    const AABB&                                    m_node_bound;
};

struct SAHBVHBinningBucket
{
    AABB m_bounding_box;
    int  m_num_entries = 0;
    int  m_is_valid = 0;
};

struct SAHBVHSplitter
{
    SAHBVHSplitter(
        const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
        const AABB&                                           node_bound);

    float cost() const;

    void sort(int&                                            axis,
              std::vector<std::shared_ptr<BoundedTriangle> >& left_branch,
              std::vector<std::shared_ptr<BoundedTriangle> >& right_branch);

    float m_min_cost;
    int   m_min_axis;  // min cost split axis
    int   m_min_bin;   // min cost split bin, index of the first bucket in the
                       // second child

    std::vector<std::shared_ptr<BoundedTriangle> > m_triangles;
    const AABB&                                    m_node_bound;
    float                                          m_edge[3];
};

struct SplitBVHBinningBucket
{
    AABB m_bounding_box;
    int  m_num_entries = 0;
    int  m_num_exits = 0;
    int  m_is_valid = 0;
};

struct SplitBVHSplitter
{
    SplitBVHSplitter(
        const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
        const AABB&                                           node_bound);

    float cost() const;

    void sort(int&                                            axis,
              std::vector<std::shared_ptr<BoundedTriangle> >& left_branch,
              std::vector<std::shared_ptr<BoundedTriangle> >& right_branch);

    float m_min_cost;
    int   m_min_axis;  // min cost split axis
    int   m_min_bin;   // min cost split bin, index of the first bucket in the
                       // second child

    const std::vector<std::shared_ptr<BoundedTriangle> >& m_triangles;
    const AABB&                                           m_node_bound;
    float                                                 m_edge[3];
};

class SplitBVH
{
public:
    SplitBVH(const std::vector<Triangle>& mesh);

    static std::shared_ptr<SplitBVHNode> build(
        const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
        uint32_t                                              level,
        const std::string&                                    direction);

    bool intersect(Ray& ray, HitInfo& hit_info) const;

private:
    std::vector<std::shared_ptr<Triangle> >        m_triangles;
    std::vector<std::shared_ptr<BoundedTriangle> > m_triangle_refs;
    std::shared_ptr<SplitBVHNode>                  m_rootNode;
};
}