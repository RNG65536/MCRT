#pragma once

// SAH BVH using binning

#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "AABB.h"
#include "Logger.h"
#include "Triangle.h"
#include "constants.h"
#include "numeric.h"

namespace debug
{
#define DEBUG_BINS 0
#define DEBUG_INFO 0

constexpr int   k_leaf_node_size = 4;
constexpr int   k_num_binning_buckets = 16;  // 256;
constexpr float k_edge_threshold = 1e-3f;

struct SAHBVHNode
{
    std::shared_ptr<AABB> bounding_box;

    uint32_t    debug_level;
    std::string debug_tag;

    // TOOD : avoid this
    virtual bool intersect(Ray& ray, HitInfo& hit_info) const = 0;
};

struct SAHBVHLeafNode : public SAHBVHNode
{
    // std::shared_ptr<BoundedTriangle> primitive[k_leaf_node_size];
    std::shared_ptr<Triangle> primitive[k_leaf_node_size];

    bool intersect(Ray& ray, HitInfo& hit_info) const override;
};

struct SAHBVHBranchNode : public SAHBVHNode
{
    std::shared_ptr<SAHBVHNode> left_child;
    std::shared_ptr<SAHBVHNode> right_child;
    int                         split_axis;

    // TODO : remove redundant logics
    bool intersect(Ray& ray, HitInfo& hit_info) const override;
};

struct SAHSplitAxisComparator
{
    SAHSplitAxisComparator(int split_axis);

    bool operator()(const std::shared_ptr<Triangle>& a,
                    const std::shared_ptr<Triangle>& b) const;

    int m_split_axis;
};

struct MedianSplitter
{
    MedianSplitter(std::vector<std::shared_ptr<Triangle> >& triangles,
                   const AABB&                              node_bound,
                   int                                      begin,
                   int                                      end);

    float cost() const;

    void sort(int& axis, int& mid);

    float m_cost = NUM_INFINITY;
    int   m_split_axis = -1;
    int   m_mid = -1;

    std::vector<std::shared_ptr<Triangle> >& m_triangles;
    const AABB&                              m_node_bound;
    int                                      m_begin;
    int                                      m_end;
};

struct SAHBinningBucket
{
    AABB m_bounding_box;
    int  m_num_entries = 0;
    int  m_is_valid = 0;
};

struct SAHSplitter
{
    SAHSplitter(std::vector<std::shared_ptr<Triangle> >& triangles,
                const AABB&                              node_bound,
                int                                      begin,
                int                                      end);

    void sort(int& axis, int& mid);

    float m_min_cost;
    int   m_min_axis;  // min cost split axis
    int   m_min_bin;   // min cost split bin, index of the first bucket in the
                       // second child

    std::vector<std::shared_ptr<Triangle> >& m_triangles;
    const AABB&                              m_node_bound;
    int                                      m_begin;
    int                                      m_end;
    float                                    m_edge[3];
};

class SAHBVH
{
public:
    SAHBVH(const std::vector<Triangle>& mesh);

    std::shared_ptr<SAHBVHNode> buildLeaf(int                begin,
                                          int                end,
                                          uint32_t           level,
                                          const std::string& direction);

    std::shared_ptr<SAHBVHBranchNode> buildBranch(int                begin,
                                                  int                end,
                                                  uint32_t           level,
                                                  const std::string& direction,
                                                  int&               mid);

    std::shared_ptr<SAHBVHNode> build(int                begin,
                                      int                end,
                                      uint32_t           level,
                                      const std::string& direction);

    bool intersect(Ray& ray, HitInfo& hit_info) const;

private:
    std::vector<std::shared_ptr<Triangle> > m_triangles;
    std::shared_ptr<SAHBVHNode>             m_rootNode;
};
}