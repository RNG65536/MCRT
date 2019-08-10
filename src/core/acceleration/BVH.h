#pragma once

// median split BVH

#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "BoundedTriangle.h"
#include "Logger.h"
#include "aabb.h"
#include "constants.h"
#include "intersection.h"
#include "numeric.h"
#include "triangle.h"

namespace debug
{
struct BVHNode
{
    std::shared_ptr<AABB> boundingBox;

    uint32_t    debug_level;
    std::string debug_tag;

    // TODO : aovid this
    virtual bool intersect(Ray& ray, HitInfo& hit_info) const = 0;
};

struct BVHLeafNode : public BVHNode
{
    // std::shared_ptr<BoundedTriangle> primitive;
    std::shared_ptr<Triangle> primitive;

    bool intersect(Ray& ray, HitInfo& hit_info) const override;
};

struct BVHBranchNode : public BVHNode
{
    std::shared_ptr<BVHNode> leftChild;
    std::shared_ptr<BVHNode> rightChild;
    int                      split_axis;

    // TODO : remove redundant logics
    bool intersect(Ray& ray, HitInfo& hit_info) const override;
};

struct SplitAxisComparator
{
    SplitAxisComparator(int split_axis);

    bool operator()(const std::shared_ptr<Triangle>& a,
                    const std::shared_ptr<Triangle>& b) const;

    int m_split_axis;
};

class BVH
{
public:
    BVH(const std::vector<Triangle>& mesh);

    std::shared_ptr<BVHNode> build(int                begin,
                                   int                end,
                                   uint32_t           level,
                                   const std::string& direction);

    bool intersect(Ray& ray, HitInfo& hit_info) const;

private:
    std::vector<std::shared_ptr<Triangle> > m_triangles;
    std::shared_ptr<BVHNode>                m_rootNode;
};
}
