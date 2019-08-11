#include <algorithm>
#include "bvh.h"

bool BVHLeafNode::intersect(Ray& ray, HitInfo& hit_info) const
{
    bool hit = primitive->intersect(ray, hit_info);
    if (hit)
    {
        // Logger::debug() << "[ " << debug_level << " ][ " << debug_tag
        //                << " ]" << std::endl;
        // Logger::debug() << "    primitive t_near = " << t_near <<
        // std::endl;

        assert(hit_info.m_t < ray.tmax);
        ray.tmax = hit_info.m_t;
    }
    return hit;  // if no hit, returns infinity by default
}

bool BVHBranchNode::intersect(Ray& ray, HitInfo& hit_info) const
{
    hit_info.reset();
    // Logger::debug() << "[ " << debug_level << " ][ " << debug_tag << " ]"
    //                << std::endl;

    // float tNear = NUM_INFINITY;  // this is necessary

    // AABB culling
    bool hitBox = boundingBox->overlap(ray);
    if (!hitBox)
    {
        return false;
    }

    bool    anyHit = false;
    HitInfo hit_info_child;

    if (leftChild)
    {
        bool hit = leftChild->intersect(ray, hit_info_child);
        if (hit)
        {
            if (hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
            }
            anyHit = true;
        }

        // speeds up but not necessary
        if (hit_info.m_t < ray.tmax)
        {
            ray.tmax = hit_info.m_t;
        }
    }

    if (rightChild)
    {
        bool hit = rightChild->intersect(ray, hit_info_child);
        if (hit)
        {
            if (hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
            }
            anyHit = true;
        }

        // speeds up but not necessary
        if (hit_info.m_t < ray.tmax)
        {
            ray.tmax = hit_info.m_t;
        }
    }

    if (anyHit)
    {
        return true;
    }
    else
    {
        hit_info.reset();
        return false;
    }
}

SplitAxisComparator::SplitAxisComparator(int split_axis)
    : m_split_axis(split_axis)
{
}

bool SplitAxisComparator::operator()(const std::shared_ptr<Triangle>& a,
                                     const std::shared_ptr<Triangle>& b) const
{
    const float3& ca = a->baryCenter();
    const float3& cb = b->baryCenter();

    if (m_split_axis == 0)
    {
        return ca.x < cb.x;
    }
    else if (m_split_axis == 1)
    {
        return ca.y < cb.y;
    }
    else if (m_split_axis == 2)
    {
        return ca.z < cb.z;
    }
    else
    {
        assert(false);
        return false;
    }
}

bool BVH::intersect(Ray& ray, HitInfo& hit_info) const
{
    // BVH traversal
    hit_info.reset();
    // tNear = NUM_INFINITY;

    bool hit = m_rootNode->intersect(ray, hit_info);

    // if (hit && tNear > ray.tmin() && tNear < ray.tmax())
    if (hit)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::shared_ptr<BVHNode> BVH::build(int                begin,
                                    int                end,
                                    uint32_t           level,
                                    const std::string& direction)
{
    assert(begin < end);

    // if (end - begin < 2)  // leaf node
    if (end - begin == 1)  // leaf node
    {
        auto node = std::make_shared<BVHLeafNode>();
        auto bounding_box = std::make_shared<AABB>();

        // TODO : check how to handle out-of-range more robustly
        assert(begin < m_triangles.size());
        node->primitive = m_triangles[begin];

        node->debug_level = level;
        node->debug_tag = std::string("leaf ") + direction;

        bounding_box->expandBy(m_triangles[begin]->boundingBox());
        node->boundingBox = bounding_box;

        return node;
    }
    else
    {
        auto node = std::make_shared<BVHBranchNode>();
        auto bounding_box = std::make_shared<AABB>();

        for (int i = begin; i < end; ++i)
        {
            const auto& box = m_triangles[i]->boundingBox();
            bounding_box->expandBy(box);
        }
        int split_axis = bounding_box->largestDimension();

        // TODO : use nth_element
        std::sort(m_triangles.begin() + begin,
                  m_triangles.begin() + end,
                  SplitAxisComparator(split_axis));

        auto mid = begin + (end - begin) / 2;
        node->leftChild = build(begin, mid, level + 1, "left");
        node->rightChild = build(mid, end, level + 1, "right");

        node->debug_level = level;
        node->debug_tag = std::string("branch ") + direction;

#if 0  // debug rebuild bounding box
	            bounding_box = std::make_shared<AABB>();
	            bounding_box->expandBy(*node->leftChild->boundingBox);
	            bounding_box->expandBy(*node->rightChild->boundingBox);
#endif
        node->boundingBox = bounding_box;

        return node;
    }
}

BVH::BVH(const std::vector<Triangle>& mesh)
{
    int total = (int)mesh.size();
    m_triangles.clear();

    for (int i = 0; i < total; i++)
    {
        // TODO : store obj ID in triangle reference
        // or make sure this is always deep copy
        // copy constructor
        auto nt = std::make_shared<Triangle>(mesh[i]);
        nt->setPrimitiveID(i);
        m_triangles.push_back(nt);
    }

    m_rootNode = build(0, (int)m_triangles.size(), 0u, "root");
}
