#include <algorithm>
#include "SplitBVH.h"
#include "numeric.h"

namespace debug
{
bool SplitBVHLeafNode::intersect(Ray& ray, HitInfo& hit_info) const
{
    hit_info.reset();

    bool hitBox = bounding_box->overlap(ray);
    if (!hitBox)
    {
        return false;
    }

    bool    any_hit = false;
    HitInfo hit_info_child;

    for (int i = 0; i < k_sbvh_leaf_node_size; i++)
    {
        if (nullptr == primitive[i])
        {
            continue;
        }

        bool hitClippedBound = primitive[i]->intersectBound(ray);
        if (!hitClippedBound)
        {
            continue;
        }

        bool hit = primitive[i]->intersect(ray, hit_info_child);
        if (hit)
        {
            //#if DEBUG_INFO
            //            Logger::debug() << "[ " << debug_level << " ][ " <<
            //            debug_tag
            //                            << " ]" << std::endl;
            //            Logger::debug() << "    primitive t_near = " << t_near
            //            << std::endl;
            //#endif

            assert(hit_info_child.m_t > ray.tmin &&
                   hit_info_child.m_t < ray.tmax);

            if (hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
                any_hit = true;
                ray.tmax = hit_info.m_t;
            }
        }
    }
    return any_hit;  // if no hit, returns infinity by default
}

bool SplitBVHBranchNode::intersect(Ray& ray, HitInfo& hit_info) const
{
    // Logger::debug() << "[ " << debug_level << " ][ " << debug_tag << " ]"
    //                << std::endl;

    hit_info.reset();

    // AABB culling
    bool hitBox = bounding_box->overlap(ray);
    if (!hitBox)
    {
        return false;
    }

    bool    anyHit = false;
    HitInfo hit_info_child;

    float ray_d[3] = {ray.dir.x, ray.dir.y, ray.dir.z};

    if (ray_d[split_axis] > 0.0f)
    {
        if (left_child)
        {
            bool hit = left_child->intersect(ray, hit_info_child);
            if (hit && hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
                anyHit = true;

                // speeds up but not necessary
                assert(hit_info.m_t < ray.tmax);
                ray.tmax = hit_info.m_t;
            }
        }

        if (right_child)
        {
            bool hit = right_child->intersect(ray, hit_info_child);
            if (hit && hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
                anyHit = true;

                // speeds up but not necessary
                assert(hit_info.m_t < ray.tmax);
                ray.tmax = hit_info.m_t;
            }
        }
    }
    else
    {
        if (right_child)
        {
            bool hit = right_child->intersect(ray, hit_info_child);
            if (hit && hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
                anyHit = true;

                // speeds up but not necessary
                assert(hit_info.m_t < ray.tmax);
                ray.tmax = hit_info.m_t;
            }
        }

        if (left_child)
        {
            bool hit = left_child->intersect(ray, hit_info_child);
            if (hit && hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
                anyHit = true;

                // speeds up but not necessary
                assert(hit_info.m_t < ray.tmax);
                ray.tmax = hit_info.m_t;
            }
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

ClippedBoundCenterComparator::ClippedBoundCenterComparator(int split_axis)
    : m_split_axis(split_axis)
{
}

bool ClippedBoundCenterComparator::operator()(
    const std::shared_ptr<BoundedTriangle>& a,
    const std::shared_ptr<BoundedTriangle>& b)
{
    const float3& ca = a->clippedBoundingBox().center();
    const float3& cb = b->clippedBoundingBox().center();

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

MedianSplitBVHSplitter::MedianSplitBVHSplitter(
    const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
    const AABB&                                           node_bound)
    : m_triangles(triangles), m_node_bound(node_bound)
{
    AABB left_bound;
    AABB right_bound;
    m_split_axis = m_node_bound.largestDimension();
    std::sort(m_triangles.begin(),
              m_triangles.end(),
              ClippedBoundCenterComparator(m_split_axis));
    m_mid = (int)m_triangles.size() / 2;
    int total = (int)m_triangles.size();
    for (int i = 0; i < m_mid; i++)
    {
        left_bound.expandBy(m_triangles[i]->clippedBoundingBox());
    }
    for (int i = m_mid; i < total; i++)
    {
        right_bound.expandBy(m_triangles[i]->clippedBoundingBox());
    }
    int left_count = m_mid;
    int right_count = total - m_mid;
    m_cost = left_bound.surfaceAreaCost() * left_count +
             right_bound.surfaceAreaCost() * right_count;

#if DEBUG_INFO
    char* axes[3] = {"x", "y", "z"};
    Logger::info() << "<median split bvh> num primitives = " << total
                   << ", min cost = " << m_cost
                   << ", min axis = " << axes[m_split_axis] << std::endl;
#endif
}

float MedianSplitBVHSplitter::cost() const
{
    return m_cost;
}

void MedianSplitBVHSplitter::sort(
    int&                                            axis,
    std::vector<std::shared_ptr<BoundedTriangle> >& left_branch,
    std::vector<std::shared_ptr<BoundedTriangle> >& right_branch)
{
    int total = (int)m_triangles.size();

    left_branch.clear();
    right_branch.clear();
    for (int i = 0; i < m_mid; i++)
    {
        left_branch.push_back(m_triangles[i]);
    }
    for (int i = m_mid; i < total; i++)
    {
        right_branch.push_back(m_triangles[i]);
    }
    axis = m_split_axis;
    return;
}

SAHBVHSplitter::SAHBVHSplitter(
    const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
    const AABB&                                           node_bound)
    : m_triangles(triangles), m_node_bound(node_bound)
{
    // AABB node_bound;
    // for (int i = begin; i < end; ++i)
    //{
    //    const auto& box = triangles[i]->boundingBox();
    //    node_bound.expandBy(box);
    //}
    float3 edge = node_bound.max() - node_bound.min();
    m_edge[0] = edge.x;
    m_edge[1] = edge.y;
    m_edge[2] = edge.z;

    m_min_cost = NUM_INFINITY;
    m_min_axis = -1;
    m_min_bin = -1;

    int total = (int)m_triangles.size();

    ////////////// debug ////////////////
#if DEBUG_BINS
    AABB  debug_x_aabb[k_sbvh_num_binning_buckets];
    AABB  debug_y_aabb[k_sbvh_num_binning_buckets];
    AABB  debug_z_aabb[k_sbvh_num_binning_buckets];
    int   debug_x_count[k_sbvh_num_binning_buckets] = {0};
    int   debug_y_count[k_sbvh_num_binning_buckets] = {0};
    int   debug_z_count[k_sbvh_num_binning_buckets] = {0};
    float debug_x_min_cost = NUM_INFINITY;
    float debug_y_min_cost = NUM_INFINITY;
    float debug_z_min_cost = NUM_INFINITY;
#endif

    /////////////////// check over x ////////////////////////
    if (edge.x > k_edge_split_threshold)
    {
        SAHBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_x = edge.x / k_sbvh_num_binning_buckets;
        // scan triangles
        for (int i = 0; i < total; i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();
            float       x0 = box.min().x;
            assert(x0 >= node_bound.min().x);
            assert(x0 <= node_bound.max().x);
            int bucket_index = i_min(int((x0 - node_bound.min().x) / delta_x),
                                     k_sbvh_num_binning_buckets - 1);

            assert(bucket_index < k_sbvh_num_binning_buckets);
            assert(bucket_index >= 0);

            bins[bucket_index].m_num_entries++;
            bins[bucket_index].m_bounding_box.expandBy(box);
            bins[bucket_index].m_is_valid = 1;
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_x_aabb[i] = bins[i].m_bounding_box;
            debug_x_count[i] = bins[i].m_num_entries;
        }
#endif

        // scan bins
        for (int i = 1; i < k_sbvh_num_binning_buckets; i++)
        {
            AABB bound_left;
            int  count_left = 0;
            for (int j = 0; j < i; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_left += bins[j].m_num_entries;
                    bound_left.expandBy(bins[j].m_bounding_box);
                }
            }

            // early skip
            float cost = 0.0f;
            if (count_left > 0)
            {
                cost += bound_left.surfaceAreaCost() * count_left;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost >= m_min_cost)
            {
                continue;
            }

            AABB bound_right;
            int  count_right = 0;
            for (int j = i; j < k_sbvh_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_entries;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceAreaCost() * count_right;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost <= m_min_cost)
            {
                m_min_cost = cost;
                m_min_axis = 0;
                m_min_bin = i;
            }

#if DEBUG_BINS
            if (cost <= debug_x_min_cost)
            {
                debug_x_min_cost = cost;
            }
#endif
        }
    }

    /////////////////// check over y ////////////////////////
    if (edge.y > k_edge_split_threshold)
    {
        SAHBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_y = edge.y / k_sbvh_num_binning_buckets;
        // scan triangles
        for (int i = 0; i < total; i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();
            float       y0 = box.min().y;
            assert(y0 >= node_bound.min().y);
            assert(y0 <= node_bound.max().y);
            int bucket_index = i_min(int((y0 - node_bound.min().y) / delta_y),
                                     k_sbvh_num_binning_buckets - 1);

            assert(bucket_index < k_sbvh_num_binning_buckets);
            assert(bucket_index >= 0);

            bins[bucket_index].m_num_entries++;
            bins[bucket_index].m_bounding_box.expandBy(box);
            bins[bucket_index].m_is_valid = 1;
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_y_aabb[i] = bins[i].m_bounding_box;
            debug_y_count[i] = bins[i].m_num_entries;
        }
#endif

        // scan bins
        for (int i = 1; i < k_sbvh_num_binning_buckets; i++)
        {
            AABB bound_left;
            int  count_left = 0;
            for (int j = 0; j < i; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_left += bins[j].m_num_entries;
                    bound_left.expandBy(bins[j].m_bounding_box);
                }
            }

            // early skip
            float cost = 0.0f;
            if (count_left > 0)
            {
                cost += bound_left.surfaceAreaCost() * count_left;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost >= m_min_cost)
            {
                continue;
            }

            AABB bound_right;
            int  count_right = 0;
            for (int j = i; j < k_sbvh_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_entries;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceAreaCost() * count_right;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost <= m_min_cost)
            {
                m_min_cost = cost;
                m_min_axis = 1;
                m_min_bin = i;
            }

#if DEBUG_BINS
            if (cost <= debug_y_min_cost)
            {
                debug_y_min_cost = cost;
            }
#endif
        }
    }

    /////////////////// check over z ////////////////////////
    if (edge.z > k_edge_split_threshold)
    {
        SAHBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_z = edge.z / k_sbvh_num_binning_buckets;
        // scan triangles
        for (int i = 0; i < total; i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();
            float       z0 = box.min().z;
            assert(z0 >= node_bound.min().z);
            assert(z0 <= node_bound.max().z);
            int bucket_index = i_min(int((z0 - node_bound.min().z) / delta_z),
                                     k_sbvh_num_binning_buckets - 1);

            assert(bucket_index < k_sbvh_num_binning_buckets);
            assert(bucket_index >= 0);

            bins[bucket_index].m_num_entries++;
            bins[bucket_index].m_bounding_box.expandBy(box);
            bins[bucket_index].m_is_valid = 1;
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_z_aabb[i] = bins[i].m_bounding_box;
            debug_z_count[i] = bins[i].m_num_entries;
        }
#endif

        // scan bins
        for (int i = 1; i < k_sbvh_num_binning_buckets; i++)
        {
            AABB bound_left;
            int  count_left = 0;
            for (int j = 0; j < i; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_left += bins[j].m_num_entries;
                    bound_left.expandBy(bins[j].m_bounding_box);
                }
            }

            // early skip
            float cost = 0.0f;
            if (count_left > 0)
            {
                cost += bound_left.surfaceAreaCost() * count_left;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost >= m_min_cost)
            {
                continue;
            }

            AABB bound_right;
            int  count_right = 0;
            for (int j = i; j < k_sbvh_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_entries;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceAreaCost() * count_right;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost <= m_min_cost)
            {
                m_min_cost = cost;
                m_min_axis = 2;
                m_min_bin = i;
            }

#if DEBUG_BINS
            if (cost <= debug_z_min_cost)
            {
                debug_z_min_cost = cost;
            }
#endif
        }
    }

        ////////////// debug ////////////////
#if DEBUG_BINS
    Logger::info() << "================================================"
                   << std::endl;
#endif

#if DEBUG_INFO
    if (m_min_axis >= 0)
    {
        char* axes[3] = {"x", "y", "z"};
        Logger::info() << "<sah bvh> num primitives = " << total
                       << ", min cost = " << m_min_cost
                       << ", min axis = " << axes[m_min_axis]
                       << ", min bin = " << m_min_bin << std::endl;
    }
    else
    {
        Logger::info() << "<sah bvh> num primitives = " << total
                       << ", cannot find min with edge split threshold "
                       << k_edge_split_threshold << std::endl;
    }
#endif

////////////// debug ////////////////
#if DEBUG_BINS
    if (m_min_bin == 1 || m_min_bin == k_sbvh_num_binning_buckets - 1)
    {
        // x
        Logger::info() << ">> x bins min cost " << debug_x_min_cost
                       << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_x_count[i] << ", "
                           << debug_x_aabb[i].toString() << std::endl;
        }
        // y
        Logger::info() << ">> y bins min cost " << debug_y_min_cost
                       << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_y_count[i] << ", "
                           << debug_y_aabb[i].toString() << std::endl;
        }
        // z
        Logger::info() << ">> z bins min cost " << debug_z_min_cost
                       << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_z_count[i] << ", "
                           << debug_z_aabb[i].toString() << std::endl;
        }
    }
#endif
}

float SAHBVHSplitter::cost() const
{
    return m_min_cost;
}

void SAHBVHSplitter::sort(
    int&                                            axis,
    std::vector<std::shared_ptr<BoundedTriangle> >& left_branch,
    std::vector<std::shared_ptr<BoundedTriangle> >& right_branch)
{
    axis = m_min_axis;

    left_branch.clear();
    right_branch.clear();
    int total = (int)m_triangles.size();

    SAHBVHBinningBucket bins[k_sbvh_num_binning_buckets];

    float delta = m_edge[m_min_axis] / k_sbvh_num_binning_buckets;
    float node_min_v[3] = {
        m_node_bound.min().x, m_node_bound.min().y, m_node_bound.min().z};

    // scan triangles
    for (int i = 0; i < total; i++)
    {
        const auto& box = m_triangles[i]->clippedBoundingBox();
        float       min_v[3] = {box.min().x, box.min().y, box.min().z};
        int         bucket_index =
            i_min(int((min_v[m_min_axis] - node_min_v[m_min_axis]) / delta),
                  k_sbvh_num_binning_buckets - 1);

        assert(bucket_index < k_sbvh_num_binning_buckets);
        assert(bucket_index >= 0);

        if (bucket_index < m_min_bin)
        {
            left_branch.push_back(m_triangles[i]);
        }
        else
        {
            right_branch.push_back(m_triangles[i]);
        }
    }
}

SplitBVHSplitter::SplitBVHSplitter(
    const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
    const AABB&                                           node_bound)
    : m_triangles(triangles), m_node_bound(node_bound)
{
    // AABB node_bound;
    // for (int i = begin; i < end; ++i)
    //{
    //    const auto& box = triangles[i]->boundingBox();
    //    node_bound.expandBy(box);
    //}
    float3 edge = m_node_bound.max() - m_node_bound.min();
    m_edge[0] = edge.x;
    m_edge[1] = edge.y;
    m_edge[2] = edge.z;

    m_min_cost = NUM_INFINITY;
    m_min_axis = -1;
    m_min_bin = -1;

    int total = (int)m_triangles.size();

    ////////////// debug ////////////////
#if DEBUG_BINS
    AABB  debug_x_aabb[k_sbvh_num_binning_buckets];
    AABB  debug_y_aabb[k_sbvh_num_binning_buckets];
    AABB  debug_z_aabb[k_sbvh_num_binning_buckets];
    int   debug_x_entries[k_sbvh_num_binning_buckets] = {0};
    int   debug_y_entries[k_sbvh_num_binning_buckets] = {0};
    int   debug_z_entries[k_sbvh_num_binning_buckets] = {0};
    int   debug_x_exits[k_sbvh_num_binning_buckets] = {0};
    int   debug_y_exits[k_sbvh_num_binning_buckets] = {0};
    int   debug_z_exits[k_sbvh_num_binning_buckets] = {0};
    float debug_x_min_cost = NUM_INFINITY;
    float debug_y_min_cost = NUM_INFINITY;
    float debug_z_min_cost = NUM_INFINITY;
    int   debug_x_min_bin = -1;
    int   debug_y_min_bin = -1;
    int   debug_z_min_bin = -1;
    int   debug_x_min_left_count = -1;
    int   debug_y_min_left_count = -1;
    int   debug_z_min_left_count = -1;
    int   debug_x_min_right_count = -1;
    int   debug_y_min_right_count = -1;
    int   debug_z_min_right_count = -1;
#endif

    /////////////////// check over x ////////////////////////
    if (edge.x > k_edge_split_threshold)
    {
        SplitBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_x = edge.x / k_sbvh_num_binning_buckets;
        float node_min = m_node_bound.min().x;
        float node_max = m_node_bound.max().x;
        // scan triangles
        for (int i = 0; i < total; i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();

            // entry bucket index
            float x0 = box.min().x;
            assert(x0 >= node_min);
            assert(x0 <= node_max);
            int i0 = i_min(int((x0 - node_min) / delta_x),
                           k_sbvh_num_binning_buckets - 1);
            assert(i0 < k_sbvh_num_binning_buckets);
            assert(i0 >= 0);

            // exit bucket index
            float x1 = box.max().x;
            assert(x1 >= node_min);
            assert(x1 <= node_max);
            int i1 = i_min(int((x1 - node_min) / delta_x),
                           k_sbvh_num_binning_buckets - 1);
            assert(i1 < k_sbvh_num_binning_buckets);
            assert(i1 >= 0);

            AABB  clipped_bbox;
            float bin_min;

            // update entry bin
            bin_min = node_min + i0 * delta_x;
            clipped_bbox =
                m_triangles[i]->clipBoundTight(bin_min, bin_min + delta_x, 0);
            if (clipped_bbox.isValid())
            {
                bins[i0].m_num_entries++;
                bins[i0].m_bounding_box.expandBy(clipped_bbox);
                bins[i0].m_is_valid = 1;
            }

            // update middle bins
            for (int j = i0 + 1; j < i1; j++)
            {
                bin_min = node_min + j * delta_x;
                clipped_bbox = m_triangles[i]->clipBoundTight(
                    bin_min, bin_min + delta_x, 0);
                if (clipped_bbox.isValid())
                {
                    bins[j].m_bounding_box.expandBy(clipped_bbox);
                    bins[j].m_is_valid = 1;
                }
            }
            // update exit bin
            bin_min = node_min + i1 * delta_x;
            clipped_bbox =
                m_triangles[i]->clipBoundTight(bin_min, bin_min + delta_x, 0);
            if (clipped_bbox.isValid())
            {
                bins[i1].m_num_exits++;
                bins[i1].m_bounding_box.expandBy(clipped_bbox);
                bins[i1].m_is_valid = 1;
            }
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_x_aabb[i] = bins[i].m_bounding_box;
            debug_x_entries[i] = bins[i].m_num_entries;
            debug_x_exits[i] = bins[i].m_num_exits;
        }
#endif

        // scan bins
        for (int i = 1; i < k_sbvh_num_binning_buckets; i++)
        {
            AABB bound_left;
            int  count_left = 0;
            for (int j = 0; j < i; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_left += bins[j].m_num_entries;
                    bound_left.expandBy(bins[j].m_bounding_box);
                }
            }

            // early skip
            float cost = 0.0f;
            if (count_left > 0)
            {
                cost += bound_left.surfaceAreaCost() * count_left;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost >= m_min_cost)
            {
                continue;
            }

            AABB bound_right;
            int  count_right = 0;
            for (int j = i; j < k_sbvh_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_exits;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceAreaCost() * count_right;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost <= m_min_cost)
            {
                m_min_cost = cost;
                m_min_axis = 0;
                m_min_bin = i;
            }

#if DEBUG_BINS
            if (cost <= debug_x_min_cost)
            {
                debug_x_min_cost = cost;
                debug_x_min_bin = i;
                debug_x_min_left_count = count_left;
                debug_x_min_right_count = count_right;
            }
#endif
        }
    }

    /////////////////// check over y ////////////////////////
    if (edge.y > k_edge_split_threshold)
    {
        SplitBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_y = edge.y / k_sbvh_num_binning_buckets;
        float node_min = m_node_bound.min().y;
        float node_max = m_node_bound.max().y;
        // scan triangles
        for (int i = 0; i < total; i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();

            // entry bucket index
            float y0 = box.min().y;
            assert(y0 >= node_min);
            assert(y0 <= node_max);
            int i0 = i_min(int((y0 - node_min) / delta_y),
                           k_sbvh_num_binning_buckets - 1);
            assert(i0 < k_sbvh_num_binning_buckets);
            assert(i0 >= 0);

            // exit bucket index
            float y1 = box.max().y;
            assert(y1 >= node_min);
            assert(y1 <= node_max);
            int i1 = i_min(int((y1 - node_min) / delta_y),
                           k_sbvh_num_binning_buckets - 1);
            assert(i1 < k_sbvh_num_binning_buckets);
            assert(i1 >= 0);

            AABB  clipped_bbox;
            float bin_min;

            // update entry bin
            bin_min = node_min + i0 * delta_y;
            clipped_bbox =
                m_triangles[i]->clipBoundTight(bin_min, bin_min + delta_y, 1);
            if (clipped_bbox.isValid())
            {
                bins[i0].m_num_entries++;
                bins[i0].m_bounding_box.expandBy(clipped_bbox);
                bins[i0].m_is_valid = 1;
            }
            // update middle bins
            for (int j = i0 + 1; j < i1; j++)
            {
                bin_min = node_min + j * delta_y;
                clipped_bbox = m_triangles[i]->clipBoundTight(
                    bin_min, bin_min + delta_y, 1);
                if (clipped_bbox.isValid())
                {
                    bins[j].m_bounding_box.expandBy(clipped_bbox);
                    bins[j].m_is_valid = 1;
                }
            }
            // update exit bin
            bin_min = node_min + i1 * delta_y;
            clipped_bbox =
                m_triangles[i]->clipBoundTight(bin_min, bin_min + delta_y, 1);
            if (clipped_bbox.isValid())
            {
                bins[i1].m_num_exits++;
                bins[i1].m_bounding_box.expandBy(clipped_bbox);
                bins[i1].m_is_valid = 1;
            }
        }

            ////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_y_aabb[i] = bins[i].m_bounding_box;
            debug_y_entries[i] = bins[i].m_num_entries;
            debug_y_exits[i] = bins[i].m_num_exits;
        }
#endif

        // scan bins
        for (int i = 1; i < k_sbvh_num_binning_buckets; i++)
        {
            AABB bound_left;
            int  count_left = 0;
            for (int j = 0; j < i; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_left += bins[j].m_num_entries;
                    bound_left.expandBy(bins[j].m_bounding_box);
                }
            }

            // early skip
            float cost = 0.0f;
            if (count_left > 0)
            {
                cost += bound_left.surfaceAreaCost() * count_left;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost >= m_min_cost)
            {
                continue;
            }

            AABB bound_right;
            int  count_right = 0;
            for (int j = i; j < k_sbvh_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_exits;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceAreaCost() * count_right;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost <= m_min_cost)
            {
                m_min_cost = cost;
                m_min_axis = 1;
                m_min_bin = i;
            }

#if DEBUG_BINS
            if (cost <= debug_y_min_cost)
            {
                debug_y_min_cost = cost;
                debug_y_min_bin = i;
                debug_y_min_left_count = count_left;
                debug_y_min_right_count = count_right;
            }
#endif
        }
    }

    /////////////////// check over z ////////////////////////
    if (edge.z > k_edge_split_threshold)
    {
        SplitBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_z = edge.z / k_sbvh_num_binning_buckets;
        float node_min = m_node_bound.min().z;
        float node_max = m_node_bound.max().z;
        // scan triangles
        for (int i = 0; i < total; i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();

            // entry bucket index
            float z0 = box.min().z;
            assert(z0 >= node_min);
            assert(z0 <= node_max);
            int i0 = i_min(int((z0 - node_min) / delta_z),
                           k_sbvh_num_binning_buckets - 1);
            assert(i0 < k_sbvh_num_binning_buckets);
            assert(i0 >= 0);

            // exit bucket index
            float z1 = box.max().z;
            assert(z1 >= node_min);
            assert(z1 <= node_max);
            int i1 = i_min(int((z1 - node_min) / delta_z),
                           k_sbvh_num_binning_buckets - 1);
            assert(i1 < k_sbvh_num_binning_buckets);
            assert(i1 >= 0);

            AABB  clipped_bbox;
            float bin_min;

            // update entry bin
            bin_min = node_min + i0 * delta_z;
            clipped_bbox =
                m_triangles[i]->clipBoundTight(bin_min, bin_min + delta_z, 2);
            if (clipped_bbox.isValid())
            {
                bins[i0].m_num_entries++;
                bins[i0].m_bounding_box.expandBy(clipped_bbox);
                bins[i0].m_is_valid = 1;
            }
            // update middle bins
            for (int j = i0; j < i1; j++)
            {
                bin_min = node_min + j * delta_z;
                clipped_bbox = m_triangles[i]->clipBoundTight(
                    bin_min, bin_min + delta_z, 2);
                if (clipped_bbox.isValid())
                {
                    bins[j].m_bounding_box.expandBy(clipped_bbox);
                    bins[j].m_is_valid = 1;
                }
            }
            // update exit bin
            bin_min = node_min + i1 * delta_z;
            clipped_bbox =
                m_triangles[i]->clipBoundTight(bin_min, bin_min + delta_z, 2);
            if (clipped_bbox.isValid())
            {
                bins[i1].m_num_exits++;
                bins[i1].m_bounding_box.expandBy(clipped_bbox);
                bins[i1].m_is_valid = 1;
            }
        }

            ////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_z_aabb[i] = bins[i].m_bounding_box;
            debug_z_entries[i] = bins[i].m_num_entries;
            debug_z_exits[i] = bins[i].m_num_exits;
        }
#endif

        // scan bins
        for (int i = 1; i < k_sbvh_num_binning_buckets; i++)
        {
            AABB bound_left;
            int  count_left = 0;
            for (int j = 0; j < i; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_left += bins[j].m_num_entries;
                    bound_left.expandBy(bins[j].m_bounding_box);
                }
            }

            // early skip
            float cost = 0.0f;
            if (count_left > 0)
            {
                cost += bound_left.surfaceAreaCost() * count_left;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost >= m_min_cost)
            {
                continue;
            }

            AABB bound_right;
            int  count_right = 0;
            for (int j = i; j < k_sbvh_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_exits;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceAreaCost() * count_right;
            }
            else
            {
                cost = NUM_INFINITY;
            }
            if (cost <= m_min_cost)
            {
                m_min_cost = cost;
                m_min_axis = 2;
                m_min_bin = i;
            }

#if DEBUG_BINS
            if (cost <= debug_z_min_cost)
            {
                debug_z_min_cost = cost;
                debug_z_min_bin = i;
                debug_z_min_left_count = count_left;
                debug_z_min_right_count = count_right;
            }
#endif
        }
    }

        ////////////// debug ////////////////
#if DEBUG_BINS
    Logger::info() << "================================================"
                   << std::endl;
#endif

#if DEBUG_INFO
    if (m_min_axis >= 0)
    {
        char* axes[3] = {"x", "y", "z"};
        Logger::info() << "<split bvh> num primitives = " << m_triangles.size()
                       << ", min cost = " << m_min_cost
                       << ", min axis = " << axes[m_min_axis]
                       << ", min bin = " << m_min_bin << std::endl;
    }
    else
    {
        Logger::info() << "<split bvh> num primitives = " << total
                       << ", cannot find min with edge split threshold "
                       << k_edge_split_threshold << std::endl;
    }
#endif

        ////////////// debug ////////////////
#if DEBUG_BINS
    // if (m_min_bin <= 2 || m_min_bin >= k_sbvh_num_binning_buckets - 2)
    {
        // x
        Logger::info() << ">> x bins min cost = " << debug_x_min_cost
                       << ", left count = " << debug_x_min_left_count
                       << ", right count = " << debug_x_min_right_count
                       << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_x_entries[i] << "/"
                           << debug_x_exits[i] << ", "
                           << debug_x_aabb[i].toString() << std::endl;
        }
        // y
        Logger::info() << ">> y bins min cost = " << debug_y_min_cost
                       << ", left count = " << debug_y_min_left_count
                       << ", right count = " << debug_y_min_right_count
                       << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_y_entries[i] << "/"
                           << debug_y_exits[i] << ", "
                           << debug_y_aabb[i].toString() << std::endl;
        }
        // z
        Logger::info() << ">> z bins min cost = " << debug_z_min_cost
                       << ", left count = " << debug_z_min_left_count
                       << ", right count = " << debug_z_min_right_count
                       << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_z_entries[i] << "/"
                           << debug_z_exits[i] << ", "
                           << debug_z_aabb[i].toString() << std::endl;
        }
    }
#endif
}

float SplitBVHSplitter::cost() const
{
    return m_min_cost;
}

void SplitBVHSplitter::sort(
    int&                                            axis,
    std::vector<std::shared_ptr<BoundedTriangle> >& left_branch,
    std::vector<std::shared_ptr<BoundedTriangle> >& right_branch)
{
    axis = m_min_axis;

    SplitBVHBinningBucket bins[k_sbvh_num_binning_buckets];

    float delta = m_edge[m_min_axis] / k_sbvh_num_binning_buckets;
    float node_min_v[3] = {
        m_node_bound.min().x, m_node_bound.min().y, m_node_bound.min().z};
    float node_max_v[3] = {
        m_node_bound.max().x, m_node_bound.max().y, m_node_bound.max().z};
    float node_min = node_min_v[m_min_axis];
    float node_max = node_max_v[m_min_axis];
    float node_clip = node_min + m_min_bin * delta;
    assert(node_clip > node_min);
    assert(node_clip < node_max);

    int count = (int)m_triangles.size();
    left_branch.clear();
    right_branch.clear();

    // scan triangles
    for (int i = 0; i < count; i++)
    {
        const auto& box = m_triangles[i]->clippedBoundingBox();
        float       min_v[3] = {box.min().x, box.min().y, box.min().z};
        float       max_v[3] = {box.max().x, box.max().y, box.max().z};

        // entry bucket index
        float v0 = min_v[m_min_axis];
        assert(v0 >= node_min);
        assert(v0 <= node_max);
        int i0 =
            i_min(int((v0 - node_min) / delta), k_sbvh_num_binning_buckets - 1);
        assert(i0 < k_sbvh_num_binning_buckets);
        assert(i0 >= 0);

        // exit bucket index
        float v1 = max_v[m_min_axis];
        assert(v1 >= node_min);
        assert(v1 <= node_max);
        int i1 =
            i_min(int((v1 - node_min) / delta), k_sbvh_num_binning_buckets - 1);
        assert(i1 < k_sbvh_num_binning_buckets);
        assert(i1 >= 0);

        // Logger::info() << "tri " << i << " : " << i0 << " / " << i1
        //               << std::endl;

        // push to the left branch only for entries
        if (i0 < m_min_bin)
        {
            AABB clipped_bbox =
                m_triangles[i]->clipBoundTight(node_min, node_clip, m_min_axis);
            if (clipped_bbox.isValid())
            {
                left_branch.push_back(std::make_shared<BoundedTriangle>(
                    m_triangles[i], clipped_bbox));
            }
        }
        // push to the right branch only for exits
        if (i1 >= m_min_bin)
        {
            AABB clipped_bbox =
                m_triangles[i]->clipBoundTight(node_clip, node_max, m_min_axis);
            if (clipped_bbox.isValid())
            {
                right_branch.push_back(std::make_shared<BoundedTriangle>(
                    m_triangles[i], clipped_bbox));
            }
        }
    }

        ////////////////////////// debug output
        ////////////////////////// debug output
        ////////////////////////// debug output
        ////////////////////////// debug output
        ////////////////////////// debug output
#if DEBUG_BINS
    {
        float3 edge(m_edge[0], m_edge[1], m_edge[2]);

        AABB             debug_z_aabb[k_sbvh_num_binning_buckets];
        int              debug_z_entries[k_sbvh_num_binning_buckets];
        std::vector<int> debug_z_entries_tri[k_sbvh_num_binning_buckets];
        int              debug_z_exits[k_sbvh_num_binning_buckets];
        std::vector<int> debug_z_exits_tri[k_sbvh_num_binning_buckets];
        int              debug_z_valids[k_sbvh_num_binning_buckets];

        SplitBVHBinningBucket bins[k_sbvh_num_binning_buckets];

        float delta_z = edge.z / k_sbvh_num_binning_buckets;
        float node_min = m_node_bound.min().z;
        float node_max = m_node_bound.max().z;
        // scan triangles
        for (int i = 0; i < m_triangles.size(); i++)
        {
            const auto& box = m_triangles[i]->clippedBoundingBox();

            // entry bucket index
            float z0 = box.min().z;
            assert(z0 >= node_min);
            assert(z0 <= node_max);
            int i0 = std::min(int((z0 - node_min) / delta_z),
                              k_sbvh_num_binning_buckets - 1);
            assert(i0 < k_sbvh_num_binning_buckets);
            assert(i0 >= 0);

            // exit bucket index
            float z1 = box.max().z;
            assert(z1 >= node_min);
            assert(z1 <= node_max);
            int i1 = std::min(int((z1 - node_min) / delta_z),
                              k_sbvh_num_binning_buckets - 1);
            assert(i1 < k_sbvh_num_binning_buckets);
            assert(i1 >= 0);

            AABB  clipped_bbox;
            float bin_min;

            // update entry bin
            bins[i0].m_num_entries++;
            bin_min = node_min + i0 * delta_z;
            clipped_bbox = box.clipZ(bin_min, bin_min + delta_z);
            bins[i0].m_bounding_box.expandBy(clipped_bbox);
            debug_z_entries_tri[i0].push_back(i);
            debug_z_valids[i0] = 1;
            // update middle bins
            for (int j = i0; j < i1; j++)
            {
                bin_min = node_min + j * delta_z;
                clipped_bbox = box.clipZ(bin_min, bin_min + delta_z);
                bins[j].m_bounding_box.expandBy(clipped_bbox);
                debug_z_valids[j] = 1;
            }
            // update exit bin
            bins[i1].m_num_exits++;
            bin_min = node_min + i1 * delta_z;
            clipped_bbox = box.clipZ(bin_min, bin_min + delta_z);
            bins[i1].m_bounding_box.expandBy(clipped_bbox);
            debug_z_exits_tri[i1].push_back(i);
            debug_z_valids[i1] = 1;
        }

            ////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            debug_z_aabb[i] = bins[i].m_bounding_box;
            debug_z_entries[i] = bins[i].m_num_entries;
            debug_z_exits[i] = bins[i].m_num_exits;
        }
#endif
        Logger::info() << "---------------------------" << std::endl;
        for (int i = 0; i < k_sbvh_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_z_entries[i] << "/"
                           << debug_z_exits[i] << "[" << debug_z_valids[i]
                           << "]"
                           << " : ";
            for (auto j : debug_z_entries_tri[i])
            {
                Logger::info() << j << ", ";
            }
            Logger::info() << "|";
            for (auto j : debug_z_exits_tri[i])
            {
                Logger::info() << j << ", ";
            }

            Logger::info() << std::endl;
        }
        Logger::info() << "---------------------------" << std::endl;
    }
#endif
    ////////////////////////// debug output
    ////////////////////////// debug output
    ////////////////////////// debug output
    ////////////////////////// debug output
    ////////////////////////// debug output

    assert(false == left_branch.empty());
    assert(false == right_branch.empty());

    // if (left_branch.size() == 6 && right_branch.size() == 14)
    //{
    //    Logger::info() << "break";
    // Logger::info() << std::endl;
    //}
}

SplitBVH::SplitBVH(const std::vector<Triangle>& mesh)
{
    m_triangles.clear();
    for (const auto& t : mesh)
    {
        // copy constructor
        m_triangles.push_back(std::make_shared<Triangle>(t));
    }

    m_triangle_refs.clear();
    for (const auto& t : m_triangles)
    {
        m_triangle_refs.push_back(
            std::make_shared<BoundedTriangle>(t, t->boundingBox()));
    }

    m_rootNode = build(m_triangle_refs, 0u, "root");
}

std::shared_ptr<debug::SplitBVHNode> SplitBVH::build(
    const std::vector<std::shared_ptr<BoundedTriangle> >& triangles,
    uint32_t                                              level,
    const std::string&                                    direction)
{
    assert(false == triangles.empty());

#if DEBUG_BINS
    if (triangles.size() < 16)
    {
        Logger::info() << "# tri = " << triangles.size() << std::endl;

        for (auto t1 : triangles)
        {
            auto t = t1->originalTriangle();

            auto a = t->a;
            auto b = t->b;
            auto c = t->c;

            Logger::info() << a.x << ", " << a.y << ", " << a.z << std::endl;
            Logger::info() << b.x << ", " << b.y << ", " << b.z << std::endl;
            Logger::info() << c.x << ", " << c.y << ", " << c.z << std::endl;
        }
    }
#endif

    if (triangles.size() <= k_sbvh_leaf_node_size)  // leaf node
    {
        auto node = std::make_shared<SplitBVHLeafNode>();
        auto bounding_box = std::make_shared<AABB>();

        node->debug_level = level;
        node->debug_tag = std::string("leaf ") + direction;

        int count = (int)triangles.size();
        for (auto i = 0; i < count; i++)
        {
            auto tri = triangles[i];

            node->primitive[i] = tri;

            bounding_box->expandBy(tri->clippedBoundingBox());
        }
        node->bounding_box = bounding_box;

        return node;
    }
    else
    {
        auto node = std::make_shared<SplitBVHBranchNode>();

        auto bounding_box = std::make_shared<AABB>();
        int  count = (int)triangles.size();
        for (int i = 0; i < count; ++i)
        {
            const auto& box = triangles[i]->clippedBoundingBox();
            bounding_box->expandBy(box);
        }

        std::vector<std::shared_ptr<BoundedTriangle> > left_branch;
        std::vector<std::shared_ptr<BoundedTriangle> > right_branch;
        int                                            split_axis;

        // debug disable splitting
#if ENABLE_MEDIAN_OBJECT_SPLIT
        MedianSplitBVHSplitter median_split(triangles, *bounding_box);
        float                  median_split_cost = median_split.cost();
#else
        float median_split_cost = NUM_INFINITY;
#endif
#if ENABLE_SAH_OBJECT_SPLIT
        SAHBVHSplitter sah_split(triangles, *bounding_box);
        float          sah_split_cost = sah_split.cost();
#else
        float sah_split_cost = NUM_INFINITY;
#endif
#if ENABLE_SAH_SPATIAL_SPLIT
        SplitBVHSplitter spatial_split(triangles, *bounding_box);
        float            spatial_split_cost = spatial_split.cost();
#else
        float spatial_split_cost = NUM_INFINITY;
#endif

        // sort low to high cost, and prioritize sah split
        int   types[3] = {1, 0, 2};
        float costs[3] = {
            sah_split_cost, median_split_cost, spatial_split_cost};
        if (costs[0] > costs[1])
        {
            std::swap(costs[0], costs[1]);
            std::swap(types[0], types[1]);
        }
        if (costs[1] > costs[2])
        {
            std::swap(costs[1], costs[2]);
            std::swap(types[1], types[2]);
        }
        if (costs[0] > costs[1])
        {
            std::swap(costs[0], costs[1]);
            std::swap(types[0], types[1]);
        }

        assert(costs[0] < NUM_INFINITY);

        switch (types[0])
        {
            case 0:
#if DEBUG_INFO
                Logger::info()
                    << "00000000 median split is better : " << median_split_cost
                    << ", " << sah_split_cost << ", " << spatial_split_cost
                    << std::endl;
#endif

#if ENABLE_MEDIAN_OBJECT_SPLIT
                median_split.sort(split_axis, left_branch, right_branch);
#else
                Logger::info()
                    << "invalid : median split is disabled" << std::endl;
                exit(-1);
#endif

                break;
            case 1:
#if DEBUG_INFO
                Logger::info()
                    << "11111111 sah split is better : " << median_split_cost
                    << ", " << sah_split_cost << ", " << spatial_split_cost
                    << std::endl;
#endif

#if ENABLE_SAH_OBJECT_SPLIT
                sah_split.sort(split_axis, left_branch, right_branch);
#else
                Logger::info()
                    << "invalid : sah object split is disabled" << std::endl;
                exit(-1);
#endif

                break;
            case 2:
                // otherwise use split bvh (may increase reference)
#if DEBUG_INFO
                Logger::info() << "22222222 spatial split is better : "
                               << median_split_cost << ", " << sah_split_cost
                               << ", " << spatial_split_cost << std::endl;
#endif

#if ENABLE_SAH_SPATIAL_SPLIT
                spatial_split.sort(split_axis, left_branch, right_branch);
#else
                Logger::info()
                    << "invalid : sah spatial split is disabled" << std::endl;
                exit(-1);
#endif

                break;
            default:
                Logger::info() << "invalid splitting technique" << std::endl;
                exit(-1);
                break;
        }

#if DEBUG_INFO
        Logger::info() << ">> " << triangles.size() << " -> "
                       << left_branch.size() << " | " << right_branch.size()
                       << " ("
                       << left_branch.size() + right_branch.size() -
                              triangles.size()
                       << " duplicates)" << std::endl;
#endif

        node->left_child = build(left_branch, level + 1, "left");
        node->right_child = build(right_branch, level + 1, "right");
        node->split_axis = split_axis;

        node->debug_level = level;
        node->debug_tag = std::string("branch ") + direction;

#if 0  // debug rebuild bounding box
	            bounding_box = std::make_shared<AABB>();
	            bounding_box->expandBy(*node->left_child->boundingBox);
	            bounding_box->expandBy(*node->right_child->boundingBox);
#endif
        node->bounding_box = bounding_box;

        return node;
    }
}

bool SplitBVH::intersect(Ray& ray, HitInfo& hit_info) const
{
    // BVH traversal
    hit_info.reset();

    bool hit = m_rootNode->intersect(ray, hit_info);

    if (hit)
    {
        return true;
    }
    else
    {
        return false;
    }
}
}