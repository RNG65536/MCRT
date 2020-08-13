#include <algorithm>
#include "sahbvh.h"

bool SAHBVHLeafNode::intersect(Ray& ray, HitInfo& hit_info) const
{
    hit_info.reset();

    bool hitBox = bounding_box->overlap(ray);
    if (!hitBox)
    {
        return false;
    }

    bool    any_hit = false;
    HitInfo hit_info_child;

    for (int i = 0; i < k_leaf_node_size; i++)
    {
        if (nullptr == primitive[i])
        {
            continue;
        }

        bool hit = primitive[i]->intersect(ray, hit_info_child);
        if (hit)
        {
            // Logger::debug() << "[ " << debug_level << " ][ " << debug_tag
            //                << " ]" << std::endl;
            // Logger::debug()
            //    << "    primitive t_near = " << t_near << std::endl;

            assert(hit_info_child.m_t > ray.tmin &&
                   hit_info_child.m_t < ray.tmax);

            if (hit_info_child < hit_info)
            {
                hit_info = hit_info_child;
                any_hit = true;
            }

            if (hit_info.m_t < ray.tmax)
            {
                ray.tmax = hit_info.m_t;
            }
        }
    }
    return any_hit;  // if no hit, returns infinity by default
}

bool SAHBVHBranchNode::intersect(Ray& ray, HitInfo& hit_info) const
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
                assert(hit_info.m_t <= ray.tmax);
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
                assert(hit_info.m_t <= ray.tmax);
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
                assert(hit_info.m_t <= ray.tmax);
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
                assert(hit_info.m_t <= ray.tmax);
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

SAHSplitAxisComparator::SAHSplitAxisComparator(int split_axis)
    : m_split_axis(split_axis)
{
}

bool SAHSplitAxisComparator::operator()(
    const std::shared_ptr<Triangle>& a,
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

MedianSplitter::MedianSplitter(
    std::vector<std::shared_ptr<Triangle> >& triangles,
    const AABB&                              node_bound,
    int                                      begin,
    int                                      end)
    : m_triangles(triangles)
    , m_node_bound(node_bound)
    , m_begin(begin)
    , m_end(end)
{
    AABB left_bound;
    AABB right_bound;
    m_split_axis = m_node_bound.largestDimension();
    std::sort(m_triangles.begin() + m_begin,
              m_triangles.begin() + m_end,
              SAHSplitAxisComparator(m_split_axis));
    m_mid = m_begin + (m_end - m_begin) / 2;
    int total = m_end - m_begin;
    for (int i = m_begin; i < m_mid; i++)
    {
        left_bound.expandBy(m_triangles[i]->boundingBox());
    }
    for (int i = m_mid; i < m_end; i++)
    {
        right_bound.expandBy(m_triangles[i]->boundingBox());
    }
    int left_count = m_mid - m_begin;
    int right_count = total - left_count;
    m_cost = left_bound.surfaceArea() * left_count +
             right_bound.surfaceArea() * right_count;

#if DEBUG_INFO
    char* axes[3] = {"x", "y", "z"};
    Logger::info() << "<median split bvh> num primitives = " << total
                   << ", min cost = " << m_cost
                   << ", min axis = " << axes[m_split_axis] << std::endl;
#endif
}

float MedianSplitter::cost() const
{
    return m_cost;
}

void MedianSplitter::sort(int& axis, int& mid)
{
    axis = m_split_axis;
    mid = m_mid;
}

SAHSplitter::SAHSplitter(std::vector<std::shared_ptr<Triangle> >& triangles,
                         const AABB&                              node_bound,
                         int                                      begin,
                         int                                      end)
    : m_triangles(triangles)
    , m_node_bound(node_bound)
    , m_begin(begin)
    , m_end(end)
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

    ////////////// debug ////////////////
#if DEBUG_BINS
    AABB  debug_x_aabb[k_num_binning_buckets];
    AABB  debug_y_aabb[k_num_binning_buckets];
    AABB  debug_z_aabb[k_num_binning_buckets];
    int   debug_x_count[k_num_binning_buckets];
    int   debug_y_count[k_num_binning_buckets];
    int   debug_z_count[k_num_binning_buckets];
    float debug_x_min_cost = NUM_INFINITY;
    float debug_y_min_cost = NUM_INFINITY;
    float debug_z_min_cost = NUM_INFINITY;
#endif

    /////////////////// check over x ////////////////////////
    if (edge.x > k_edge_threshold)
    {
        SAHBinningBucket bins[k_num_binning_buckets];

        float delta_x = edge.x / k_num_binning_buckets;
        // scan triangles
        for (int i = m_begin; i < m_end; i++)
        {
            const auto& box = m_triangles[i]->boundingBox();
            float       x0 = box.min().x;
            assert(x0 >= node_bound.min().x);
            assert(x0 <= node_bound.max().x);
            int bucket_index = i_min(int((x0 - node_bound.min().x) / delta_x),
                                     k_num_binning_buckets - 1);

            assert(bucket_index < k_num_binning_buckets);
            assert(bucket_index >= 0);

            bins[bucket_index].m_num_entries++;
            bins[bucket_index].m_bounding_box.expandBy(box);
            bins[bucket_index].m_is_valid = 1;
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_num_binning_buckets; i++)
        {
            debug_x_aabb[i] = bins[i].m_bounding_box;
            debug_x_count[i] = bins[i].m_num_entries;
        }
#endif

        // scan bins
        for (int i = 1; i < k_num_binning_buckets; i++)
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
                cost += bound_left.surfaceArea() * count_left;
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
            for (int j = i; j < k_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_entries;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceArea() * count_right;
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
    if (edge.y > k_edge_threshold)
    {
        SAHBinningBucket bins[k_num_binning_buckets];

        float delta_y = edge.y / k_num_binning_buckets;
        // scan triangles
        for (int i = m_begin; i < m_end; i++)
        {
            const auto& box = m_triangles[i]->boundingBox();
            float       y0 = box.min().y;
            assert(y0 >= node_bound.min().y);
            assert(y0 <= node_bound.max().y);
            int bucket_index = i_min(int((y0 - node_bound.min().y) / delta_y),
                                     k_num_binning_buckets - 1);

            assert(bucket_index < k_num_binning_buckets);
            assert(bucket_index >= 0);

            bins[bucket_index].m_num_entries++;
            bins[bucket_index].m_bounding_box.expandBy(box);
            bins[bucket_index].m_is_valid = 1;
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_num_binning_buckets; i++)
        {
            debug_y_aabb[i] = bins[i].m_bounding_box;
            debug_y_count[i] = bins[i].m_num_entries;
        }
#endif

        // scan bins
        for (int i = 1; i < k_num_binning_buckets; i++)
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
                cost += bound_left.surfaceArea() * count_left;
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
            for (int j = i; j < k_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_entries;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceArea() * count_right;
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
    if (edge.z > k_edge_threshold)
    {
        SAHBinningBucket bins[k_num_binning_buckets];

        float delta_z = edge.z / k_num_binning_buckets;
        // scan triangles
        for (int i = m_begin; i < m_end; i++)
        {
            const auto& box = m_triangles[i]->boundingBox();
            float       z0 = box.min().z;
            assert(z0 >= node_bound.min().z);
            assert(z0 <= node_bound.max().z);
            int bucket_index = i_min(int((z0 - node_bound.min().z) / delta_z),
                                     k_num_binning_buckets - 1);

            assert(bucket_index < k_num_binning_buckets);
            assert(bucket_index >= 0);

            bins[bucket_index].m_num_entries++;
            bins[bucket_index].m_bounding_box.expandBy(box);
            bins[bucket_index].m_is_valid = 1;
        }

////////////// debug ////////////////
#if DEBUG_BINS
        for (int i = 0; i < k_num_binning_buckets; i++)
        {
            debug_z_aabb[i] = bins[i].m_bounding_box;
            debug_z_count[i] = bins[i].m_num_entries;
        }
#endif

        // scan bins
        for (int i = 1; i < k_num_binning_buckets; i++)
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
                cost += bound_left.surfaceArea() * count_left;
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
            for (int j = i; j < k_num_binning_buckets; j++)
            {
                if (bins[j].m_is_valid)
                {
                    count_right += bins[j].m_num_entries;
                    bound_right.expandBy(bins[j].m_bounding_box);
                }
            }

            if (count_right > 0)
            {
                cost += bound_right.surfaceArea() * count_right;
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
    Logger::info() << "num primitives = " << m_end - m_begin
                   << ", min cost = " << m_min_cost
                   << ", min axis = " << m_min_axis
                   << ", min bin = " << m_min_bin << std::endl;
#endif

////////////// debug ////////////////
#if DEBUG_BINS
    // if (m_min_bin == 1 || m_min_bin == k_num_binning_buckets - 1)
    {
        // x
        Logger::info() << ">> x bins min cost " << debug_x_min_cost
                       << std::endl;
        for (int i = 0; i < k_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_x_count[i] << ", "
                           << debug_x_aabb[i].toString() << std::endl;
        }
        // y
        Logger::info() << ">> y bins min cost " << debug_y_min_cost
                       << std::endl;
        for (int i = 0; i < k_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_y_count[i] << ", "
                           << debug_y_aabb[i].toString() << std::endl;
        }
        // z
        Logger::info() << ">> z bins min cost " << debug_z_min_cost
                       << std::endl;
        for (int i = 0; i < k_num_binning_buckets; i++)
        {
            Logger::info() << i << " : " << debug_z_count[i] << ", "
                           << debug_z_aabb[i].toString() << std::endl;
        }
    }
#endif
}

void SAHSplitter::sort(int& axis, int& mid)
{
    axis = m_min_axis;

    SAHBinningBucket bins[k_num_binning_buckets];

    float delta = m_edge[m_min_axis] / k_num_binning_buckets;
    float node_min_v[3] = {
        m_node_bound.min().x, m_node_bound.min().y, m_node_bound.min().z};

    std::vector<std::shared_ptr<Triangle> > temp(m_end - m_begin);
    int                                     begin_offset = 0;
    int                                     end_offset = m_end - m_begin - 1;

    // scan triangles
    for (int i = m_begin; i < m_end; i++)
    {
        const auto& box = m_triangles[i]->boundingBox();
        float       min_v[3] = {box.min().x, box.min().y, box.min().z};
        int         bucket_index =
            i_min(int((min_v[m_min_axis] - node_min_v[m_min_axis]) / delta),
                  k_num_binning_buckets - 1);

        assert(bucket_index < k_num_binning_buckets);
        assert(bucket_index >= 0);

        assert(begin_offset <= end_offset);

        // push to either end of the temp array
        // TODO : replace this with in-place swapping
        if (bucket_index < m_min_bin)
        {
            temp[begin_offset++] = m_triangles[i];
        }
        else
        {
            temp[end_offset--] = m_triangles[i];
        }
    }

    // write back the sorted result
    for (int i = 0; i < temp.size(); i++)
    {
        m_triangles[m_begin + i] = temp[i];
    }

    assert(begin_offset > 0);
    assert(begin_offset < m_end - m_begin);
    mid = m_begin + begin_offset;
}

SAHBVH::SAHBVH(const std::vector<Triangle>& mesh)
{
    m_triangles.clear();
    for (const auto& t : mesh)
    {
        // copy constructor
        m_triangles.push_back(std::make_shared<Triangle>(t));
    }

    m_rootNode = build(0, (int)m_triangles.size(), 0u, "root");
}

std::shared_ptr<SAHBVHNode> SAHBVH::buildLeaf(
    int begin, int end, uint32_t level, const std::string& direction)
{
    auto node = std::make_shared<SAHBVHLeafNode>();
    auto bounding_box = std::make_shared<AABB>();

    node->debug_level = level;
    node->debug_tag = std::string("leaf ") + direction;

    int count = end - begin;
    for (auto i = 0; i < count; i++)
    {
        assert(begin + i < m_triangles.size());
        auto tri = m_triangles[begin + i];

        node->primitive[i] = tri;

        bounding_box->expandBy(tri->boundingBox());
    }
    node->bounding_box = bounding_box;

    return node;
}

std::shared_ptr<SAHBVHBranchNode> SAHBVH::buildBranch(
    int begin, int end, uint32_t level, const std::string& direction, int& mid)
{
    auto node = std::make_shared<SAHBVHBranchNode>();

    auto bounding_box = std::make_shared<AABB>();
    for (int i = begin; i < end; ++i)
    {
        const auto& box = m_triangles[i]->boundingBox();
        bounding_box->expandBy(box);
    }

#if 0  // debug rebuild bounding box
	            bounding_box = std::make_shared<AABB>();
	            bounding_box->expandBy(*node->left_child->boundingBox);
	            bounding_box->expandBy(*node->right_child->boundingBox);
#endif
    node->bounding_box = bounding_box;

    node->debug_level = level;
    node->debug_tag = std::string("branch ") + direction;

    int         split_axis;
    SAHSplitter sah(m_triangles, *bounding_box, begin, end);

    float cost = sah.m_min_cost;
    if (cost == NUM_INFINITY)
    {
#if DEBUG_INFO
        Logger::info() << "resort to median split" << std::endl;
#endif
        MedianSplitter median_split(m_triangles, *bounding_box, begin, end);
        median_split.sort(split_axis, mid);
    }
    else
    {
        sah.sort(split_axis, mid);
    }
    node->split_axis = split_axis;

    return node;
}

std::shared_ptr<SAHBVHNode> SAHBVH::build(int                begin,
                                                 int                end,
                                                 uint32_t           level,
                                                 const std::string& direction)
{
    assert(begin < end);

    if (end - begin <= k_leaf_node_size)  // leaf node
    {
        return buildLeaf(begin, end, level, direction);
    }
    else
    {
        int  mid = -1;
        auto node = buildBranch(begin, end, level, direction, mid);
        node->left_child = build(begin, mid, level + 1, "left");
        node->right_child = build(mid, end, level + 1, "right");

        return node;
    }
}

bool SAHBVH::intersect(Ray& ray, HitInfo& hit_info) const
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
