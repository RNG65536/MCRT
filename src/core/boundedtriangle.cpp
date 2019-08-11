#include "boundedtriangle.h"

AABB BoundedTriangle::clipBoundTight(float min,
                                                   float max,
                                                   int   axis) const
{
#if 0
	        const auto& clipped_bound = m_clipped_bound;
	        auto more_clipped_bound = clipped_bound.clipAxis(min, max, axis);
	        return more_clipped_bound;
#else
    // because the number of triangles is often larger than the number of
    // bins it is better to compute a tight bound for each bin, rather than
    // scanning the triangles and computing a loose bound for each side when
    // shifting the split plane
    const auto& clipped_bound = m_clipped_bound;
    const auto& full_bound = m_triangle->boundingBox();

    if (clipped_bound.min()[axis] > max || clipped_bound.max()[axis] < min)
    {
        return AABB();
    }
    else
    {
        const float3& a = m_triangle->a;
        const float3& b = m_triangle->b;
        const float3& c = m_triangle->c;

        // compute intersection of the three edges of the triangle with the
        // clipping planes (no intersection, or one point, or two points)
        // the tight bound of the portion of the triangle in the slab
        // is the bound on the vertices in the slab and the intersection
        // points.
        // FIXED : use inclusive borders (>= and <=) to avoid ray leaking
        // when testing.
        AABB tight_bound;

        // solve a + t * (b - a) = pmin and a + t * (b - a) = pmax if ab
        // intersects min or max plane
        {
            float3 ba = b - a;
            float  tmin = (min - a[axis]) / ba[axis];
            float  tmax = (max - a[axis]) / ba[axis];
            if (tmin >= 0.0f && tmin <= 1.0f)
            {
                tight_bound.expandBy(a + ba * tmin);
            }
            if (tmax >= 0.0f && tmax <= 1.0f)
            {
                tight_bound.expandBy(a + ba * tmax);
            }
        }
        // solve b + t * (c - b) = pmin and b + t * (c - b) = pmax if bc
        // intersects min or max plane
        {
            float3 cb = c - b;
            float  tmin = (min - b[axis]) / cb[axis];
            float  tmax = (max - b[axis]) / cb[axis];
            if (tmin >= 0.0f && tmin <= 1.0f)
            {
                tight_bound.expandBy(b + cb * tmin);
            }
            if (tmax >= 0.0f && tmax <= 1.0f)
            {
                tight_bound.expandBy(b + cb * tmax);
            }
        }
        // solve c + t * (a - c) = pmin and c + t * (a - c) = pmax if ca
        // intersects min or max plane
        {
            float3 ac = a - c;
            float  tmin = (min - c[axis]) / ac[axis];
            float  tmax = (max - c[axis]) / ac[axis];
            if (tmin >= 0.0f && tmin <= 1.0f)
            {
                tight_bound.expandBy(c + ac * tmin);
            }
            if (tmax >= 0.0f && tmax <= 1.0f)
            {
                tight_bound.expandBy(c + ac * tmax);
            }
        }
        if (a[axis] >= min && a[axis] <= max)
        {
            tight_bound.expandBy(a);
        }
        if (b[axis] >= min && b[axis] <= max)
        {
            tight_bound.expandBy(b);
        }
        if (c[axis] >= min && c[axis] <= max)
        {
            tight_bound.expandBy(c);
        }

        auto more_clipped_bound = clipped_bound.clipAxis(min, max, axis);

#if 1
        // compute intersection of the tight bound and the triangle's
        // clipped bound
        tight_bound.intersectBy(more_clipped_bound);
        if (tight_bound.isValid())
        {
            return more_clipped_bound.intersect(tight_bound);
        }
#else
        // this is a fallback solution for reference
        tight_bound = full_bound.clipAxis(min, max, axis);
        tight_bound.intersectBy(more_clipped_bound);
        if (tight_bound.isValid())
        {
            return more_clipped_bound.intersect(tight_bound);
        }
#endif

        if (more_clipped_bound.isValid())
        {
            return more_clipped_bound;
        }

        return AABB();
    }
#endif
}

bool BoundedTriangle::intersectBound(Ray& ray) const
{
    return m_clipped_bound.overlap(ray);
}

bool BoundedTriangle::intersect(Ray& ray, HitInfo& hit_info) const
{
    if (false == m_clipped_bound.overlap(ray))
    {
        hit_info.reset();
        return false;
    }
    return m_triangle->intersect(ray, hit_info);
}

const std::shared_ptr<Triangle>& BoundedTriangle::originalTriangle()
    const
{
    return m_triangle;
}

const AABB& BoundedTriangle::clippedBoundingBox() const
{
    return m_clipped_bound;
}

BoundedTriangle::BoundedTriangle(
    const std::shared_ptr<BoundedTriangle>& triangle, const AABB& aabb)
    : m_triangle(triangle->m_triangle), m_clipped_bound(aabb)
{
}

BoundedTriangle::BoundedTriangle(
    const std::shared_ptr<Triangle>& triangle, const AABB& aabb)
    : m_triangle(triangle), m_clipped_bound(aabb)
{
}
