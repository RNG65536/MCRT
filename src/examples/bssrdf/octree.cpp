#include "octree.h"
#include "numeric.h"

vec3 OctreeSample::getRadiantExitance(vec3 position, vec3 normal, vec3 direction /* view ray */,
    std::function<vec3(const vec3&, const vec3&, const vec3&, const vec3&, const vec3&, const vec3&)> kernel)
{
    vec3& sample_loc = centroid;
    //         real r = (sample_loc - position).length();

    //         return sumPower[channel] * Rd(r, channel);
    //         return sumPower[channel] MULTIPLY_BY_FOUR_PI * bssrdf(sample_loc, 

    //         const real T12 = 1.0 - fresnel(f_max(0, dot(weightedNormal, weightedLightDir)), ETA);
    //         return mult(mult((F_dt MULTIPLY_BY_FOUR_PI), sumPower), bssrdf(sample_loc, 
    //         return mult(mult(F_dt, sumPower), bssrdf(sample_loc, 
    return meanIrradiance * sumArea * kernel(sample_loc, weightedNormal, weightedLightDir, position, normal, direction);
}

OctreeSample::OctreeSample(const vec3& position, const vec3& normal, const vec3& direction, const vec3& rad, float area) : centroid(position), meanIrradiance(rad), weightedNormal(normal), weightedLightDir(direction), sumArea(area)
{

}

float Octree::meanIrradianceRecursive() const
{
    if (isLeafNode())
    {
        return luminance(this->payload->meanIrradiance);
    }
    else
    {
        float sum_power = 0;
        float sum_area = 0;
        for (int i = 0; i < 8; ++i)
        {
            if (children[i]->payload)
            {
                sum_power += children[i]->payload->sumArea * children[i]->meanIrradianceRecursive();
                sum_area += children[i]->payload->sumArea;
            }
        }
        if (sum_area > 0)
        {
            sum_power /= sum_area;
        }
        return sum_power;
    }
}

void Octree::finalize()
{
    this->buildAggregationRecursive();

//     float mean_irr = luminance(this->payload->meanIrradiance);
//     float mean_irr2 = this->meanIrradianceRecursive();
//     cout << "!!!!!!!!!!!!! " << mean_irr << ", " << mean_irr2 << endl;
}

void Octree::buildAggregationRecursive()
{
    // Leaf nodes do nothing
    if (isLeafNode()) {
        return;
    }

    // Build children before computing things based on their results
    for (int i = 0; i < 8; ++i) {
        children[i]->buildAggregationRecursive();
    }

    // Compute the weighted center of mass and total power
    vec3 CoM(0, 0, 0), meanIrradiance(0, 0, 0);
    vec3 meanNormal(0, 0, 0), meanDirection(0, 0, 0);
    float weightSum = 0.0f;
    float totalArea = 0;
    for (int i = 0; i < 8; ++i) {
        if (children[i]->payload == nullptr) continue;
        float weight_i = luminance(children[i]->payload->meanIrradiance) * children[i]->payload->sumArea;
        weightSum += weight_i;
        CoM += children[i]->payload->centroid * weight_i; // power weighted mean position -> center of mass
        meanIrradiance += children[i]->payload->meanIrradiance * children[i]->payload->sumArea;
        totalArea += children[i]->payload->sumArea;
        meanNormal += children[i]->payload->weightedNormal * weight_i;
        meanDirection += children[i]->payload->weightedLightDir * weight_i;
    }

    if (weightSum > 0) {
        CoM = CoM * (1.0f / weightSum);
        meanNormal = meanNormal.normalize();         // ??
        meanDirection = meanDirection.normalize();   // ??
    }

    if (totalArea > 0) {
        meanIrradiance /= totalArea;
    }


    // why zero normal and direction ??
    if (0 == meanNormal.length()) {
        //             printf("zero normal\n");
        //             printf("weight sum: %f\n", weightSum);
        meanNormal.x = 1;
    }
    if (0 == meanDirection.length()) {
        //             printf("zero direction\n");
        //             printf("weight sum: %f\n", weightSum);
        meanDirection.x = 1;
    }


    // Create aggregate
    if (nullptr != payload) delete payload; //PREVENTS MEMORY LEAKAGE
    payload = new OctreeSample(CoM, meanNormal, meanDirection, meanIrradiance, totalArea);
}

float Octree::luminance(const vec3& c) const
{
    return ::luminance(c);
}

float Octree::weightingFucntion(const vec3& position, const vec3& normal, const vec3& direction, const OctreeSample* sample) const
{
    vec3 d = sample->centroid - position;
//     return 1.0 / d.lengthSquared();

    //////////////////////////////////////////////////////////////////////////
    vec3 z = vec3(0, 1, 0);
    float ret = /*luminance(sample->meanIrradiance) * */sample->sumArea / d.lengthSquared();
//     ret *= luminance(sample->meanIrradiance);
//     ret *= luminance(bssrdf(sample->centroid, sample->weightedNormal, sample->weightedLightDir, position, normal, direction));

//     cout << ret << endl;
    return ret;
}

vec3 Octree::fastsumRadiance(const vec3& p, const vec3& nl, const vec3& dir, float traversal_threshold, std::function<vec3(const vec3&, const vec3&, const vec3&, const vec3&, const vec3&, const vec3&)> kernel) const
{
    int debug_count = 0;
    static int debug_max_count = 0;
    float inv_overall_mean_irradiance = luminance(this->payload->meanIrradiance);

    vec3 f(0, 0, 0);

    const Octree *todo[64];
    int stack = 0; // points to the topmost valid element

    // push root
    todo[stack] = (this);

    while (stack >= 0)
    {
        debug_count++;

        const Octree *item = todo[stack--]; // pop item

        // Leaf?
        if (item->isLeafNode())
        {
            if (item->payload != nullptr)
            {
                f += item->payload->getRadiantExitance(p, nl, dir, kernel);
            }
        }
        else
        {
            for (int n = 0; n < 8; ++n)
            {
                Octree *Child = item->children[n]; //C cannot possibly be NULL, because item is interior

                vec3 &cmin(Child->bmin), &cmax(Child->bmax);
                bool inside = p.x < cmax.x && p.y < cmax.y && p.z < cmax.z && p.x >= cmin.x && p.y >= cmin.y && p.z >= cmin.z;

                if (inside)
                {
                    todo[++stack] = Child; // note that the for-loop for the children may continue
                }
                else if (Child->payload != nullptr)
                {
                    // We're not in this quadrant. See if we should recurse anyway, or just
                    // take the aggregate value. if distance_to_CoM * theta > L then use aggregate!
                    //                         const float L_square = (sq(Child->halfDimension.x) +
                    //                                                 sq(Child->halfDimension.y) +
                    //                                                 sq(Child->halfDimension.z)) * 4.0f;

                    // note that this must be dimensionless
                    float weight = luminance(Child->payload->meanIrradiance) * inv_overall_mean_irradiance * weightingFucntion(p, nl, dir, Child->payload);
                    if (weight < traversal_threshold)
                    {
                        f += Child->payload->getRadiantExitance(p, nl, dir, kernel);
                    }
                    else
                    {
                        todo[++stack] = Child;
                    }
                }
            }
        }
    }

    //     if (debug_count > debug_max_count)
    //     {
    //         debug_max_count = debug_count;
    //         cout << "[ " << debug_max_count << " ]" << endl;
    //     }

    return f;
}

void Octree::insert(OctreeSample* point)
{
    // If this node doesn't have a point yet assigned 
    // and it is a leaf, then we're done!
    if (isLeafNode()) {
        if (payload == nullptr) {
            payload = point;
            return;

        }
        else {
            // We're at a leaf, but there's already something here
            // We will split this node and then insert the old data and this new node

            // Save this data point that was here for a later re-insert
            OctreeSample *oldPoint = payload;

            // Split the node
            for (int i = 0; i < 8; ++i) {
                // Compute new bounding box
                vec3 newOrigin = origin;
                newOrigin.x += halfDimension.x * (i & 4 ? .5f : -.5f);
                newOrigin.y += halfDimension.y * (i & 2 ? .5f : -.5f);
                newOrigin.z += halfDimension.z * (i & 1 ? .5f : -.5f);
                children[i] = new Octree(newOrigin, halfDimension*.5f);
            }
            payload = nullptr;

            // re-insert the old point, and insert this new point
            insert(oldPoint);
            insert(point);
        }

    }
    else {
        // We are at an interior node. Insert recursively into a child
        int octant = getOctantContainingPoint(point->centroid);
        children[octant]->insert(point);
    }
}

int Octree::getOctantContainingPoint(const vec3& point) const
{
    int oct = 0;
    if (point.x >= origin.x) oct |= 4;
    if (point.y >= origin.y) oct |= 2;
    if (point.z >= origin.z) oct |= 1;
    return oct;
}

bool Octree::isLeafNode() const
{
    return children[0] == nullptr;
}

Octree::~Octree()
{
    // SOLVES MEMORY LEAKAGE
    if (nullptr != payload) {
        delete payload;
    }
    for (int i = 0; i < 8; ++i) {
        if (children[i] != nullptr)
            delete children[i];
    }
}

Octree::Octree(const vec3& origin, const vec3& halfDimension) : origin(origin), halfDimension(halfDimension)
{
    for (int i = 0; i < 8; ++i)
        children[i] = nullptr;
    payload = nullptr;

    const vec3 margin = vec3(1, 1, 1) * .0f;
    bmin = origin - halfDimension - margin;
    bmax = origin + halfDimension + margin;
}
