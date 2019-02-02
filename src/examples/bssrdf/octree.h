#ifndef octree_h__
#define octree_h__

#include <functional>
#include "vectors.h"

// TODO: to actually store area for all levels of irradiance samples,
// and use sumArea / L^2 to approximate the solid angle subtended by the levels

// #define DIVIDE_BY_FOUR_PI /(4.0*M_PI)
// #define MULTIPLY_BY_FOUR_PI *4*M_PI // * 4 //* (4 * PI)
// #define MULTIPLY_BY_FOUR_PI //*4*M_PI // * 4 //* (4 * PI)

//////////////////////////////////////////////////////////////////////////
// assuming uniform sampling, so area weighting need not be considered
struct OctreeSample {
    // to evaluate dirpole bssrdf, need { position, normal, direction }

    vec3 centroid;          // position
    vec3 meanIrradiance;  // average irradiance
    float sumArea;

    vec3 weightedNormal;    // normal
    vec3 weightedLightDir;  // direction

    OctreeSample(const vec3& position, const vec3& normal, const vec3& direction, const vec3& rad, float area);

    vec3 getRadiantExitance(vec3 position, vec3 normal, vec3 direction /* view ray */, std::function<vec3(const vec3&, const vec3&, const vec3&, const vec3&, const vec3&, const vec3&)> kernel);
};

class Octree {
    vec3 origin, halfDimension;
    vec3 bmin, bmax;
    Octree *children[8]; // Children
    OctreeSample* payload; // Possible use : 1) instances stored at a leaf 2) aggregation at interior nodes

    /*
    child:	0 1 2 3 4 5 6 7
    x:      - - - - + + + +
    y:      - - + + - - + +
    z:      - + - + - + - +
    */

    Octree(const Octree& tree);
    int getOctantContainingPoint(const vec3& point) const;
    bool isLeafNode() const;
    float luminance(const vec3& c) const;
    float weightingFucntion(const vec3& position, const vec3& normal, const vec3& direction, const OctreeSample* sample) const;

    void buildAggregationRecursive();
    float meanIrradianceRecursive() const;

public:
    Octree(const vec3& origin, const vec3& halfDimension);
    ~Octree();
    void insert(OctreeSample* point);
    void finalize();

    //Barnes-Hut algorithm
    vec3 fastsumRadiance(const vec3& p, const vec3& nl, const vec3& dir, float traversal_threshold, std::function<vec3(const vec3&, const vec3&, const vec3&, const vec3&, const vec3&, const vec3&)> kernel) const;
};

#endif // octree_h__
