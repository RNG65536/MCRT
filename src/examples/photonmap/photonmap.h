#ifndef photonmap_h__
#define photonmap_h__

//------------------------------------------------------------------
// photonmap.cpp
// An example implementation of the photon map data structure
//
// Henrik Wann Jensen - February 2001
//------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

//
// This is the photon
// The power is not compressed so the
// size is 28 bytes
//
typedef struct Photon
{
    float pos[3];               // photon position
    short plane;                // splitting palne for kd-tree
    unsigned char theta, phi;   // incoming direction
    float dir[3];
    float power[3];             // photon power (uncompressed)
} Photon;

//
// This structure is used only to locate the
// nearest photons
//
typedef struct NearestPhotons
{
    int max;
    int found;
    int got_heap;
    float pos[3];
    float *dist2;
    const Photon **index;
} NearestPhotons;

//
// This is the Photon_map class
//
class PhotonMap
{
public:

    // This is the constructor for the photon map.
    // To create the photon map it is necessary to specify the
    // maximum number of photons that will be stored
    //
    PhotonMap(const int max_phot);
    ~PhotonMap();

    // store puts a photon into the flat array that will form
    // the final kd-tree
    // Call this function to store a photon
    //
    void store(
        const float power[3],       // photon power
        const float pos[3],         // photon position
        const float dir[3]);        // photon direction

    // scale_photon_power is used to scale the power of all
    // photons once they have been emitted from the light
    // source. scale = 1 / (#emitted photons)
    // Call this function after each light source is processed
    //
    void scale_photon_power(
        const float scale);         // 1/(number of emitted photons)

    // balance creates a left-balanced kd-tree from the flat photon array
    // This function should be called before the photon map
    // is used for rendering
    //
    void balance(void);             // balance the kd-tree (before use!)

    // irradiance_estimate computes an irradiance estimate
    // at a given surface position
    //
    void irradiance_estimate(
        float irrad[3],             // returned irradiance
        const float pos[3],         // surface position
        const float normal[3],      // surface normal at pos
        const float max_dist,       // max distance to look for photons
        const int nphotons) const;  // number of photons to use

    // locate_photons finds the nearest photons in the
    // photon map given the parameters in np
    //
    void locate_photons(
        NearestPhotons *const np,   // np is used to locate the photons
        const int index) const;     // the photons

    // photon_dir returns the direction of a photon
    //
    void photon_dir(
        float *dir,                 // direction of photon (returned)
        const Photon *p) const;     // the photon

private:

    // See "Realistic Image Synthesis using Photon Mapping" Chapter 6
    // for an explanation of this function
    //
    void balance_segment(
        Photon **pbal,
        Photon **porg,
        const int index,
        const int start,
        const int end);

    // median_split spits the photon array into two separate
    // pieces around the median, with all photons below
    // the median in the lower half and all photons above
    // the median in the upper half. The comparison
    // criteria is the axis (indicated by the axis parameter)
    // (inspired by routine in "Algorithm in C++" by Sedgewick)
    //
    void median_split(
        Photon **p,
        const int start,
        const int end,
        const int median,
        const int axis);

    Photon *photons;

    int stored_photons;
    int half_stored_photons;
    int max_photons;
    int prev_scale;

    float costheta[256];
    float sintheta[256];
    float cosphi[256];
    float sinphi[256];

    float bbox_min[3];      // use bbox_min
    float bbox_max[3];      // use bbox_max
};

#endif // photonmap_h__
