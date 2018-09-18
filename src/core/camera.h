#pragma once

#include "ray.h"
#include "vectors.h"

// pinhole camera model

class Vert;

class Camera
{
public:
    Camera(const vec3& origin,
           const vec3& lookat,
           const vec3& up,
           int         film_width,
           int         film_height,
           float       fov_y);

    // map from unit [0, 1]^2 to camera rays
    Ray makeRay(float mx, float my) const;

    // with perspective correction ??
    float measurementFunction(const Vert& current, const Vert& next) const;
    float measurementFunction(const Vert& v) const;
    float measurementFunction(const vec3& position, const vec3& normal) const;
    bool rasterizePrimaryRay(const vec3& direction, float& px, float& py) const;
    bool rasterizePoint(const vec3& position, float& px, float& py) const;

    int   filmWidth() const;
    int   filmHeight() const;
    float filmDistance() const;
    vec3  forwardDirection() const;
    vec3  position() const;
    float filmPdfA() const;
    float pixelPdfA() const;

private:
    int   m_film_width;
    int   m_film_height;
    float m_film_world_width;
    float m_film_world_height;
    float m_dist;
    vec3  m_origin, m_frame_x, m_frame_y, m_forward;
};
