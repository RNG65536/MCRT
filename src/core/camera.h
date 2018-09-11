#ifndef camera_h__
#define camera_h__

#include <glm/glm.hpp>
#include "vectors.h"
#include "ray.h"

// camera for emitting primary rays

class Vert;

class CCamera
{
private:
    glm::mat4 m_viewMatrix;
    glm::mat4 m_projectionMatrix;
    glm::mat4 m_invViewMatrix;
    glm::mat4 m_invProjectionMatrix;

    int m_film_width, m_film_height;
    float m_film_world_width;
    float m_film_world_height;
    float m_dist;
    vec3 m_origin, m_frame_x, m_frame_y, m_forward;

public:
    CCamera(const vec3& origin, const vec3& lookat, const vec3& up, int film_width, int film_height, float fov_y);
    Ray makeRay(float mx, float my) const;
    Ray sample(float rnd1, float rnd2) const;

    // with perspective correction ??
    float measurementFunction(const Vert& current, const Vert& next) const;
    float measurementFunction(const Vert& v) const;
    float measurementFunction(const vec3& position, const vec3& normal) const;
    bool rasterizePrimaryRay(const vec3& direction, float& px, float& py) const;
    bool rasterizePoint(const vec3& position, float& px, float& py) const;
    int getFilmWidth() const;
    int getFilmHeight() const;
    float getFilmDistance() const;
    vec3 getForwardDirection() const;
    vec3 getPosition() const;
    float filmPdfA() const;
    float pixelPdfA() const;
};

#endif // camera_h__
