#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "camera.h"
#include "geometry.h"
#include "lightpath.h"
#include "numeric.h"

Camera::Camera(const vec3& origin,
               const vec3& lookat,
               const vec3& up,
               int         film_width,
               int         film_height,
               float       fov_y)
{
    glm::mat4 view = glm::lookAt(glm::vec3(origin.x, origin.y, origin.z),
                                 glm::vec3(lookat.x, lookat.y, lookat.z),
                                 glm::vec3(up.x, up.y, up.z));
    glm::mat4 projection =
        glm::perspective(glm::radians(fov_y),
                         float(film_width) / float(film_height),
                         0.01f,
                         1000.0f);

    m_film_width = film_width;
    m_film_height = film_height;
    m_film_world_width = float(film_width) / float(film_height);
    m_film_world_height = 1.0f;  // this is arbitrary

    m_origin = origin;
    m_dist = m_film_world_height / (2.0f * tan(toRadians(fov_y * 0.5f)));
    m_forward = normalize(lookat - origin);
    m_frame_x = normalize(cross(m_forward, up));
    m_frame_y = cross(m_frame_x, m_forward);
}

Ray Camera::makeRay(float mx, float my) const
{
    const vec3 su = m_frame_x * (mx - 0.5f) * m_film_world_width;
    const vec3 sv = m_frame_y * (my - 0.5f) * m_film_world_height;
    const vec3 sw = m_forward * m_dist;
    return Ray(m_origin, normalize(su + sv + sw));
}

bool Camera::rasterizePrimaryRay(const vec3& direction,
                                 float&      px,
                                 float&      py) const
{
    vec3 screen_center = m_origin + m_forward * m_dist;
    vec3 screen_position = m_origin +
                           (direction * (m_dist / dot(direction, m_forward))) -
                           screen_center;

    // in pixel space
    px = dot(m_frame_x, screen_position) * (m_film_width / m_film_world_width) +
         (m_film_width * 0.5f);
    py = dot(m_frame_y, screen_position) *
             (m_film_height / m_film_world_height) +
         (m_film_height * 0.5f);

    return px >= 0 && px < m_film_width && py >= 0 && py < m_film_height;
}

bool Camera::rasterizePoint(const vec3& position, float& px, float& py) const
{
    vec3 direction = normalize(position - m_origin);

    return rasterizePrimaryRay(direction, px, py);
}

float Camera::measurementFunction(const Vert& current, const Vert& next) const
{
    float w = 1.0f / (float(m_film_world_width) * float(m_film_world_height));

    vec3  ray_dir = normalize(next.p - current.p);
    float cos_theta = dot(ray_dir, m_forward);
    float distance_to_screen_2 = m_dist / cos_theta;
    distance_to_screen_2 = distance_to_screen_2 * distance_to_screen_2;
    w *= 1.0f / (cos_theta / distance_to_screen_2);

    w *= fabs(directionToArea(current, next));

    return w;
}

float Camera::measurementFunction(const Vert& v) const
{
#if 0
    vec3 ray_dir = normalize(v.p - m_origin);
    float cos_theta = dot(ray_dir, m_forward);
    float dist2 = m_dist / cos_theta; // distance from lens to pixel
    dist2 *= dist2;

    float film_pdfA = 1.0 / (float(m_film_world_width) * float(m_film_world_height));
    float w = __pdfAtoW(film_pdfA, dist2, cos_theta);

    w = __pdfWtoA(w, dot(m_origin - v.p, m_origin - v.p),
        dot(v.n, normalize(m_origin - v.p)));

    return w;
#else
    vec3  ray_dir = normalize(v.p - m_origin);
    float cos_theta = dot(ray_dir, m_forward);
    float dist2 = m_dist / cos_theta;  // distance from lens to pixel
    dist2 *= dist2;

    float film_pdfA =
        1.0f / (float(m_film_world_width) * float(m_film_world_height));
    float w = __pdfAtoW(film_pdfA, dist2, cos_theta);

    w = __pdfWtoA(w,
                  dot(m_origin - v.p, m_origin - v.p),
                  f_max(dot(v.n, normalize(m_origin - v.p)), 0));

    return w;
#endif
}

float Camera::measurementFunction(const vec3& position,
                                  const vec3& normal) const
{
    vec3  ray_dir = normalize(position - m_origin);
    float cos_theta = dot(ray_dir, m_forward);
    float dist2 = m_dist / cos_theta;  // distance from lens to pixel
    dist2 *= dist2;

    float film_pdfA =
        1.0f / (float(m_film_world_width) * float(m_film_world_height));
    float w = __pdfAtoW(film_pdfA, dist2, cos_theta);

    w *= pdfWtoA(m_origin, position, normal);

    return w;
}

vec3 Camera::forwardDirection() const
{
    return m_forward;
}

int Camera::filmHeight() const
{
    return m_film_height;
}

int Camera::filmWidth() const
{
    return m_film_width;
}

float Camera::filmDistance() const
{
    return m_dist;
}

vec3 Camera::position() const
{
    return m_origin;
}

float Camera::pixelPdfA() const
{
    return m_film_width * m_film_height * filmPdfA();
}

float Camera::filmPdfA() const
{
    return 1.0f / (float(m_film_world_width) * float(m_film_world_height));
}
