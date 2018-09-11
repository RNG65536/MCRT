#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "camera.h"
#include "geometry.h"
#include "lightpath.h"
#include "numeric.h"

Ray CCamera::makeRay(float mx, float my) const
{
    mx = mx * 2.0f - 1.0f;
    my = my * 2.0f - 1.0f;
    glm::vec4 a = glm::vec4(mx, my, -1.0f, 1.0f); // near
    glm::vec4 b = glm::vec4(mx, my, 1.0f, 1.0f); // far
    a = m_invProjectionMatrix * a;
    b = m_invProjectionMatrix * b;
    a *= 1.0f / a.w;
    b *= 1.0f / b.w;
    a = m_invViewMatrix * a;
    b = m_invViewMatrix * b;
    return Ray(vec3(a.x, a.y, a.z), normalize(vec3(b.x - a.x, b.y - a.y, b.z - a.z)));
}

CCamera::CCamera(const vec3& origin, const vec3& lookat, const vec3& up, int film_width, int film_height, float fov_y)
{
    glm::mat4 view = glm::lookAt(
        glm::vec3(origin.x, origin.y, origin.z),
        glm::vec3(lookat.x, lookat.y, lookat.z),
        glm::vec3(up.x, up.y, up.z));
    glm::mat4 projection = glm::perspective(glm::radians(fov_y), float(film_width) / float(film_height), 0.01f, 1000.0f);

    m_viewMatrix = view;
    m_projectionMatrix = projection;
    m_invViewMatrix = glm::inverse(m_viewMatrix);
    m_invProjectionMatrix = glm::inverse(m_projectionMatrix);

    //////////////////////////////////////////////////////////////////////////

    m_film_width = film_width;
    m_film_height = film_height;
    m_film_world_width = float(film_width) / float(film_height);
    m_film_world_height = 1.0f;

    m_origin = origin;
    m_dist = m_film_world_height / (2.0 * tan(toRadians(fov_y * 0.5)));
    m_forward = normalize(lookat - origin);
    m_frame_x = normalize(cross(m_forward, up));
    m_frame_y = cross(m_frame_x, m_forward);
}

bool CCamera::rasterizePrimaryRay(const vec3& direction, float& px, float& py) const
{
    vec3 screen_center = m_origin + m_forward * m_dist;
    vec3 screen_position = m_origin + (direction *
        (m_dist / dot(direction, m_forward))) - screen_center;

    // in pixel space
    px = dot(m_frame_x, screen_position) * (m_film_width / m_film_world_width) + (m_film_width * 0.5);
    py = dot(m_frame_y, screen_position) * (m_film_height / m_film_world_height) + (m_film_height * 0.5);

    return px >= 0 && px < m_film_width && py >= 0 && py < m_film_height;
}

bool CCamera::rasterizePoint(const vec3& position, float& px, float& py) const
{
    vec3 direction = normalize(position - m_origin);

    return rasterizePrimaryRay(direction, px, py);
}

float CCamera::measurementFunction(const Vert& current, const Vert& next) const
{
    float w = 1.0 / (float(m_film_world_width) * float(m_film_world_height));

    vec3 ray_dir = normalize(next.p - current.p);
    float cos_theta = dot(ray_dir, m_forward);
    float distance_to_screen_2 = m_dist / cos_theta;
    distance_to_screen_2 = distance_to_screen_2 * distance_to_screen_2;
    w *= 1.0 / (cos_theta / distance_to_screen_2);

    w *= fabs(directionToArea(current, next));

    return w;
}

float CCamera::measurementFunction(const Vert& v) const
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
    vec3 ray_dir = normalize(v.p - m_origin);
    float cos_theta = dot(ray_dir, m_forward);
    float dist2 = m_dist / cos_theta; // distance from lens to pixel
    dist2 *= dist2;

    float film_pdfA = 1.0 / (float(m_film_world_width) * float(m_film_world_height));
    float w = __pdfAtoW(film_pdfA, dist2, cos_theta);

    w = __pdfWtoA(w, dot(m_origin - v.p, m_origin - v.p),
        f_max(dot(v.n, normalize(m_origin - v.p)), 0));

    return w;
#endif
}

float CCamera::measurementFunction(const vec3& position, const vec3& normal) const
{
    vec3 ray_dir = normalize(position - m_origin);
    float cos_theta = dot(ray_dir, m_forward);
    float dist2 = m_dist / cos_theta; // distance from lens to pixel
    dist2 *= dist2;

    float film_pdfA = 1.0 / (float(m_film_world_width) * float(m_film_world_height));
    float w = __pdfAtoW(film_pdfA, dist2, cos_theta);

    w *= pdfWtoA(m_origin, position, normal);

    return w;
}

Ray CCamera::sample(float rnd1, float rnd2) const
{
    const vec3d su = vec3d(m_frame_x * -(0.5 - rnd1) * m_film_world_width);
    const vec3d sv = vec3d(m_frame_y * (0.5 - rnd2) * m_film_world_height);
    const vec3d sw = vec3d(m_forward * m_dist);
    return Ray(m_origin, normalize(su + sv + sw).toFloat());
}

vec3 CCamera::getForwardDirection() const
{
    return m_forward;
}

int CCamera::getFilmHeight() const
{
    return m_film_height;
}

int CCamera::getFilmWidth() const
{
    return m_film_width;
}

float CCamera::getFilmDistance() const
{
    return m_dist;
}

vec3 CCamera::getPosition() const
{
    return m_origin;
}

float CCamera::pixelPdfA() const
{
    return getFilmWidth() * getFilmHeight() * filmPdfA();
}

float CCamera::filmPdfA() const
{
    return 1.0 / (float(m_film_world_width) * float(m_film_world_height));
}
