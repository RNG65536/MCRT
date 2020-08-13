#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include "film.h"
#include "numeric.h"

void FrameBuffer::dumpPPM(const char* filename)
{
    // Save result to a PPM image (keep these flags if you compile under
    // Windows) NOTE::especially std:ios::binary which is equivalent to "wb" in
    // fprintf()
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    if (!ofs.is_open())
    {
        std::cout << "cannot write to " << filename << std::endl;
        return;
    }

    ofs << "P6\n" << m_width << " " << m_height << "\n255\n";
    for (int j = m_height - 1; j >= 0; --j)
    {
        for (int i = 0; i < m_width; ++i)
        {
            int offset = (i + j * m_width) * 3;
            ofs << (unsigned char)(f_min(1.0f, m_buffer[offset + 0]) * 255)
                << (unsigned char)(f_min(1.0f, m_buffer[offset + 1]) * 255)
                << (unsigned char)(f_min(1.0f, m_buffer[offset + 2]) * 255);
        }
    }
    ofs.close();
}

struct HDRPixel
{
    uint8_t r, g, b, e;
    HDRPixel() = default;
    HDRPixel(uint8_t r, uint8_t g, uint8_t b, uint8_t e)
        : r(r), g(g), b(b), e(e)
    {
    }
    uint8_t operator[](int i) const
    {
        return (&r)[i];
    }
};

HDRPixel toRGBE(const vec3& c)
{
    float d = f_max(f_max(c.x, c.y), c.z);
    if (d <= 1e-32)
    {
        return HDRPixel(0, 0, 0, 0);
    }
    int   e;
    float m = frexp(d, &e);  // d = m * 2^e
    d = m * 256.0f / d;
    return HDRPixel(static_cast<uint8_t>(c.x * d),
                    static_cast<uint8_t>(c.y * d),
                    static_cast<uint8_t>(c.z * d),
                    static_cast<uint8_t>(e + 128));
}

void FrameBuffer::dumpHDR(const char* filename)
{
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    if (!ofs.is_open())
    {
        std::cout << "cannot write to " << filename << std::endl;
    }

    ofs << "#?RADIANCE" << std::endl;
    ofs << "# Made with custom writer" << std::endl;
    ofs << "FORMAT=32-bit_rle_rgbe" << std::endl;
    ofs << "EXPOSURE=1.0" << std::endl;
    ofs << std::endl;

    ofs << "-Y " << m_height << " +X " << m_width << std::endl;

    for (int j = m_height - 1; j >= 0; --j)
    {
        std::vector<HDRPixel> line(m_width);
        for (int i = 0; i < m_width; i++)
        {
            int      offset = (i + j * m_width) * 3;
            HDRPixel p = toRGBE(vec3(
                m_buffer[offset], m_buffer[offset + 1], m_buffer[offset + 2]));
            line[i] = p;
        }
        ofs << uint8_t(2) << uint8_t(2);
        ofs << uint8_t((m_width >> 8) & 0xFF) << uint8_t(m_width & 0xFF);
        for (int k = 0; k < 4; k++)
        {
            for (int cursor = 0; cursor < m_width;)
            {
                const int cursor_move = std::min(127, m_width - cursor);
                ofs << uint8_t(cursor_move);
                for (int i = cursor; i < cursor + cursor_move; i++)
                {
                    ofs << uint8_t(line[i][k]);
                }
                cursor += cursor_move;
            }
        }
    }
    ofs.close();
}

void FrameBuffer::tonemapReinhard()
{
    int    pixel_count = m_total / 3;
    float* lum = new float[pixel_count];
    float* fb = &m_buffer[0];

    float lum_eps = 1e-7f;

    for (int n = 0; n < pixel_count; n++)
    {
        lum[n] = ::luminance(vec3(fb[n * 3], fb[n * 3 + 1], fb[n * 3 + 2]));
        if (lum[n] < lum_eps) lum[n] = lum_eps;
    }

    float lum_min = std::numeric_limits<float>::max();
    float lum_max = -std::numeric_limits<float>::max();
    for (int n = 0; n < pixel_count; n++)
    {
        lum_min = lum_min < lum[n] ? lum_min : lum[n];
        lum_max = lum_max > lum[n] ? lum_max : lum[n];
    }

    float l_logmean = 0;
    float l_mean = 0;
    float r_mean = 0;
    float g_mean = 0;
    float b_mean = 0;
    for (int n = 0; n < pixel_count; n++)
    {
        l_logmean += logf(lum[n]);
        l_mean += lum[n];
        r_mean += fb[n * 3];
        g_mean += fb[n * 3 + 1];
        b_mean += fb[n * 3 + 2];
    }
    l_logmean /= pixel_count;
    l_mean /= pixel_count;
    r_mean /= pixel_count;
    g_mean /= pixel_count;
    b_mean /= pixel_count;

    float lmin = logf(lum_min);
    float lmax = logf(lum_max);
    float k = (lmax - l_logmean) / (lmax - lmin);
    float m0 = 0.3f + 0.7f * powf(k, 1.4f);  // contrast
    m0 = 0.77f;                              // hdrsee default

    float m = m0;  // Contrast [0.3f, 1.0f]
    //     printf("contrast: %f\n", m);

    float c = 0.5;  // Chromatic Adaptation  [0.0f, 1.0f]
    float a = 0;    // Light Adaptation  [0.0f, 1.0f]
    float f = 0;  // Intensity [-35.0f, 10.0f] (void*)func = intuitiveintensity
                  // specify by log scale

    f = expf(-f);

    for (int n = 0; n < pixel_count; n++)
    {
        float r(fb[n * 3]), g(fb[n * 3 + 1]), b(fb[n * 3 + 2]);

        float r_lc = c * r + (1.0f - c) * lum[n];       // local adaptation
        float r_gc = c * r_mean + (1.0f - c) * l_mean;  // global adaptation
        float r_ca = a * r_lc + (1.0f - a) * r_gc;      // pixel adaptation

        float g_lc = c * g + (1.0f - c) * lum[n];       // local adaptation
        float g_gc = c * g_mean + (1.0f - c) * l_mean;  // global adaptation
        float g_ca = a * g_lc + (1.0f - a) * g_gc;      // pixel adaptation

        float b_lc = c * b + (1.0f - c) * lum[n];       // local adaptation
        float b_gc = c * b_mean + (1.0f - c) * l_mean;  // global adaptation
        float b_ca = a * b_lc + (1.0f - a) * b_gc;      // pixel adaptation

        r = r / (r + powf(f * r_ca, m));
        g = g / (g + powf(f * g_ca, m));
        b = b / (b + powf(f * b_ca, m));

        fb[n * 3] = r;
        fb[n * 3 + 1] = g;
        fb[n * 3 + 2] = b;
    }

    delete[] lum;
}

void FrameBuffer::tonemapGamma(float gamma)
{
    float inv_gamma = 1.0f / gamma;
    for (int n = 0; n < m_total; n++)
    {
        m_buffer[n] = powf(clampf(m_buffer[n], 0.0f, 1.0f), inv_gamma);
    }
}

void FrameBuffer::accumulatePixel(int i, int j, const vec3& c)
{
    if (i < 0 || i >= m_width || j < 0 || j >= m_height)
    {
        return;
    }

    int offset = (i + j * m_width) * 3;
    m_buffer[offset + 0] += c.x;
    m_buffer[offset + 1] += c.y;
    m_buffer[offset + 2] += c.z;
}

void FrameBuffer::accumulateBuffer(const FrameBuffer& f)
{
    std::transform(m_buffer.begin(),
                   m_buffer.end(),
                   f.m_buffer.begin(),
                   m_buffer.begin(),
                   std::plus<float>());
}

void FrameBuffer::resize(int w, int h)
{
    m_width = w;
    m_height = h;
    m_total = m_width * m_height * 3;
    m_buffer.resize(m_total, 0);
}

FrameBuffer::FrameBuffer(int w, int h)
{
    this->resize(w, h);
}

FrameBuffer::FrameBuffer() : FrameBuffer(0, 0)
{
}

FrameBuffer::~FrameBuffer()
{
}

vec3 FrameBuffer::pixel(int i, int j) const
{
    return vec3(m_buffer[(i + j * m_width) * 3 + 0],
                m_buffer[(i + j * m_width) * 3 + 1],
                m_buffer[(i + j * m_width) * 3 + 2]);
}

void FrameBuffer::scale(float s)
{
    for (auto& v : m_buffer)
    {
        v *= s;
    }
}

int FrameBuffer::height() const
{
    return m_height;
}

int FrameBuffer::width() const
{
    return m_width;
}
