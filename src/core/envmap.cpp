#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include "envmap.h"
#include "numeric.h"
#include "sample.h"

#define NUM_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923

// #define LIGHT_COUNT 100000
// #define LIGHT_COUNT 10000
// #define LIGHT_COUNT 1000 //++++++++++++++++++ best practice for full summation
// #define LIGHT_COUNT 100
#define  LIGHT_COUNT 0

static const bool is_pre_warped = true;

static inline float f_max(float a, float b) { return (((a) > (b)) ? (a) : (b)); }
static inline float f_min(float a, float b) { return (((a) < (b)) ? (a) : (b)); }

// static float randf() { return rand() / (RAND_MAX + 1.0f); }
// static float randf(){
//     return (float(rand()) + float(rand()) * float(RAND_MAX)) / (float(RAND_MAX) * float(RAND_MAX));
// }
static int clampi(int x, int a, int b) { return x < a ? a : x > b ? b : x; }
float EnvmapLoader::luminance(float r, float g, float b) {
    return 0.299 * r + 0.587 * g + 0.114 * b;
}

// DirLight EnvmapLoader::getLightRandom() {
//     return m_lights[clampi(randf() * m_lights.size(), 0, m_lights.size() - 1)];
// }

// DirLight EnvmapLoader::getLight(int n) {
//     return m_lights[n];
// }

static float dir_to_theta(float dir[3]) {
    //     float theta = atan(dir[0] / dir[2]) + M_PI_2;
    //     if ( dir[2] < 0 ) theta += M_PI;
    float theta = atan(dir[2] / dir[0]) + M_PI_2;
    if (dir[0] < 0) theta += NUM_PI;
    //     float theta = atan(dir[2] / dir[0]) + M_PI_2;
    //     if ( dir[0] > 0 ) theta += M_PI;
    return theta;
}

void EnvmapLoader::analyze()
{
    // visual debugging
    if (0)
    {
        // clear with black
        m_buffer.assign(m_header.width * m_header.height * 3, 0);

        for (int n = 0; n < 100000; n++) {
            float u, v;
            sampleByPdf(randf(), randf(), u, v);

            int u_i = u * m_header.width;
            int v_i = v * m_header.height;

            u_i = clampi(u_i, 0, m_header.width - 1);
            v_i = clampi(v_i, 0, m_header.height - 1);

            int offset = (u_i + v_i * m_header.width) * 3;
            m_buffer[offset + 0] += 0.1;
            m_buffer[offset + 1] += 0.1;
            m_buffer[offset + 2] += 0.1;
        }

        toPPM("__analyze.ppm");
    }

    // test dir_to_theta()
    //     for (int n = 0; n < 1000; n++) {
    // 
    //         float phi = M_PI_2;
    //         float theta = randf() * 2 * M_PI;
    // 
    //         float dir[3];
    //         dir[0] = sin(phi) * sin(theta);
    //         dir[1] = cos(phi);
    //         dir[2] = sin(phi) * -cos(theta);
    // 
    //         printf("theta: %f -> %f\n", theta, dir_to_theta(dir));
    //     }
}

EnvmapLoader::EnvmapLoader(const std::string filepath)
{
    FILE *fp = fopen(filepath.c_str(), "rb");
    if (!fp)
    {
        printf("Cannot load envmap %s\n", filepath.c_str());
        return;
    }
    printf("Loaded envmap %s\n", filepath.c_str());

    fread(&m_header, 1, sizeof(EnvmapHeader), fp);
    int total = m_header.width * m_header.height * 3;
    m_buffer.resize(total);
    fread(&m_buffer[0], sizeof(float), total, fp);

    fclose(fp);


    buildCdf();
//     generate_light_samples();

    analyze();
}

EnvmapLoader::EnvmapLoader(const std::string filepath, int width, int height)
{
    FILE *fp = fopen(filepath.c_str(), "rb");
    if (!fp)
    {
        printf("Cannot load envmap %s\n", filepath.c_str());
        return;
    }
    printf("Loaded envmap %s\n", filepath.c_str());

    m_header.width = width;
    m_header.height = height;
    int total = m_header.width * m_header.height * 3;
    m_buffer.resize(total);
    fread(&m_buffer[0], sizeof(float), total, fp);

    fclose(fp);


    buildCdf();
//     generate_light_samples();

    analyze();
}

void EnvmapLoader::scale(float x)
{
    for (int n = 0; n < m_buffer.size(); n++) {
        m_buffer[n] *= x;
    }

    std::vector<float> lum(m_buffer.size() / 3);
    for (int i = 0; i < m_buffer.size(); i += 3)
    {
        lum[i / 3] = luminance(m_buffer[i + 0], m_buffer[i + 1], m_buffer[i + 2]);
    }
    std::cout << "envmap max = " << *std::max_element(lum.begin(), lum.end()) << std::endl;
}

void EnvmapLoader::toPPM(const char *filename)
{
    // Save result to a PPM image (keep these flags if you compile under Windows)    
    // NOTE::especially std:ios::binary which is equivalent to "wb" in fprintf()
    std::ofstream ofs(filename, std::ios::out | std::ios::binary);
    ofs << "P6\n" << m_header.width << " " << m_header.height << "\n255\n";
    int total = m_header.width * m_header.height * 3;
    for (int n = 0; n < total; n++) {
        ofs << (unsigned char)(f_min(1.0f, m_buffer[n]) * 255);
    }
    ofs.close();
}

// void EnvmapLoader::generate_light_samples()
// {
//     for (int n = 0; n < LIGHT_COUNT; n++) {
//         DirLight dl = sample();
//         m_lights.push_back(dl);
//     }
// 
//     struct compare_lights {
//         compare_lights() {}
//         bool operator() (DirLight& a, DirLight& b) {
//             return luminance(a.rad[0], a.rad[1], a.rad[2]) >
//                 luminance(b.rad[0], b.rad[1], b.rad[2]);
//         }
//     };
//     std::sort(m_lights.begin(), m_lights.end(), compare_lights());
// 
//     printf("#samples on envmap: %d\n", m_lights.size());
// }

struct multiply_by {
    float s;
    multiply_by(float _s) : s(_s) {}
    float operator() (float x) {
        return x * s;
    }
};

void EnvmapLoader::buildCdf()
// the whole algorithm is written compactly here to implement the corresponding matlab code test_marginal_pdf.m
{
    int M = m_header.height;
    int N = m_header.width;

    std::vector<float> lum(M * N);
    for (int n = 0; n < lum.size(); n++) {
        lum[n] = luminance(getR(n), getG(n), getB(n));
    }
    m_luminance = lum;

    // pre-warp with sin(phi) such that samples are driven away from the poles
    if (is_pre_warped)
    {
        for (int m = 0; m < M; m++) {
            float phi = NUM_PI * (m + 0.5f) / M;
            for (int n = 0; n < N; n++) {
                lum[n + m * N] *= sin(phi); // should be M_PI - phi, but does not matter due to symmetry
            }
        }
        //         m_luminance_warped = lum;
    }

    // P(uv) (not P(v|u)) is proportional to lum
    m_normalizationFactor = 1.0 / std::accumulate(lum.begin(), lum.end(), 0.0f); // MUST use 0.0 or 0.0f instead of 0 !!!!!!!!!!
    //////////////////////////////////////////////////////////////////////////

//     m_pdf_u.resize(N);
//     m_pdf_uv.resize(M * N);
//     m_cdf_u.resize(N);
//     m_cdf_uv.resize(M * N);
// 
//     // prepare pdf_u
//     for (int n = 0; n < N; n++) {
//         float sum = 0;
//         for (int m = 0; m < M; m++) {
//             sum += lum[n + m * N];
//         }
//         m_pdf_u[n] = sum;
//     }
// 
//     // prepare pdf_uv and normalized
//     for (int m = 0; m < M; m++) {
//         for (int n = 0; n < N; n++) {
//             m_pdf_uv[n + m * N] = lum[n + m * N] / m_pdf_u[n];
//         }
//     }
// 
//     // then normalize pdf_u
//     float m_pdf_u_norm = 1.0f / std::accumulate(m_pdf_u.begin(), m_pdf_u.end(), 0.0f);
//     std::transform(m_pdf_u.begin(), m_pdf_u.end(), m_pdf_u.begin(), multiply_by(m_pdf_u_norm));
// 
//     // prepare cdf_u
//     float acc = 0;
//     for (int n = 0; n < N; n++) {
//         acc += m_pdf_u[n];
//         m_cdf_u[n] = acc;
//     }
// 
//     // prepare cdf_uv
//     for (int n = 0; n < N; n++) {
//         float acc = 0;
//         for (int m = 0; m < M; m++) {
//             acc += m_pdf_uv[n + m * N];
//             m_cdf_uv[n + m * N] = acc;
//         }
//     }
// 
//     printf("norm factor: %f\n", m_normalizationFactor);
//     float sum = 0;
//     double sum2 = 0;
//     for (int n = 0; n < N; n++) {
//         for (int m = 0; m < M; m++) {
//             sum += m_normalizationFactor * m_luminance[n + m * N];
//             sum2 += m_normalizationFactor * m_luminance[n + m * N] * float(m_header.width) * float(m_header.height);
// 
//             //             printf("pdf1 = %f\n", lum[n + m * N] * m_normalizationFactor);
//             //             printf("pdf2 = %f\n", m_pdf_u[n] * m_pdf_uv[n + m * N]);
//         }
//     }
//     printf("sum = %f\n", sum);
//     //     printf("sum2 = %f\n", sum2);
    //////////////////////////////////////////////////////////////////////////

    // M -- height, N -- width
    std::vector<std::vector<float>> columns(N);
    for (int n = 0; n < N; n++)
    {
        columns[n].clear();
        columns[n].resize(M);
        for (int m = 0; m < M; m++)
        {
            columns[n][m] = lum[n + m * N];
        }
    }

    std::vector<float> column_sum(N);
    for (int n = 0; n < N; n++)
    {
        column_sum[n] = std::accumulate(columns[n].begin(), columns[n].end(), 0.0f);
    }
    m_sampleU = new DiscreteSampler(column_sum);

    m_sampleUV.resize(N);
    for (int n = 0; n < N; n++)
    {
        m_sampleUV[n] = new DiscreteSampler(columns[n]);
    }
}

void EnvmapLoader::sampleByPdf(float rnd_u, float rnd_v, float& u_ref, float& v_ref) const
{
    int M = m_header.height;
    int N = m_header.width;

    float& u = rnd_u;
    float& v = rnd_v;

//     int count;
// 
//     // TODO: use binary search
// #if 0
//     for (int n = 0; n < N; n++) {
//         count = n;
//         if (m_cdf_u[n] > u) {
//             break;
//         }
//     }
// #else
//     {
//         int begin = 0;
//         int end = N - 1;
//         while (end > begin) {
//             int mid = begin + (end - begin) / 2;
//             float c = m_cdf_u[mid];
//             if (c >= u) {
//                 end = mid;
//             }
//             else {
//                 begin = mid + 1;
//             }
//         }
//         count = begin;
//     }
// #endif
// 
//     //     float u_sample = float(count) / N;
//     //     int u_i = clampi(floor(u_sample * N + 0.5), 0, N - 1);
//     int u_i = count;

    int u_i = m_sampleU->sample(u);

// #if 0
//     for (int m = 0; m < M; m++) {
//         count = m;
//         if (m_cdf_uv[u_i + m * N] > v) {
//             break;
//         }
//     }
// #else
//     {
//         int begin = 0;
//         int end = M - 1;
//         while (end > begin) {
//             int mid = begin + (end - begin) / 2;
//             float c = m_cdf_uv[u_i + mid * N];
//             if (c >= v) {
//                 end = mid;
//             }
//             else {
//                 begin = mid + 1;
//             }
//         }
//         count = begin;
//     }
// #endif
// 
//     //     float v_sample = float(count) / M;
//     int v_i = count;

    int v_i = m_sampleUV[u_i]->sample(v);

    //     u_ref = u_sample + randf() / N;
    //     v_ref = v_sample + randf() / M;
    u_ref = (u_i + randf()) / N;
    v_ref = (v_i + randf()) / M;
}

DirLight EnvmapLoader::sample(float r1, float r2) const
{
    float u, v;
    sampleByPdf(r1, r2, u, v);

    int u_i = clampi(u * m_header.width, 0, m_header.width - 1);
    int v_i = clampi(v * m_header.height, 0, m_header.height - 1);

    int offset = (u_i + v_i * m_header.width) * 3;

    //      float phi = (1 - v) * M_PI;
    float phi = v * NUM_PI;
    float theta = u * 2 * NUM_PI;

    DirLight dl;
    dl.dir[0] = sin(phi) * sin(theta);
    dl.dir[1] = cos(phi);
    dl.dir[2] = sin(phi) * -cos(theta);
    dl.pdf = m_normalizationFactor * m_luminance[u_i + v_i * m_header.width] *
        float(m_header.width) * float(m_header.height) / (2 * NUM_PI * NUM_PI /* * sin(phi) */ /* this one is cancelled */);

    if (!is_pre_warped) {
        dl.pdf /= sin(phi);
    }

    {
        //         if (dl.pdf == 0) {
        //             printf("!!!sampled singular light\n");
        //             continue;
        //         }
    }

    //     dl.pdf = 1;
    //     dl.rad[0] = m_buffer[offset + 0] / dl.pdf;
    //     dl.rad[1] = m_buffer[offset + 1] / dl.pdf;
    //     dl.rad[2] = m_buffer[offset + 2] / dl.pdf;
    dl.rad[0] = m_buffer[offset + 0];
    dl.rad[1] = m_buffer[offset + 1];
    dl.rad[2] = m_buffer[offset + 2];

    return dl;
}

float EnvmapLoader::sampleRadiance(float dir[3], int channel) const {
    float phi = acos(dir[1]);
    float theta = dir_to_theta(dir);

    int m = clampi(int(phi / NUM_PI * m_header.height), 0, m_header.height - 1);
    int n = clampi(int(theta / (NUM_PI * 2) * m_header.width), 0, m_header.width - 1);

    return m_buffer[(n + m * m_header.width) * 3 + channel];
}

void EnvmapLoader::sampleRadiance(float dir[3], float rad[3]) const {
    float phi = acos(dir[1]);
    float theta = dir_to_theta(dir);

    int m = clampi(int(phi / NUM_PI * m_header.height), 0, m_header.height - 1);
    int n = clampi(int(theta / (NUM_PI * 2) * m_header.width), 0, m_header.width - 1);

    rad[0] = m_buffer[(n + m * m_header.width) * 3];
    rad[1] = m_buffer[(n + m * m_header.width) * 3 + 1];
    rad[2] = m_buffer[(n + m * m_header.width) * 3 + 2];
}

// int EnvmapLoader::getLightCount() {
//     return m_lights.size();
// }

float EnvmapLoader::getPdf(float dir[3]) const
{
    float phi = acos(dir[1]);
    float theta = dir_to_theta(dir);
    int u_i = clampi(int(theta / (NUM_PI * 2) * m_header.width), 0, m_header.width - 1);
    int v_i = clampi(int(phi / NUM_PI * m_header.height), 0, m_header.height - 1);

    float pdf = m_normalizationFactor * m_luminance[u_i + v_i * m_header.width] *
        float(m_header.width) * float(m_header.height) / (2 * NUM_PI * NUM_PI /* * sin(phi) */ /* this one is cancelled */);

    if (!is_pre_warped) {
        pdf /= sin(phi);
    }
    return pdf;
}

float EnvmapLoader::getB(int i)
{
    return m_buffer[i * 3 + 2];
}

float EnvmapLoader::getB(int m, int n)
{
    return m_buffer[(n + m * m_header.width) * 3 + 2];
}

float EnvmapLoader::getG(int i)
{
    return m_buffer[i * 3 + 1];
}

float EnvmapLoader::getG(int m, int n)
{
    return m_buffer[(n + m * m_header.width) * 3 + 1];
}

float EnvmapLoader::getR(int i)
{
    return m_buffer[i * 3 + 0];
}

float EnvmapLoader::getR(int m, int n)
{
    return m_buffer[(n + m * m_header.width) * 3 + 0];
}

EnvmapLoader::~EnvmapLoader()
{
    if (m_sampleU)
    {
        delete m_sampleU;
    }
    for (auto p : m_sampleUV)
    {
        if (p)
        {
            delete p;
        }
    }
}

EnvmapLoader::EnvmapLoader()
{
}
