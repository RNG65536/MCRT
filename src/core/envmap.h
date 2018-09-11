#ifndef envmap_h__
#define envmap_h__

#include <vector>

// envmap lighting and sampling

struct DirLight
{
    float dir[3];
    float rad[3]; // pre-multiplied with 1/pdf
    float pdf;
};

class DiscreteSampler;

class EnvmapLoader
{
    // u for selecting column, v for selecting row

    struct EnvmapHeader 
    {
        int width, height; // u for width-wise, v for height-wise
    };

    std::vector<float> m_buffer; // row-major, top-down storage
    EnvmapHeader m_header;

//     std::vector<float> m_pdf_u, m_cdf_u; // width-wise [1 * width]              , P(u)   (apriori prob)
//     std::vector<float> m_pdf_uv, m_cdf_uv; // height-wise [height * width]      , P(v|u) (conditional prob)
//     std::vector<DirLight> m_lights; // if full summation is used, then pdf for each sample should not be applied ??
    // at least the bias is smaller than using dl.pdf at small sampling rates ??
    float m_normalizationFactor;
    DiscreteSampler *m_sampleU = nullptr;
    std::vector<DiscreteSampler*> m_sampleUV;

    std::vector<float> m_luminance;
    //     std::vector<float> m_luminance_warped;

    void buildCdf();
    void sampleByPdf(float rnd_u, float rnd_v, float& u_ref, float& v_ref) const;
    void analyze();

public:
    EnvmapLoader();
    EnvmapLoader(const std::string filepath);
    EnvmapLoader(const std::string filepath, int width, int height);
    ~EnvmapLoader();

    void toPPM(const char *filename);

    float getR(int m, int n);
    float getG(int m, int n);
    float getB(int m, int n);
    float getR(int i);
    float getG(int i);
    float getB(int i);

//     void generate_light_samples();
//     DirLight getLightRandom();
//     DirLight getLight(int n);
//     int getLightCount();

    float sampleRadiance(float dir[3], int channel) const;
    void sampleRadiance(float dir[3], float rad[3]) const;
    static float luminance(float r, float g, float b);
    void scale(float x);
    DirLight sample(float r1, float r2) const;
    float getPdf(float dir[3]) const;
};

#endif // envmap_h__
