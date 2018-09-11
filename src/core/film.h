#ifndef film_h__
#define film_h__

#include <vector>
#include "vectors.h"

class FrameBuffer
{
private:
    std::vector<float> m_buffer;
    int m_width, m_height;
    int m_total;

public:
    FrameBuffer();
    FrameBuffer(int w, int h);
    ~FrameBuffer();
    void resize(int w, int h);
    void scale(float s);

    void accumulatePixel(int i, int j, const vec3& c);
    void accumulateBuffer(const FrameBuffer& f);

    void tonemapGamma(float gamma);
    void tonemapReinhard();
    void dumpPPM(const char *filename);
    void dumpHDR(const char *filename);

    vec3 pixel(int i, int j) const;
    int width() const;
    int height() const;
};
#endif // film_h__
