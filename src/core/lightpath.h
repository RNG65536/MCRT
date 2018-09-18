#pragma once

#include "vectors.h"

#define MinPathLength 1
#define MaxPathLength 11  // 9 // 10
#define MaxEvents (MaxPathLength + 1)

typedef enum { FromEye, FromLight } PathType;
typedef enum { Emitter, Diffuse, Specular, Lens } VertexType;

class TriangleObject;

// path data
class Vert
{
private:
    char   m_delta_flag;
    double m_total_pdfA;
    vec3   m_orienting_normal;
    vec3   m_throughput;

public:
    vec3                  p;  // position
    vec3                  n;  // intrinsic geometry normal
    const TriangleObject* obj = nullptr;

    vec3       weight;  // bsdf_weight returned by bsdf.sample()
    float      pWo;     // forward pdfW, returned by bsdf.sample()
    float      pWi;     // backward pdfW
    float      cos_wo;
    vec3       debug_wi;
    vec3       debug_wo;
    PathType   debug_type;
    VertexType debug_vertex_type;

    Vert();
    Vert(vec3 _p, vec3 _n, const TriangleObject* _obj = nullptr);
    bool        canConnectTo(const vec3& position) const;
    void        setDeltaFlag(char flag);
    char        isDelta() const;
    void        setTotalPdfA(double pdf);
    double      getTotalPdfA() const;
    void        setOrientingNormal(const vec3& n);
    const vec3& getOrientingNormal() const;
    void        setThroughput(const vec3& t);
    const vec3& getThroughput() const;
};

// the light path, represented as a vector of vertices
class Path
{
    Vert x[MaxEvents];
    int  n_verts = 0;  // path length

public:
    Path();

    const Vert& operator[](int i) const;
    Vert&       operator[](int i);
    const Vert& begin() const;
    Vert&       begin();
    const Vert& end() const;
    Vert&       end();
    void        insert(const Vert& v);
    int         size() const;
    void        resize(int _n);
    Path        subpath(int n) const;
};

class CombinedPath
{
    Path path;
    int  num_eye_verts;
    int  num_light_verts;
    int  path_length;

public:
    CombinedPath(const Path& eye, const Path& light);
    Vert&       operator[](int i);
    const Vert& operator[](int i) const;
    int         size() const;
    int         eyeSize() const;
    int         lightSize() const;
    int         toLightIndex(int i) const;
    int         toEyeIndex(int i) const;
    const Vert& eye(int i) const;
    Vert&       eye(int i);
    const Vert& light(int i) const;
    Vert&       light(int i);
    const Vert& eyeBegin() const;
    const Vert& eyeEnd() const;
    const Vert& lightBegin() const;
    const Vert& lightEnd() const;
    Vert&       eyeBegin();
    Vert&       eyeEnd();
    Vert&       lightBegin();
    Vert&       lightEnd();
    int         numNonSpecularEdges() const;
};

// pixel contribution
class Contribution
{
public:
    float x, y;  // pixel coordinates
    vec3  c;     // pixel color
    int   s, t;

    Contribution();
    Contribution(float x_, float y_, vec3 c_, int s_ = -1, int t_ = -1);
};

// pixel contributions of subpaths of different combinations
class PathContribution
{
public:
    Contribution c[MaxEvents * MaxEvents];
    int          n;
    float        sc;  // scalar contribution ??

    PathContribution();
};
