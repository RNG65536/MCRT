#include <cassert>
#include "lightpath.h"
#include "triangle.h"

Vert::Vert(vec3 _p, vec3 _n, const TriangleObject* _obj)
    : p(_p), n(_n), obj(_obj)
{
}

Vert::Vert()
{
}

bool Vert::canConnectTo(const vec3& position) const
{
    const vec3 nl = dot(debug_wo, n) > 0 ? n : -n;
    return dot(position - p, nl) > 0;
}

void Vert::setDeltaFlag(char flag)
{
    m_delta_flag = flag;
}

char Vert::isDelta() const
{
    return m_delta_flag;
}

void Vert::setTotalPdfA(double pdf)
{
    m_total_pdfA = pdf;
}

double Vert::getTotalPdfA() const
{
    return m_total_pdfA;
}

void Vert::setOrientingNormal(const vec3& n)
{
    m_orienting_normal = n;
}

const vec3& Vert::getOrientingNormal() const
{
    return m_orienting_normal;
}

void Vert::setThroughput(const vec3& t)
{
    m_throughput = t;
}

const vec3& Vert::getThroughput() const
{
    return m_throughput;
}

Path::Path()
{
    n_verts = 0;
    for (int i = 0; i < MaxEvents; i++)
    {
        x[i].obj = nullptr;
    }
}

Path Path::subpath(int n) const
{
    Path path;

    //         for (int i = 0; i < n; i++)
    //         {
    //             path.insert(x[i]);
    //         }

    for (int i = 0; i < n; i++)
    {
        path.x[i] = x[i];
    }
    path.n_verts = n;

    return path;
}

void Path::resize(int _n)
{
    n_verts = _n;
}

int Path::size() const
{
    return n_verts;
}

void Path::insert(const Vert& v)
{
    assert(n_verts < MaxEvents);
    x[n_verts++] = v;
}

Vert& Path::end()
{
    return x[n_verts - 1];
}

const Vert& Path::end() const
{
    return x[n_verts - 1];
}

Vert& Path::begin()
{
    return x[0];
}

const Vert& Path::begin() const
{
    return x[0];
}

Vert& Path::operator[](int i)
{
    assert(i < n_verts);
    return x[i];
}

const Vert& Path::operator[](int i) const
{
    assert(i < n_verts);
    return x[i];
}

Contribution::Contribution(float x_, float y_, vec3 c_, int s_, int t_)
    : x(x_), y(y_), c(c_), s(s_), t(t_)
{
}

Contribution::Contribution()
{
}

PathContribution::PathContribution()
{
    n = 0;
    sc = 0.0;
}

static Path combinePaths(const Path& eye, const Path& light)
{
    Path SampledPath;

    int es = eye.size();
    int ls = light.size();

    // count of segments
    int path_length = ls + es - 1;

    SampledPath.resize(es + ls);
    for (int i = 0; i < es; i++)
    {
        SampledPath[i] = eye[i];
    }
    for (int i = 0; i < ls; i++)
    {
        SampledPath[path_length - i] = light[i];
    }

    return SampledPath;
}

int CombinedPath::toEyeIndex(int i) const
{
    return path_length - i;
}

int CombinedPath::toLightIndex(int i) const
{
    return path_length - i;
}

const Vert& CombinedPath::operator[](int i) const
{
    return path[i];
}

Vert& CombinedPath::operator[](int i)
{
    return path[i];
}

CombinedPath::CombinedPath(const Path& eye, const Path& light)
{
    path = combinePaths(eye, light);
    num_eye_verts = eye.size();
    num_light_verts = light.size();
    path_length = path.size() - 1;
    assert(path_length + 1 == num_eye_verts + num_light_verts);
}

int CombinedPath::size() const
{
    return path.size();
}

int CombinedPath::lightSize() const
{
    return num_light_verts;
}

int CombinedPath::eyeSize() const
{
    return num_eye_verts;
}

Vert& CombinedPath::eye(int i)
{
    return path[i];
}

Vert& CombinedPath::light(int i)
{
    return path[path_length - i];
}

const Vert& CombinedPath::eye(int i) const
{
    return path[i];
}

const Vert& CombinedPath::light(int i) const
{
    return path[path_length - i];
}

Vert& CombinedPath::lightEnd()
{
    return path[num_eye_verts];
}

const Vert& CombinedPath::lightEnd() const
{
    return path[num_eye_verts];
}

Vert& CombinedPath::lightBegin()
{
    return path[path_length];
}

const Vert& CombinedPath::lightBegin() const
{
    return path[path_length];
}

Vert& CombinedPath::eyeEnd()
{
    return path[num_eye_verts - 1];
}

const Vert& CombinedPath::eyeEnd() const
{
    return path[num_eye_verts - 1];
}

Vert& CombinedPath::eyeBegin()
{
    return path[0];
}

const Vert& CombinedPath::eyeBegin() const
{
    return path[0];
}

int CombinedPath::numNonSpecularEdges() const
{
    int num = 0;
    for (int i = 1; i < eyeSize(); i++)
    {
        if (!eye(i).isDelta() && !eye(i - 1).isDelta())
        {
            num++;
        }
    }
    for (int j = 1; j < lightSize(); j++)
    {
        if (!light(j).isDelta() && !light(j - 1).isDelta())
        {
            ++num;
        }
    }
    if (lightSize() > 0 && !lightEnd().isDelta() && eyeSize() > 0 &&
        !eyeEnd().isDelta())
    {
        num++;
    }
    return num;
}
