#include <fstream>
#include <iostream>
#include <sstream>
#include "mesh.h"
#include "numeric.h"

using namespace std;

template <typename T>
static void dumpVector(const vector<T>& data, ofstream& ofs)
{
    int len = sizeof(T) * data.size();
    ofs.write((const char *)&len, sizeof(int));
    ofs.write((const char *)&data[0], len);
}

template <typename T>
static void loadVector(vector<T>& data, ifstream& ifs) {
    int len;
    ifs.read((char *)&len, sizeof(int));
    data.resize(len / sizeof(T));
    ifs.read((char *)&data[0], len);
}

CMesh::CMesh(const string& filepath)
{
    m_verts.clear();
    m_faces.clear();
    m_vertex_normals.clear();

    ifstream ifs(filepath.c_str(), ios::binary);

    if (false == ifs.good()) {
        std::cout << "cannot open " << filepath << std::endl;
        return;
    }

    stringstream reader;
    reader << ifs.rdbuf();

    string token;
    while (false == reader.eof()) {
        reader >> token;
        if (token == "v") {
            float a, b, c;
            reader >> a >> b >> c;
            m_verts.push_back(vec3(a, b, c));
        }
        else
            if (token == "f") {
                unsigned int a, b, c;
                reader >> a >> b >> c;
                m_faces.push_back(int3(a - 1, b - 1, c - 1));
            }
    }

    genNormals();
}

void CMesh::genNormals()
{
    m_vertex_normals.resize(m_verts.size());
    m_vertex_normals.assign(m_verts.size(), vec3(0, 0, 0));

    std::vector<vec3> temp_face_normals(m_faces.size(), vec3(0, 0, 0));

    for (int i = 0; i < m_faces.size(); i++)
    {
        int a = m_faces[i].x;
        int b = m_faces[i].y;
        int c = m_faces[i].z;
        vec3 faceNormal = cross(m_verts[b] - m_verts[a], m_verts[c] - m_verts[a]);
        temp_face_normals[i] = normalize(faceNormal);

        m_vertex_normals[a] += faceNormal;
        m_vertex_normals[b] += faceNormal;
        m_vertex_normals[c] += faceNormal;
    }

    for (int i = 0; i < m_vertex_normals.size(); i++)
    {
        m_vertex_normals[i] = normalize(m_vertex_normals[i]);
    }

    m_vertex_alpha.resize(m_verts.size());
    m_vertex_alpha.assign(m_verts.size(), 1); // this is important, since genNormals() may be called multiple times

    for (int i = 0; i < m_faces.size(); i++)
    {
        int a = m_faces[i].x;
        int b = m_faces[i].y;
        int c = m_faces[i].z;

        m_vertex_alpha[a] = f_min(m_vertex_alpha[a],
            dot(temp_face_normals[i], m_vertex_normals[a]));
        m_vertex_alpha[b] = f_min(m_vertex_alpha[b],
            dot(temp_face_normals[i], m_vertex_normals[b]));
        m_vertex_alpha[c] = f_min(m_vertex_alpha[c],
            dot(temp_face_normals[i], m_vertex_normals[c]));
    }

    const float w = 0.03632;
    for (auto& a : m_vertex_alpha)
    {
        assert(a >= -1 && a <= 1);
        a = acos(a) * (1 + w * sq(1 - a));
//         a = f_max(a, 0);
    }
}

void CMesh::dumpBinary(string binfile) const
{
    ofstream ofs(binfile.c_str(), ios::binary | ios::out);
    dumpVector(m_verts, ofs);
    dumpVector(m_faces, ofs);
    dumpVector(m_vertex_normals, ofs);
    ofs.close();
}

void CMesh::loadBinary(string binfile)
{
    ifstream ifs(binfile.c_str(), ios::binary | ios::in);
    loadVector(m_verts, ifs);
    loadVector(m_faces, ifs);
    loadVector(m_vertex_normals, ifs);
    ifs.close();
}

void CMesh::flipSide()
{
    for (int i = 0; i < m_faces.size(); i++)
    {
        m_faces[i] = int3(m_faces[i].x, m_faces[i].z, m_faces[i].y);
    }
}

void CMesh::transform(const glm::mat4& mat)
{
    for (int i = 0; i < m_verts.size(); i++) {
        glm::vec4 v = glm::vec4(m_verts[i].x, m_verts[i].y, m_verts[i].z, 1.0f);
        v = mat * v;
        m_verts[i] = vec3(v.x, v.y, v.z);
    }
}
