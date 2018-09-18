#pragma once

#include <glm/glm.hpp>
#include <string>
#include <vector>
#include "vectors.h"

// obj mesh loader

class Mesh
{
public:
    Mesh();

    std::vector<vec3>  m_verts;
    std::vector<int3>  m_faces;
    std::vector<vec3>  m_vertex_normals;
    std::vector<float> m_vertex_alpha;

    Mesh(const std::string& filename);
    void genNormals();
    void flipSide();
    void transform(const glm::mat4& mat);
    void dumpBinary(const std::string& binfile) const;
    void loadBinary(const std::string& binfile);
};
