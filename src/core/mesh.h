#ifndef mesh_h__
#define mesh_h__

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "vectors.h"

// obj mesh loader

class CMesh {
public:
    CMesh() {}

    std::vector< vec3 > m_verts;
    std::vector< int3 > m_faces;
    std::vector< vec3 > m_vertex_normals;
    std::vector<float> m_vertex_alpha;

    CMesh(const std::string& filename);
    void genNormals();
    void flipSide();
    void transform(const glm::mat4& mat);
    void dumpBinary(std::string binfile) const;
    void loadBinary(std::string binfile);
};

#endif // mesh_h__
