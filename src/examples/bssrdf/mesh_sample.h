#ifndef meshsample_h__
#define meshsample_h__

#include <vector>

struct Position {
    float x, y, z;
    struct Payload { int triangle_idx; float u; float v; } payload;

    Position();
    Position(float x, float y, float z);
    //     Position(float x, float y, float z, const Payload& p) : x(x), y(y), z(z), payload(p) {}
    //     Position(const Position& a) : x(a.x), y(a.y), z(a.z), payload(a.payload) {}

    static float distanceSquared(const Position& a, const Position& b);

};

int poissonSample(const std::vector<Position>& p_in, std::vector<Position>& p_out, float min_dist);
void check_min_distance(std::vector<Position>& pp);

#endif // meshsample_h__
