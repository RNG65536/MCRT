#include "constants.h"
#include "debug_mesh.h"
#include "numeric.h"

typedef float3 vector3f;

DebugMesh::DebugMesh()
{
    // mesh with even faces, split bvh is not better than SAH  object split bvh
    initMesh();

    // for this mesh, split bvh is far better than SAH object split bvh
    // initRandomMesh();
}

void DebugMesh::initMesh()
{
    std::vector<vector3f> vertices;

    vertices.push_back(vector3f(0.184326, -0.566187, 2.238917));
    vertices.push_back(vector3f(-0.690674, -0.566187, 2.238917));
    vertices.push_back(vector3f(0.246826, -0.636499, 2.160792));
    vertices.push_back(vector3f(-0.753174, -0.636499, 2.160792));
    vertices.push_back(vector3f(0.293701, -0.675562, 2.051417));
    vertices.push_back(vector3f(-0.800049, -0.675562, 2.051417));
    vertices.push_back(vector3f(0.098389, -0.753687, 2.090480));
    vertices.push_back(vector3f(-0.604736, -0.753687, 2.090480));
    vertices.push_back(vector3f(0.098389, -0.698999, 2.192042));
    vertices.push_back(vector3f(-0.604736, -0.698999, 2.192042));
    vertices.push_back(vector3f(0.098389, -0.597437, 2.254542));
    vertices.push_back(vector3f(-0.604736, -0.597437, 2.254542));
    vertices.push_back(vector3f(0.020264, -0.566187, 2.270167));
    vertices.push_back(vector3f(-0.526611, -0.566187, 2.270167));
    vertices.push_back(vector3f(-0.050049, -0.636499, 2.215480));
    vertices.push_back(vector3f(-0.456299, -0.636499, 2.215480));
    vertices.push_back(vector3f(-0.096924, -0.675562, 2.121730));
    vertices.push_back(vector3f(-0.409424, -0.675562, 2.121730));
    vertices.push_back(vector3f(-0.175049, -0.488062, 2.129542));
    vertices.push_back(vector3f(-0.331299, -0.488062, 2.129542));
    vertices.push_back(vector3f(-0.112549, -0.488062, 2.215480));
    vertices.push_back(vector3f(-0.393799, -0.488062, 2.215480));
    vertices.push_back(vector3f(-0.010986, -0.488062, 2.270167));
    vertices.push_back(vector3f(-0.495361, -0.488062, 2.270167));
    vertices.push_back(vector3f(0.020264, -0.402124, 2.270167));
    vertices.push_back(vector3f(-0.526611, -0.402124, 2.270167));
    vertices.push_back(vector3f(-0.050049, -0.339624, 2.215480));
    vertices.push_back(vector3f(-0.456299, -0.339624, 2.215480));
    vertices.push_back(vector3f(-0.096924, -0.292749, 2.121730));
    vertices.push_back(vector3f(-0.409424, -0.292749, 2.121730));
    vertices.push_back(vector3f(0.098389, -0.214624, 2.090480));
    vertices.push_back(vector3f(-0.604736, -0.214624, 2.090480));
    vertices.push_back(vector3f(0.098389, -0.277124, 2.192042));
    vertices.push_back(vector3f(-0.604736, -0.277124, 2.192042));
    vertices.push_back(vector3f(0.098389, -0.370874, 2.254542));
    vertices.push_back(vector3f(-0.604736, -0.370874, 2.254542));
    vertices.push_back(vector3f(0.184326, -0.402124, 2.238917));
    vertices.push_back(vector3f(-0.690674, -0.402124, 2.238917));
    vertices.push_back(vector3f(0.246826, -0.339624, 2.160792));
    vertices.push_back(vector3f(-0.753174, -0.339624, 2.160792));
    vertices.push_back(vector3f(0.293701, -0.292749, 2.051417));
    vertices.push_back(vector3f(-0.800049, -0.292749, 2.051417));
    vertices.push_back(vector3f(0.371826, -0.488062, 2.035792));
    vertices.push_back(vector3f(-0.878174, -0.488062, 2.035792));
    vertices.push_back(vector3f(0.309326, -0.488062, 2.145167));
    vertices.push_back(vector3f(-0.815674, -0.488062, 2.145167));
    vertices.push_back(vector3f(0.215576, -0.488062, 2.231105));
    vertices.push_back(vector3f(-0.721924, -0.488062, 2.231105));
    vertices.push_back(vector3f(0.223389, -0.488062, 2.246730));
    vertices.push_back(vector3f(-0.729736, -0.488062, 2.246730));
    vertices.push_back(vector3f(0.192139, -0.394312, 2.254542));
    vertices.push_back(vector3f(-0.698486, -0.394312, 2.254542));
    vertices.push_back(vector3f(0.098389, -0.355249, 2.277980));
    vertices.push_back(vector3f(-0.604736, -0.355249, 2.277980));
    vertices.push_back(vector3f(0.012451, -0.394312, 2.293605));
    vertices.push_back(vector3f(-0.518799, -0.394312, 2.293605));
    vertices.push_back(vector3f(-0.026611, -0.488062, 2.293605));
    vertices.push_back(vector3f(-0.479736, -0.488062, 2.293605));
    vertices.push_back(vector3f(0.012451, -0.573999, 2.293605));
    vertices.push_back(vector3f(-0.518799, -0.573999, 2.293605));
    vertices.push_back(vector3f(0.098389, -0.488062, 2.301417));
    vertices.push_back(vector3f(-0.604736, -0.488062, 2.301417));
    vertices.push_back(vector3f(0.098389, -0.613062, 2.277980));
    vertices.push_back(vector3f(-0.604736, -0.613062, 2.277980));
    vertices.push_back(vector3f(0.192139, -0.573999, 2.254542));
    vertices.push_back(vector3f(-0.698486, -0.573999, 2.254542));
    vertices.push_back(vector3f(-0.253174, -0.300562, 2.215480));
    vertices.push_back(vector3f(-0.253174, -0.378687, 2.293605));
    vertices.push_back(vector3f(-0.253174, -1.409937, 2.207667));
    vertices.push_back(vector3f(-0.253174, -1.050562, 2.254542));
    vertices.push_back(vector3f(-0.253174, -0.917749, 2.270167));
    vertices.push_back(vector3f(-0.253174, -1.503687, 2.192042));
    vertices.push_back(vector3f(-0.253174, -0.323999, 2.074855));
    vertices.push_back(vector3f(-0.253174, -0.159937, 2.043605));
    vertices.push_back(vector3f(-0.253174, 0.168188, 0.926417));
    vertices.push_back(vector3f(-0.253174, -0.167749, 0.621730));
    vertices.push_back(vector3f(-0.253174, -0.659937, 0.645167));
    vertices.push_back(vector3f(-0.253174, -1.113062, 1.121730));
    vertices.push_back(vector3f(-0.050049, -0.917749, 2.035792));
    vertices.push_back(vector3f(-0.456299, -0.917749, 2.035792));
    vertices.push_back(vector3f(0.059326, -1.167749, 2.043605));
    vertices.push_back(vector3f(-0.565674, -1.167749, 2.043605));
    vertices.push_back(vector3f(0.098389, -1.425562, 2.043605));
    vertices.push_back(vector3f(-0.604736, -1.425562, 2.043605));
    vertices.push_back(vector3f(0.114014, -1.620874, 2.004542));
    vertices.push_back(vector3f(-0.620361, -1.620874, 2.004542));
    vertices.push_back(vector3f(0.074951, -1.675562, 1.996730));
    vertices.push_back(vector3f(-0.581299, -1.675562, 1.996730));
    vertices.push_back(vector3f(-0.073486, -1.698999, 2.027980));
    vertices.push_back(vector3f(-0.432861, -1.698999, 2.027980));
    vertices.push_back(vector3f(-0.253174, -1.714624, 2.051417));
    vertices.push_back(vector3f(0.184326, -0.870874, 2.004542));
    vertices.push_back(vector3f(-0.690674, -0.870874, 2.004542));
    vertices.push_back(vector3f(0.379639, -0.769312, 2.012355));
    vertices.push_back(vector3f(-0.885986, -0.769312, 2.012355));
    vertices.push_back(vector3f(0.574951, -0.581812, 1.918605));
    vertices.push_back(vector3f(-1.081299, -0.581812, 1.918605));
    vertices.push_back(vector3f(0.606201, -0.300562, 2.067042));
    vertices.push_back(vector3f(-1.112549, -0.300562, 2.067042));
    vertices.push_back(vector3f(0.457764, -0.245874, 2.098292));
    vertices.push_back(vector3f(-0.964111, -0.245874, 2.098292));
    vertices.push_back(vector3f(0.239014, -0.128687, 2.160792));
    vertices.push_back(vector3f(-0.745361, -0.128687, 2.160792));
    vertices.push_back(vector3f(0.067139, 0.027563, 2.207667));
    vertices.push_back(vector3f(-0.573486, 0.027563, 2.207667));
    vertices.push_back(vector3f(-0.096924, -0.011499, 2.231105));
    vertices.push_back(vector3f(-0.409424, -0.011499, 2.231105));
    vertices.push_back(vector3f(-0.190674, -0.238062, 2.223292));
    vertices.push_back(vector3f(-0.315674, -0.238062, 2.223292));
    vertices.push_back(vector3f(-0.089111, -0.316187, 2.246730));
    vertices.push_back(vector3f(-0.417236, -0.316187, 2.246730));
    vertices.push_back(vector3f(-0.128174, -0.425562, 2.238917));
    vertices.push_back(vector3f(-0.378174, -0.425562, 2.238917));
    vertices.push_back(vector3f(-0.050049, -0.636499, 2.215480));
    vertices.push_back(vector3f(-0.456299, -0.636499, 2.215480));
    vertices.push_back(vector3f(0.121826, -0.714624, 2.176417));
    vertices.push_back(vector3f(-0.628174, -0.714624, 2.176417));
    vertices.push_back(vector3f(0.239014, -0.667749, 2.145167));
    vertices.push_back(vector3f(-0.745361, -0.667749, 2.145167));
    vertices.push_back(vector3f(0.371826, -0.542749, 2.121730));
    vertices.push_back(vector3f(-0.878174, -0.542749, 2.121730));
    vertices.push_back(vector3f(0.387451, -0.433374, 2.121730));
    vertices.push_back(vector3f(-0.893799, -0.433374, 2.121730));
    vertices.push_back(vector3f(0.348389, -0.355249, 2.137355));
    vertices.push_back(vector3f(-0.854736, -0.355249, 2.137355));
    vertices.push_back(vector3f(0.176514, -0.292749, 2.192042));
    vertices.push_back(vector3f(-0.682861, -0.292749, 2.192042));
    vertices.push_back(vector3f(-0.003174, -0.261499, 2.231105));
    vertices.push_back(vector3f(-0.503174, -0.261499, 2.231105));
    vertices.push_back(vector3f(-0.253174, -1.495874, 2.207667));
    vertices.push_back(vector3f(-0.143799, -1.448999, 2.207667));
    vertices.push_back(vector3f(-0.362549, -1.448999, 2.207667));
    vertices.push_back(vector3f(-0.135986, -1.566187, 2.184230));
    vertices.push_back(vector3f(-0.370361, -1.566187, 2.184230));
    vertices.push_back(vector3f(-0.190674, -1.613062, 2.168605));
    vertices.push_back(vector3f(-0.315674, -1.613062, 2.168605));
    vertices.push_back(vector3f(-0.253174, -1.620874, 2.160792));
    vertices.push_back(vector3f(-0.253174, -0.925562, 2.223292));
    vertices.push_back(vector3f(-0.253174, -0.870874, 2.215480));
    vertices.push_back(vector3f(-0.151611, -0.878687, 2.215480));
    vertices.push_back(vector3f(-0.354736, -0.878687, 2.215480));
    vertices.push_back(vector3f(-0.128174, -0.956812, 2.223292));
    vertices.push_back(vector3f(-0.378174, -0.956812, 2.223292));
    vertices.push_back(vector3f(-0.167236, -1.019312, 2.215480));
    vertices.push_back(vector3f(-0.339111, -1.019312, 2.215480));
    vertices.push_back(vector3f(0.145264, -0.777124, 2.145167));
    vertices.push_back(vector3f(-0.651611, -0.777124, 2.145167));
    vertices.push_back(vector3f(0.364014, -0.675562, 2.098292));
    vertices.push_back(vector3f(-0.870361, -0.675562, 2.098292));
    vertices.push_back(vector3f(0.473389, -0.527124, 2.074855));
    vertices.push_back(vector3f(-0.979736, -0.527124, 2.074855));
    vertices.push_back(vector3f(0.489014, -0.355249, 2.129542));
    vertices.push_back(vector3f(-0.995361, -0.355249, 2.129542));
    vertices.push_back(vector3f(0.434326, -0.316187, 2.199855));
    vertices.push_back(vector3f(-0.940674, -0.316187, 2.199855));
    vertices.push_back(vector3f(0.184326, -0.183374, 2.270167));
    vertices.push_back(vector3f(-0.690674, -0.183374, 2.270167));
    vertices.push_back(vector3f(0.059326, -0.089624, 2.309230));
    vertices.push_back(vector3f(-0.565674, -0.089624, 2.309230));
    vertices.push_back(vector3f(-0.050049, -0.113062, 2.324855));
    vertices.push_back(vector3f(-0.456299, -0.113062, 2.324855));
    vertices.push_back(vector3f(-0.151611, -0.300562, 2.317042));
    vertices.push_back(vector3f(-0.354736, -0.300562, 2.317042));
    vertices.push_back(vector3f(-0.128174, -0.831812, 2.285792));
    vertices.push_back(vector3f(-0.378174, -0.831812, 2.285792));
    vertices.push_back(vector3f(-0.042236, -1.175562, 2.184230));
    vertices.push_back(vector3f(-0.464111, -1.175562, 2.184230));
    vertices.push_back(vector3f(-0.003174, -1.433374, 2.160792));
    vertices.push_back(vector3f(-0.503174, -1.433374, 2.160792));
    vertices.push_back(vector3f(0.012451, -1.550562, 2.137355));
    vertices.push_back(vector3f(-0.518799, -1.550562, 2.137355));
    vertices.push_back(vector3f(-0.018799, -1.644312, 2.106105));
    vertices.push_back(vector3f(-0.487549, -1.644312, 2.106105));
    vertices.push_back(vector3f(-0.089111, -1.659937, 2.106105));
    vertices.push_back(vector3f(-0.417236, -1.659937, 2.106105));
    vertices.push_back(vector3f(-0.253174, -1.675562, 2.113917));
    vertices.push_back(vector3f(-0.253174, -0.683374, 2.199855));
    vertices.push_back(vector3f(-0.253174, -0.519312, 2.238917));
    vertices.push_back(vector3f(0.074951, -0.253687, 2.215480));
    vertices.push_back(vector3f(-0.581299, -0.253687, 2.215480));
    vertices.push_back(vector3f(-0.089111, -0.589624, 2.223292));
    vertices.push_back(vector3f(-0.417236, -0.589624, 2.223292));
    vertices.push_back(vector3f(-0.120361, -0.519312, 2.231105));
    vertices.push_back(vector3f(-0.385986, -0.519312, 2.231105));
    vertices.push_back(vector3f(-0.135986, -1.417749, 2.207667));
    vertices.push_back(vector3f(-0.370361, -1.417749, 2.207667));
    vertices.push_back(vector3f(-0.175049, -1.175562, 2.223292));
    vertices.push_back(vector3f(-0.331299, -1.175562, 2.223292));
    vertices.push_back(vector3f(-0.253174, -1.175562, 2.223292));
    vertices.push_back(vector3f(-0.253174, -1.058374, 2.215480));
    vertices.push_back(vector3f(-0.159424, -1.003687, 2.254542));
    vertices.push_back(vector3f(-0.346924, -1.003687, 2.254542));
    vertices.push_back(vector3f(-0.120361, -0.956812, 2.270167));
    vertices.push_back(vector3f(-0.385986, -0.956812, 2.270167));
    vertices.push_back(vector3f(-0.143799, -0.863062, 2.254542));
    vertices.push_back(vector3f(-0.362549, -0.863062, 2.254542));
    vertices.push_back(vector3f(-0.214111, -0.855249, 2.254542));
    vertices.push_back(vector3f(-0.292236, -0.855249, 2.254542));
    vertices.push_back(vector3f(-0.253174, -0.933374, 2.301417));
    vertices.push_back(vector3f(-0.206299, -0.878687, 2.285792));
    vertices.push_back(vector3f(-0.300049, -0.878687, 2.285792));
    vertices.push_back(vector3f(-0.159424, -0.886499, 2.285792));
    vertices.push_back(vector3f(-0.346924, -0.886499, 2.285792));
    vertices.push_back(vector3f(-0.143799, -0.956812, 2.301417));
    vertices.push_back(vector3f(-0.362549, -0.956812, 2.301417));
    vertices.push_back(vector3f(-0.175049, -0.980249, 2.277980));
    vertices.push_back(vector3f(-0.331299, -0.980249, 2.277980));
    vertices.push_back(vector3f(-0.253174, -1.019312, 2.277980));
    vertices.push_back(vector3f(0.004639, -1.042749, 2.027980));
    vertices.push_back(vector3f(-0.510986, -1.042749, 2.027980));
    vertices.push_back(vector3f(-0.089111, -0.972437, 2.184230));
    vertices.push_back(vector3f(-0.417236, -0.972437, 2.184230));
    vertices.push_back(vector3f(-0.073486, -1.042749, 2.184230));
    vertices.push_back(vector3f(-0.432861, -1.042749, 2.184230));
    vertices.push_back(vector3f(-0.018799, -0.980249, 2.027980));
    vertices.push_back(vector3f(-0.487549, -0.980249, 2.027980));
    vertices.push_back(vector3f(-0.253174, -1.605249, 2.160792));
    vertices.push_back(vector3f(-0.206299, -1.597437, 2.160792));
    vertices.push_back(vector3f(-0.300049, -1.597437, 2.160792));
    vertices.push_back(vector3f(-0.159424, -1.550562, 2.184230));
    vertices.push_back(vector3f(-0.346924, -1.550562, 2.184230));
    vertices.push_back(vector3f(-0.159424, -1.472437, 2.199855));
    vertices.push_back(vector3f(-0.346924, -1.472437, 2.199855));
    vertices.push_back(vector3f(-0.253174, -1.511499, 2.129542));
    vertices.push_back(vector3f(-0.159424, -1.480249, 2.137355));
    vertices.push_back(vector3f(-0.346924, -1.480249, 2.137355));
    vertices.push_back(vector3f(-0.159424, -1.542749, 2.113917));
    vertices.push_back(vector3f(-0.346924, -1.542749, 2.113917));
    vertices.push_back(vector3f(-0.206299, -1.581812, 2.106105));
    vertices.push_back(vector3f(-0.300049, -1.581812, 2.106105));
    vertices.push_back(vector3f(-0.253174, -1.589624, 2.106105));
    vertices.push_back(vector3f(-0.081299, -0.511499, 2.254542));
    vertices.push_back(vector3f(-0.425049, -0.511499, 2.254542));
    vertices.push_back(vector3f(-0.065674, -0.573999, 2.246730));
    vertices.push_back(vector3f(-0.440674, -0.573999, 2.246730));
    vertices.push_back(vector3f(0.082764, -0.300562, 2.231105));
    vertices.push_back(vector3f(-0.589111, -0.300562, 2.231105));
    vertices.push_back(vector3f(0.020264, -0.308374, 2.246730));
    vertices.push_back(vector3f(-0.526611, -0.308374, 2.246730));
    vertices.push_back(vector3f(0.168701, -0.331812, 2.246730));
    vertices.push_back(vector3f(-0.675049, -0.331812, 2.246730));
    vertices.push_back(vector3f(0.309326, -0.378687, 2.168605));
    vertices.push_back(vector3f(-0.815674, -0.378687, 2.168605));
    vertices.push_back(vector3f(0.332764, -0.441187, 2.160792));
    vertices.push_back(vector3f(-0.839111, -0.441187, 2.160792));
    vertices.push_back(vector3f(0.324951, -0.534937, 2.152980));
    vertices.push_back(vector3f(-0.831299, -0.534937, 2.152980));
    vertices.push_back(vector3f(0.223389, -0.628687, 2.192042));
    vertices.push_back(vector3f(-0.729736, -0.628687, 2.192042));
    vertices.push_back(vector3f(0.121826, -0.667749, 2.215480));
    vertices.push_back(vector3f(-0.628174, -0.667749, 2.215480));
    vertices.push_back(vector3f(-0.026611, -0.620874, 2.254542));
    vertices.push_back(vector3f(-0.479736, -0.620874, 2.254542));
    vertices.push_back(vector3f(-0.073486, -0.433374, 2.254542));
    vertices.push_back(vector3f(-0.432861, -0.433374, 2.254542));
    vertices.push_back(vector3f(-0.042236, -0.355249, 2.254542));
    vertices.push_back(vector3f(-0.464111, -0.355249, 2.254542));
    vertices.push_back(vector3f(-0.018799, -0.370874, 2.231105));
    vertices.push_back(vector3f(-0.487549, -0.370874, 2.231105));
    vertices.push_back(vector3f(-0.057861, -0.433374, 2.231105));
    vertices.push_back(vector3f(-0.448486, -0.433374, 2.231105));
    vertices.push_back(vector3f(-0.010986, -0.605249, 2.231105));
    vertices.push_back(vector3f(-0.495361, -0.605249, 2.231105));
    vertices.push_back(vector3f(0.121826, -0.644312, 2.199855));
    vertices.push_back(vector3f(-0.628174, -0.644312, 2.199855));
    vertices.push_back(vector3f(0.207764, -0.613062, 2.176417));
    vertices.push_back(vector3f(-0.714111, -0.613062, 2.176417));
    vertices.push_back(vector3f(0.293701, -0.519312, 2.145167));
    vertices.push_back(vector3f(-0.800049, -0.519312, 2.145167));
    vertices.push_back(vector3f(0.301514, -0.448999, 2.145167));
    vertices.push_back(vector3f(-0.807861, -0.448999, 2.145167));
    vertices.push_back(vector3f(0.278076, -0.394312, 2.152980));
    vertices.push_back(vector3f(-0.784424, -0.394312, 2.152980));
    vertices.push_back(vector3f(0.160889, -0.339624, 2.223292));
    vertices.push_back(vector3f(-0.667236, -0.339624, 2.223292));
    vertices.push_back(vector3f(0.028076, -0.331812, 2.238917));
    vertices.push_back(vector3f(-0.534424, -0.331812, 2.238917));
    vertices.push_back(vector3f(0.082764, -0.323999, 2.223292));
    vertices.push_back(vector3f(-0.589111, -0.323999, 2.223292));
    vertices.push_back(vector3f(-0.050049, -0.558374, 2.223292));
    vertices.push_back(vector3f(-0.456299, -0.558374, 2.223292));
    vertices.push_back(vector3f(-0.057861, -0.503687, 2.223292));
    vertices.push_back(vector3f(-0.448486, -0.503687, 2.223292));
    vertices.push_back(vector3f(-0.143799, -0.269312, 2.082667));
    vertices.push_back(vector3f(-0.362549, -0.269312, 2.082667));
    vertices.push_back(vector3f(-0.057861, -0.066187, 2.090480));
    vertices.push_back(vector3f(-0.448486, -0.066187, 2.090480));
    vertices.push_back(vector3f(0.082764, -0.042749, 2.067042));
    vertices.push_back(vector3f(-0.589111, -0.042749, 2.067042));
    vertices.push_back(vector3f(0.231201, -0.175562, 2.027980));
    vertices.push_back(vector3f(-0.737549, -0.175562, 2.027980));
    vertices.push_back(vector3f(0.426514, -0.277124, 1.965480));
    vertices.push_back(vector3f(-0.932861, -0.277124, 1.965480));
    vertices.push_back(vector3f(0.543701, -0.323999, 1.934230));
    vertices.push_back(vector3f(-1.050049, -0.323999, 1.934230));
    vertices.push_back(vector3f(0.520264, -0.566187, 1.848292));
    vertices.push_back(vector3f(-1.026611, -0.566187, 1.848292));
    vertices.push_back(vector3f(0.348389, -0.730249, 1.887355));
    vertices.push_back(vector3f(-0.854736, -0.730249, 1.887355));
    vertices.push_back(vector3f(0.184326, -0.823999, 1.942042));
    vertices.push_back(vector3f(-0.690674, -0.823999, 1.942042));
    vertices.push_back(vector3f(-0.253174, 0.168188, 1.762355));
    vertices.push_back(vector3f(-0.253174, 0.254126, 1.395167));
    vertices.push_back(vector3f(-0.253174, -0.925562, 0.801417));
    vertices.push_back(vector3f(-0.253174, -1.191187, 1.660792));
    vertices.push_back(vector3f(-0.253174, -1.706812, 1.934230));
    vertices.push_back(vector3f(-0.253174, -1.534937, 1.817042));
    vertices.push_back(vector3f(-0.253174, -1.300562, 1.793605));
    vertices.push_back(vector3f(-0.253174, -1.214624, 1.754542));
    vertices.push_back(vector3f(0.598389, -0.495874, 1.527980));
    vertices.push_back(vector3f(-1.104736, -0.495874, 1.527980));
    vertices.push_back(vector3f(0.606201, -0.409937, 1.426417));
    vertices.push_back(vector3f(-1.112549, -0.409937, 1.426417));
    vertices.push_back(vector3f(0.520264, -0.464624, 1.035792));
    vertices.push_back(vector3f(-1.026611, -0.464624, 1.035792));
    vertices.push_back(vector3f(0.207764, -0.292749, 0.770167));
    vertices.push_back(vector3f(-0.714111, -0.292749, 0.770167));
    vertices.push_back(vector3f(0.481201, -0.777124, 1.543605));
    vertices.push_back(vector3f(-0.987549, -0.777124, 1.543605));
    vertices.push_back(vector3f(0.340576, -0.855249, 1.309230));
    vertices.push_back(vector3f(-0.846924, -0.855249, 1.309230));
    vertices.push_back(vector3f(0.387451, -0.738062, 1.043605));
    vertices.push_back(vector3f(-0.893799, -0.738062, 1.043605));
    vertices.push_back(vector3f(0.082764, -0.675562, 0.809230));
    vertices.push_back(vector3f(-0.589111, -0.675562, 0.809230));
    vertices.push_back(vector3f(-0.018799, -1.081812, 1.879542));
    vertices.push_back(vector3f(-0.487549, -1.081812, 1.879542));
    vertices.push_back(vector3f(-0.073486, -1.144312, 1.731105));
    vertices.push_back(vector3f(-0.432861, -1.144312, 1.731105));
    vertices.push_back(vector3f(0.035889, -1.441187, 1.856105));
    vertices.push_back(vector3f(-0.542236, -1.441187, 1.856105));
    vertices.push_back(vector3f(-0.003174, -1.230249, 1.863917));
    vertices.push_back(vector3f(-0.503174, -1.230249, 1.863917));
    vertices.push_back(vector3f(0.074951, -1.644312, 1.871730));
    vertices.push_back(vector3f(-0.581299, -1.644312, 1.871730));
    vertices.push_back(vector3f(-0.112549, -1.488062, 1.840480));
    vertices.push_back(vector3f(-0.393799, -1.488062, 1.840480));
    vertices.push_back(vector3f(-0.128174, -1.269312, 1.832667));
    vertices.push_back(vector3f(-0.378174, -1.269312, 1.832667));
    vertices.push_back(vector3f(-0.089111, -1.675562, 1.910792));
    vertices.push_back(vector3f(-0.417236, -1.675562, 1.910792));
    vertices.push_back(vector3f(-0.034424, -1.011499, 1.902980));
    vertices.push_back(vector3f(-0.471924, -1.011499, 1.902980));
    vertices.push_back(vector3f(-0.042236, -0.956812, 1.942042));
    vertices.push_back(vector3f(-0.464111, -0.956812, 1.942042));
    vertices.push_back(vector3f(-0.050049, -0.902124, 1.973292));
    vertices.push_back(vector3f(-0.456299, -0.902124, 1.973292));
    vertices.push_back(vector3f(-0.042236, -1.120874, 1.637355));
    vertices.push_back(vector3f(-0.464111, -1.120874, 1.637355));
    vertices.push_back(vector3f(0.043701, -1.042749, 1.207667));
    vertices.push_back(vector3f(-0.550049, -1.042749, 1.207667));
    vertices.push_back(vector3f(0.090576, -0.878687, 0.934230));
    vertices.push_back(vector3f(-0.596924, -0.878687, 0.934230));
    vertices.push_back(vector3f(0.199951, 0.136938, 1.090480));
    vertices.push_back(vector3f(-0.706299, 0.136938, 1.090480));
    vertices.push_back(vector3f(0.199951, 0.199438, 1.402980));
    vertices.push_back(vector3f(-0.706299, 0.199438, 1.402980));
    vertices.push_back(vector3f(0.199951, 0.121313, 1.707667));
    vertices.push_back(vector3f(-0.706299, 0.121313, 1.707667));
    vertices.push_back(vector3f(0.207764, -0.206812, 1.902980));
    vertices.push_back(vector3f(-0.714111, -0.206812, 1.902980));
    vertices.push_back(vector3f(0.473389, -0.323999, 1.809230));
    vertices.push_back(vector3f(-0.979736, -0.323999, 1.809230));
    vertices.push_back(vector3f(0.379639, -0.277124, 1.754542));
    vertices.push_back(vector3f(-0.885986, -0.277124, 1.754542));
    vertices.push_back(vector3f(0.387451, -0.027124, 1.527980));
    vertices.push_back(vector3f(-0.893799, -0.027124, 1.527980));
    vertices.push_back(vector3f(0.543701, -0.167749, 1.598292));
    vertices.push_back(vector3f(-1.050049, -0.167749, 1.598292));
    vertices.push_back(vector3f(0.543701, -0.113062, 1.356105));
    vertices.push_back(vector3f(-1.050049, -0.113062, 1.356105));
    vertices.push_back(vector3f(0.387451, 0.019751, 1.277980));
    vertices.push_back(vector3f(-0.893799, 0.019751, 1.277980));
    vertices.push_back(vector3f(0.387451, -0.050562, 1.027980));
    vertices.push_back(vector3f(-0.893799, -0.050562, 1.027980));
    vertices.push_back(vector3f(0.543701, -0.191187, 1.113917));
    vertices.push_back(vector3f(-1.050049, -0.191187, 1.113917));
    vertices.push_back(vector3f(0.364014, -0.402124, 0.887355));
    vertices.push_back(vector3f(-0.870361, -0.402124, 0.887355));
    vertices.push_back(vector3f(0.231201, -0.706812, 0.926417));
    vertices.push_back(vector3f(-0.737549, -0.706812, 0.926417));
    vertices.push_back(vector3f(0.567139, -0.402124, 1.270167));
    vertices.push_back(vector3f(-1.073486, -0.402124, 1.270167));
    vertices.push_back(vector3f(0.153076, -0.902124, 1.621730));
    vertices.push_back(vector3f(-0.659424, -0.902124, 1.621730));
    vertices.push_back(vector3f(0.176514, -0.925562, 1.262355));
    vertices.push_back(vector3f(-0.682861, -0.925562, 1.262355));
    vertices.push_back(vector3f(0.637451, -0.323999, 1.238917));
    vertices.push_back(vector3f(-1.143799, -0.323999, 1.238917));
    vertices.push_back(vector3f(0.520264, -0.870874, 1.348292));
    vertices.push_back(vector3f(-1.026611, -0.870874, 1.348292));
    vertices.push_back(vector3f(0.785889, -0.831812, 1.145167));
    vertices.push_back(vector3f(-1.292236, -0.831812, 1.145167));
    vertices.push_back(vector3f(1.028076, -0.675562, 1.043605));
    vertices.push_back(vector3f(-1.534424, -0.675562, 1.043605));
    vertices.push_back(vector3f(1.098389, -0.409937, 1.051417));
    vertices.push_back(vector3f(-1.604736, -0.409937, 1.051417));
    vertices.push_back(vector3f(0.981201, -0.222437, 1.051417));
    vertices.push_back(vector3f(-1.487549, -0.222437, 1.051417));
    vertices.push_back(vector3f(0.770264, -0.253687, 1.160792));
    vertices.push_back(vector3f(-1.276611, -0.253687, 1.160792));
    vertices.push_back(vector3f(0.762451, -0.316187, 1.184230));
    vertices.push_back(vector3f(-1.268799, -0.316187, 1.184230));
    vertices.push_back(vector3f(0.934326, -0.292749, 1.082667));
    vertices.push_back(vector3f(-1.440674, -0.292749, 1.082667));
    vertices.push_back(vector3f(1.012451, -0.441187, 1.067042));
    vertices.push_back(vector3f(-1.518799, -0.441187, 1.067042));
    vertices.push_back(vector3f(0.957764, -0.652124, 1.067042));
    vertices.push_back(vector3f(-1.464111, -0.652124, 1.067042));
    vertices.push_back(vector3f(0.778076, -0.769312, 1.168605));
    vertices.push_back(vector3f(-1.284424, -0.769312, 1.168605));
    vertices.push_back(vector3f(0.574951, -0.800562, 1.340480));
    vertices.push_back(vector3f(-1.081299, -0.800562, 1.340480));
    vertices.push_back(vector3f(0.668701, -0.370874, 1.254542));
    vertices.push_back(vector3f(-1.175049, -0.370874, 1.254542));
    vertices.push_back(vector3f(0.692139, -0.425562, 1.184230));
    vertices.push_back(vector3f(-1.198486, -0.425562, 1.184230));
    vertices.push_back(vector3f(0.629639, -0.753687, 1.262355));
    vertices.push_back(vector3f(-1.135986, -0.753687, 1.262355));
    vertices.push_back(vector3f(0.785889, -0.730249, 1.106105));
    vertices.push_back(vector3f(-1.292236, -0.730249, 1.106105));
    vertices.push_back(vector3f(0.934326, -0.636499, 1.027980));
    vertices.push_back(vector3f(-1.440674, -0.636499, 1.027980));
    vertices.push_back(vector3f(0.981201, -0.480249, 1.027980));
    vertices.push_back(vector3f(-1.487549, -0.480249, 1.027980));
    vertices.push_back(vector3f(0.918701, -0.370874, 1.035792));
    vertices.push_back(vector3f(-1.425049, -0.370874, 1.035792));
    vertices.push_back(vector3f(0.770264, -0.386499, 1.113917));
    vertices.push_back(vector3f(-1.276611, -0.386499, 1.113917));
    vertices.push_back(vector3f(0.590576, -0.441187, 1.262355));
    vertices.push_back(vector3f(-1.096924, -0.441187, 1.262355));
    vertices.push_back(vector3f(0.582764, -0.558374, 1.199855));
    vertices.push_back(vector3f(-1.089111, -0.558374, 1.199855));
    vertices.push_back(vector3f(0.504639, -0.636499, 1.199855));
    vertices.push_back(vector3f(-1.010986, -0.636499, 1.199855));
    vertices.push_back(vector3f(0.567139, -0.644312, 1.199855));
    vertices.push_back(vector3f(-1.073486, -0.644312, 1.199855));
    vertices.push_back(vector3f(0.590576, -0.714624, 1.199855));
    vertices.push_back(vector3f(-1.096924, -0.714624, 1.199855));
    vertices.push_back(vector3f(0.559326, -0.745874, 1.199855));
    vertices.push_back(vector3f(-1.065674, -0.745874, 1.199855));
    vertices.push_back(vector3f(0.473389, -0.730249, 1.402980));
    vertices.push_back(vector3f(-0.979736, -0.730249, 1.402980));
    vertices.push_back(vector3f(0.465576, -0.753687, 1.301417));
    vertices.push_back(vector3f(-0.971924, -0.753687, 1.301417));
    vertices.push_back(vector3f(0.465576, -0.691187, 1.285792));
    vertices.push_back(vector3f(-0.971924, -0.691187, 1.285792));
    vertices.push_back(vector3f(0.543701, -0.527124, 1.262355));
    vertices.push_back(vector3f(-1.050049, -0.527124, 1.262355));
    vertices.push_back(vector3f(0.637451, -0.488062, 1.207667));
    vertices.push_back(vector3f(-1.143799, -0.488062, 1.207667));
    vertices.push_back(vector3f(0.637451, -0.495874, 1.152980));
    vertices.push_back(vector3f(-1.143799, -0.495874, 1.152980));
    vertices.push_back(vector3f(0.559326, -0.745874, 1.152980));
    vertices.push_back(vector3f(-1.065674, -0.745874, 1.152980));
    vertices.push_back(vector3f(0.598389, -0.714624, 1.152980));
    vertices.push_back(vector3f(-1.104736, -0.714624, 1.152980));
    vertices.push_back(vector3f(0.574951, -0.652124, 1.152980));
    vertices.push_back(vector3f(-1.081299, -0.652124, 1.152980));
    vertices.push_back(vector3f(0.512451, -0.636499, 1.152980));
    vertices.push_back(vector3f(-1.018799, -0.636499, 1.152980));
    vertices.push_back(vector3f(0.590576, -0.558374, 1.152980));
    vertices.push_back(vector3f(-1.096924, -0.558374, 1.152980));
    vertices.push_back(vector3f(0.785889, -0.402124, 1.059230));
    vertices.push_back(vector3f(-1.292236, -0.402124, 1.059230));
    vertices.push_back(vector3f(0.934326, -0.386499, 0.988917));
    vertices.push_back(vector3f(-1.440674, -0.386499, 0.988917));
    vertices.push_back(vector3f(1.004639, -0.488062, 0.981105));
    vertices.push_back(vector3f(-1.510986, -0.488062, 0.981105));
    vertices.push_back(vector3f(0.957764, -0.644312, 0.988917));
    vertices.push_back(vector3f(-1.464111, -0.644312, 0.988917));
    vertices.push_back(vector3f(0.793701, -0.730249, 1.051417));
    vertices.push_back(vector3f(-1.300049, -0.730249, 1.051417));
    vertices.push_back(vector3f(0.629639, -0.745874, 1.207667));
    vertices.push_back(vector3f(-1.135986, -0.745874, 1.207667));
    vertices.push_back(vector3f(0.699951, -0.441187, 1.129542));
    vertices.push_back(vector3f(-1.206299, -0.441187, 1.129542));
    vertices.push_back(vector3f(0.637451, -0.620874, 1.145167));
    vertices.push_back(vector3f(-1.143799, -0.620874, 1.145167));
    vertices.push_back(vector3f(0.684326, -0.667749, 1.137355));
    vertices.push_back(vector3f(-1.190674, -0.667749, 1.137355));
    vertices.push_back(vector3f(0.746826, -0.605249, 1.106105));
    vertices.push_back(vector3f(-1.253174, -0.605249, 1.106105));
    vertices.push_back(vector3f(0.707764, -0.558374, 1.121730));
    vertices.push_back(vector3f(-1.214111, -0.558374, 1.121730));
    vertices.push_back(vector3f(0.762451, -0.495874, 1.098292));
    vertices.push_back(vector3f(-1.268799, -0.495874, 1.098292));
    vertices.push_back(vector3f(0.801514, -0.542749, 1.090480));
    vertices.push_back(vector3f(-1.307861, -0.542749, 1.090480));
    vertices.push_back(vector3f(0.856201, -0.519312, 1.082667));
    vertices.push_back(vector3f(-1.362549, -0.519312, 1.082667));
    vertices.push_back(vector3f(0.832764, -0.456812, 1.082667));
    vertices.push_back(vector3f(-1.339111, -0.456812, 1.082667));
    vertices.push_back(vector3f(0.770264, -0.292749, 0.988917));
    vertices.push_back(vector3f(-1.276611, -0.292749, 0.988917));
    vertices.push_back(vector3f(0.996826, -0.261499, 0.926417));
    vertices.push_back(vector3f(-1.503174, -0.261499, 0.926417));
    vertices.push_back(vector3f(1.114014, -0.433374, 0.973292));
    vertices.push_back(vector3f(-1.620361, -0.433374, 0.973292));
    vertices.push_back(vector3f(1.059326, -0.675562, 0.942042));
    vertices.push_back(vector3f(-1.565674, -0.675562, 0.942042));
    vertices.push_back(vector3f(0.785889, -0.816187, 0.981105));
    vertices.push_back(vector3f(-1.292236, -0.816187, 0.981105));
    vertices.push_back(vector3f(0.535889, -0.855249, 1.145167));
    vertices.push_back(vector3f(-1.042236, -0.855249, 1.145167));
    vertices.push_back(vector3f(0.606201, -0.347437, 1.090480));
    vertices.push_back(vector3f(-1.112549, -0.347437, 1.090480));

    std::vector<int> indices = {
        47,  3,   45,  4,   48,  46,  45,  5,   43,  6,   46,  44,  3,   7,
        5,   8,   4,   6,   1,   9,   3,   10,  2,   4,   11,  15,  9,   16,
        12,  10,  9,   17,  7,   18,  10,  8,   21,  17,  15,  22,  18,  20,
        13,  21,  15,  22,  14,  16,  23,  27,  21,  28,  24,  22,  27,  19,
        21,  28,  20,  30,  33,  29,  27,  34,  30,  32,  35,  27,  25,  36,
        28,  34,  37,  33,  35,  38,  34,  40,  39,  31,  33,  40,  32,  42,
        45,  41,  39,  46,  42,  44,  47,  39,  37,  48,  40,  46,  37,  49,
        47,  38,  50,  52,  35,  51,  37,  36,  52,  54,  25,  53,  35,  26,
        54,  56,  23,  55,  25,  24,  56,  58,  23,  59,  57,  60,  24,  58,
        13,  63,  59,  64,  14,  60,  11,  65,  63,  66,  12,  64,  1,   49,
        65,  50,  2,   66,  61,  65,  49,  50,  66,  62,  63,  65,  61,  62,
        66,  64,  61,  59,  63,  64,  60,  62,  61,  57,  59,  60,  58,  62,
        61,  55,  57,  58,  56,  62,  61,  53,  55,  56,  54,  62,  61,  51,
        53,  54,  52,  62,  61,  49,  51,  52,  50,  62,  174, 91,  89,  175,
        91,  176, 172, 89,  87,  173, 90,  175, 85,  172, 87,  173, 86,  88,
        83,  170, 85,  171, 84,  86,  81,  168, 83,  169, 82,  84,  79,  146,
        164, 147, 80,  165, 94,  146, 92,  95,  147, 149, 94,  150, 148, 151,
        95,  149, 98,  150, 96,  99,  151, 153, 100, 152, 98,  101, 153, 155,
        102, 154, 100, 103, 155, 157, 102, 158, 156, 159, 103, 157, 106, 158,
        104, 107, 159, 161, 108, 160, 106, 109, 161, 163, 67,  162, 108, 67,
        163, 68,  128, 162, 110, 129, 163, 161, 128, 158, 160, 159, 129, 161,
        156, 179, 126, 157, 180, 159, 154, 126, 124, 155, 127, 157, 152, 124,
        122, 153, 125, 155, 150, 122, 120, 151, 123, 153, 148, 120, 118, 149,
        121, 151, 146, 118, 116, 147, 119, 149, 164, 116, 114, 165, 117, 147,
        114, 177, 164, 177, 115, 165, 162, 112, 110, 163, 113, 68,  112, 178,
        183, 178, 113, 184, 181, 178, 177, 182, 178, 184, 135, 176, 174, 176,
        136, 175, 133, 174, 172, 175, 134, 173, 133, 170, 131, 134, 171, 173,
        166, 185, 168, 186, 167, 169, 131, 168, 185, 169, 132, 186, 190, 187,
        144, 190, 188, 189, 187, 69,  185, 188, 69,  189, 131, 69,  130, 132,
        69,  186, 142, 191, 144, 192, 143, 145, 140, 193, 142, 194, 141, 143,
        197, 140, 139, 198, 141, 196, 71,  139, 138, 71,  139, 198, 144, 70,
        190, 145, 70,  192, 191, 208, 70,  192, 208, 207, 71,  200, 197, 201,
        71,  198, 197, 202, 195, 203, 198, 196, 202, 193, 195, 203, 194, 205,
        193, 206, 191, 207, 194, 192, 204, 200, 199, 205, 201, 203, 199, 206,
        204, 207, 199, 205, 139, 164, 177, 165, 139, 177, 140, 211, 164, 212,
        141, 165, 144, 211, 142, 145, 212, 214, 187, 213, 144, 188, 214, 167,
        209, 166, 81,  210, 167, 214, 215, 213, 209, 216, 214, 212, 79,  211,
        215, 212, 80,  216, 130, 222, 131, 130, 223, 72,  133, 222, 220, 223,
        134, 221, 135, 220, 218, 221, 136, 219, 137, 218, 217, 219, 137, 217,
        218, 231, 217, 219, 231, 230, 218, 227, 229, 228, 219, 230, 220, 225,
        227, 226, 221, 228, 72,  225, 222, 72,  226, 224, 224, 229, 225, 230,
        224, 226, 225, 229, 227, 228, 230, 226, 183, 234, 232, 235, 184, 233,
        112, 232, 254, 233, 113, 255, 112, 256, 110, 113, 257, 255, 114, 234,
        181, 115, 235, 253, 114, 250, 252, 251, 115, 253, 116, 248, 250, 249,
        117, 251, 118, 246, 248, 247, 119, 249, 120, 244, 246, 245, 121, 247,
        124, 244, 122, 125, 245, 243, 126, 242, 124, 127, 243, 241, 126, 236,
        240, 237, 127, 241, 179, 238, 236, 239, 180, 237, 128, 256, 238, 257,
        129, 239, 256, 276, 238, 257, 277, 259, 236, 276, 278, 277, 237, 279,
        236, 274, 240, 237, 275, 279, 240, 272, 242, 241, 273, 275, 244, 272,
        270, 273, 245, 271, 244, 268, 246, 245, 269, 271, 248, 268, 266, 269,
        249, 267, 248, 264, 250, 249, 265, 267, 250, 262, 252, 251, 263, 265,
        234, 262, 280, 263, 235, 281, 256, 260, 258, 261, 257, 259, 254, 282,
        260, 283, 255, 261, 232, 280, 282, 281, 233, 283, 67,  284, 73,  285,
        67,  73,  108, 286, 284, 287, 109, 285, 104, 286, 106, 105, 287, 289,
        102, 288, 104, 103, 289, 291, 100, 290, 102, 101, 291, 293, 100, 294,
        292, 295, 101, 293, 96,  294, 98,  97,  295, 297, 96,  298, 296, 299,
        97,  297, 94,  300, 298, 301, 95,  299, 309, 338, 308, 309, 339, 329,
        308, 336, 307, 308, 337, 339, 307, 340, 306, 307, 341, 337, 89,  306,
        340, 306, 90,  341, 87,  340, 334, 341, 88,  335, 85,  334, 330, 335,
        86,  331, 83,  330, 332, 331, 84,  333, 330, 338, 332, 339, 331, 333,
        334, 336, 330, 335, 337, 341, 332, 328, 326, 333, 329, 339, 81,  332,
        326, 333, 82,  327, 342, 215, 209, 343, 216, 345, 326, 209, 81,  327,
        210, 343, 215, 346, 79,  216, 347, 345, 346, 92,  79,  347, 93,  301,
        324, 304, 77,  325, 304, 353, 352, 78,  304, 353, 78,  351, 78,  348,
        305, 349, 78,  305, 305, 328, 309, 329, 305, 309, 328, 342, 326, 329,
        343, 349, 296, 318, 310, 319, 297, 311, 316, 77,  76,  317, 77,  325,
        358, 303, 302, 359, 303, 357, 303, 354, 75,  355, 303, 75,  75,  316,
        76,  317, 75,  76,  292, 362, 364, 363, 293, 365, 364, 368, 366, 369,
        365, 367, 366, 370, 372, 371, 367, 373, 372, 376, 374, 377, 373, 375,
        378, 376, 314, 379, 377, 375, 316, 374, 378, 375, 317, 379, 354, 372,
        374, 373, 355, 375, 356, 366, 372, 367, 357, 373, 358, 364, 366, 365,
        359, 367, 292, 360, 290, 293, 361, 365, 360, 302, 74,  361, 302, 359,
        284, 288, 290, 289, 285, 291, 284, 360, 74,  361, 285, 74,  73,  284,
        74,  74,  285, 73,  296, 362, 294, 297, 363, 311, 310, 368, 362, 369,
        311, 363, 312, 370, 368, 371, 313, 369, 376, 382, 314, 377, 383, 371,
        350, 384, 348, 351, 385, 387, 384, 320, 318, 385, 321, 387, 298, 384,
        318, 385, 299, 319, 300, 342, 384, 343, 301, 385, 342, 348, 384, 385,
        349, 343, 300, 346, 344, 345, 347, 301, 322, 378, 314, 323, 379, 381,
        378, 324, 316, 379, 325, 381, 386, 322, 320, 387, 323, 381, 352, 386,
        350, 353, 387, 381, 324, 380, 352, 353, 381, 325, 388, 402, 400, 389,
        403, 415, 400, 404, 398, 405, 401, 399, 404, 396, 398, 405, 397, 407,
        406, 394, 396, 407, 395, 409, 408, 392, 394, 409, 393, 411, 392, 412,
        390, 413, 393, 391, 410, 418, 412, 419, 411, 413, 408, 420, 410, 421,
        409, 411, 424, 408, 406, 425, 409, 423, 426, 406, 404, 427, 407, 425,
        428, 404, 402, 429, 405, 427, 402, 416, 428, 417, 403, 429, 320, 442,
        318, 321, 443, 445, 390, 444, 320, 391, 445, 413, 310, 442, 312, 443,
        311, 313, 382, 414, 388, 415, 383, 389, 412, 440, 444, 441, 413, 445,
        446, 440, 438, 447, 441, 445, 434, 438, 436, 439, 435, 437, 448, 434,
        432, 449, 435, 447, 448, 450, 430, 449, 451, 433, 430, 416, 414, 431,
        417, 451, 312, 430, 382, 431, 313, 383, 442, 448, 312, 443, 449, 447,
        442, 444, 446, 447, 445, 443, 416, 452, 476, 453, 417, 477, 432, 452,
        450, 433, 453, 463, 432, 460, 462, 461, 433, 463, 436, 460, 434, 437,
        461, 459, 438, 458, 436, 439, 459, 457, 438, 454, 456, 455, 439, 457,
        440, 474, 454, 475, 441, 455, 428, 476, 464, 477, 429, 465, 426, 464,
        466, 465, 427, 467, 424, 466, 468, 467, 425, 469, 424, 470, 422, 425,
        471, 469, 422, 472, 420, 423, 473, 471, 420, 474, 418, 421, 475, 473,
        456, 478, 458, 457, 479, 481, 480, 484, 478, 481, 485, 483, 484, 488,
        486, 489, 485, 487, 488, 492, 486, 489, 493, 491, 464, 486, 492, 487,
        465, 493, 484, 476, 452, 485, 477, 487, 462, 484, 452, 463, 485, 479,
        458, 462, 460, 463, 459, 461, 474, 456, 454, 475, 457, 481, 472, 480,
        474, 481, 473, 475, 488, 472, 470, 489, 473, 483, 490, 470, 468, 491,
        471, 489, 466, 490, 468, 491, 467, 469, 464, 492, 466, 467, 493, 465,
        392, 504, 502, 505, 393, 503, 394, 502, 500, 503, 395, 501, 394, 498,
        396, 395, 499, 501, 396, 496, 398, 397, 497, 499, 398, 494, 400, 399,
        495, 497, 400, 506, 388, 401, 507, 495, 502, 506, 494, 503, 507, 505,
        494, 500, 502, 501, 495, 503, 496, 498, 500, 501, 499, 497, 382, 506,
        314, 383, 507, 389, 314, 504, 322, 505, 315, 323, 320, 504, 390, 505,
        321, 391, 47,  1,   3,   4,   2,   48,  45,  3,   5,   6,   4,   46,
        3,   9,   7,   8,   10,  4,   1,   11,  9,   10,  12,  2,   11,  13,
        15,  16,  14,  12,  9,   15,  17,  18,  16,  10,  21,  19,  17,  22,
        16,  18,  13,  23,  21,  22,  24,  14,  23,  25,  27,  28,  26,  24,
        27,  29,  19,  28,  22,  20,  33,  31,  29,  34,  28,  30,  35,  33,
        27,  36,  26,  28,  37,  39,  33,  38,  36,  34,  39,  41,  31,  40,
        34,  32,  45,  43,  41,  46,  40,  42,  47,  45,  39,  48,  38,  40,
        37,  51,  49,  38,  48,  50,  35,  53,  51,  36,  38,  52,  25,  55,
        53,  26,  36,  54,  23,  57,  55,  24,  26,  56,  23,  13,  59,  60,
        14,  24,  13,  11,  63,  64,  12,  14,  11,  1,   65,  66,  2,   12,
        1,   47,  49,  50,  48,  2,   174, 176, 91,  175, 90,  91,  172, 174,
        89,  173, 88,  90,  85,  170, 172, 173, 171, 86,  83,  168, 170, 171,
        169, 84,  81,  166, 168, 169, 167, 82,  79,  92,  146, 147, 93,  80,
        94,  148, 146, 95,  93,  147, 94,  96,  150, 151, 97,  95,  98,  152,
        150, 99,  97,  151, 100, 154, 152, 101, 99,  153, 102, 156, 154, 103,
        101, 155, 102, 104, 158, 159, 105, 103, 106, 160, 158, 107, 105, 159,
        108, 162, 160, 109, 107, 161, 67,  68,  162, 67,  109, 163, 128, 160,
        162, 129, 111, 163, 128, 179, 158, 159, 180, 129, 156, 158, 179, 157,
        127, 180, 154, 156, 126, 155, 125, 127, 152, 154, 124, 153, 123, 125,
        150, 152, 122, 151, 121, 123, 148, 150, 120, 149, 119, 121, 146, 148,
        118, 147, 117, 119, 164, 146, 116, 165, 115, 117, 114, 181, 177, 177,
        182, 115, 162, 68,  112, 163, 111, 113, 112, 68,  178, 178, 68,  113,
        181, 183, 178, 182, 177, 178, 135, 137, 176, 176, 137, 136, 133, 135,
        174, 175, 136, 134, 133, 172, 170, 134, 132, 171, 166, 187, 185, 186,
        188, 167, 131, 170, 168, 169, 171, 132, 190, 189, 187, 190, 145, 188,
        187, 189, 69,  188, 186, 69,  131, 185, 69,  132, 130, 69,  142, 193,
        191, 192, 194, 143, 140, 195, 193, 194, 196, 141, 197, 195, 140, 198,
        139, 141, 71,  197, 139, 71,  138, 139, 144, 191, 70,  145, 190, 70,
        191, 206, 208, 192, 70,  208, 71,  199, 200, 201, 199, 71,  197, 200,
        202, 203, 201, 198, 202, 204, 193, 203, 196, 194, 193, 204, 206, 207,
        205, 194, 204, 202, 200, 205, 199, 201, 199, 208, 206, 207, 208, 199,
        139, 140, 164, 165, 141, 139, 140, 142, 211, 212, 143, 141, 144, 213,
        211, 145, 143, 212, 187, 166, 213, 188, 145, 214, 209, 213, 166, 210,
        82,  167, 215, 211, 213, 216, 210, 214, 79,  164, 211, 212, 165, 80,
        130, 72,  222, 130, 132, 223, 133, 131, 222, 223, 132, 134, 135, 133,
        220, 221, 134, 136, 137, 135, 218, 219, 136, 137, 218, 229, 231, 219,
        217, 231, 218, 220, 227, 228, 221, 219, 220, 222, 225, 226, 223, 221,
        72,  224, 225, 72,  223, 226, 224, 231, 229, 230, 231, 224, 183, 181,
        234, 235, 182, 184, 112, 183, 232, 233, 184, 113, 112, 254, 256, 113,
        111, 257, 114, 252, 234, 115, 182, 235, 114, 116, 250, 251, 117, 115,
        116, 118, 248, 249, 119, 117, 118, 120, 246, 247, 121, 119, 120, 122,
        244, 245, 123, 121, 124, 242, 244, 125, 123, 245, 126, 240, 242, 127,
        125, 243, 126, 179, 236, 237, 180, 127, 179, 128, 238, 239, 129, 180,
        128, 110, 256, 257, 111, 129, 256, 258, 276, 257, 239, 277, 236, 238,
        276, 277, 239, 237, 236, 278, 274, 237, 241, 275, 240, 274, 272, 241,
        243, 273, 244, 242, 272, 273, 243, 245, 244, 270, 268, 245, 247, 269,
        248, 246, 268, 269, 247, 249, 248, 266, 264, 249, 251, 265, 250, 264,
        262, 251, 253, 263, 234, 252, 262, 263, 253, 235, 256, 254, 260, 261,
        255, 257, 254, 232, 282, 283, 233, 255, 232, 234, 280, 281, 235, 233,
        67,  108, 284, 285, 109, 67,  108, 106, 286, 287, 107, 109, 104, 288,
        286, 105, 107, 287, 102, 290, 288, 103, 105, 289, 100, 292, 290, 101,
        103, 291, 100, 98,  294, 295, 99,  101, 96,  296, 294, 97,  99,  295,
        96,  94,  298, 299, 95,  97,  94,  92,  300, 301, 93,  95,  309, 328,
        338, 309, 308, 339, 308, 338, 336, 308, 307, 337, 307, 336, 340, 307,
        306, 341, 89,  91,  306, 306, 91,  90,  87,  89,  340, 341, 90,  88,
        85,  87,  334, 335, 88,  86,  83,  85,  330, 331, 86,  84,  330, 336,
        338, 339, 337, 331, 334, 340, 336, 335, 331, 337, 332, 338, 328, 333,
        327, 329, 81,  83,  332, 333, 84,  82,  342, 344, 215, 343, 210, 216,
        326, 342, 209, 327, 82,  210, 215, 344, 346, 216, 80,  347, 346, 300,
        92,  347, 80,  93,  324, 352, 304, 325, 77,  304, 352, 350, 78,  353,
        304, 78,  78,  350, 348, 349, 351, 78,  305, 348, 328, 329, 349, 305,
        328, 348, 342, 329, 327, 343, 296, 298, 318, 319, 299, 297, 316, 324,
        77,  317, 76,  77,  358, 356, 303, 359, 302, 303, 303, 356, 354, 355,
        357, 303, 75,  354, 316, 317, 355, 75,  292, 294, 362, 363, 295, 293,
        364, 362, 368, 369, 363, 365, 366, 368, 370, 371, 369, 367, 372, 370,
        376, 377, 371, 373, 378, 374, 376, 379, 315, 377, 316, 354, 374, 375,
        355, 317, 354, 356, 372, 373, 357, 355, 356, 358, 366, 367, 359, 357,
        358, 360, 364, 365, 361, 359, 292, 364, 360, 293, 291, 361, 360, 358,
        302, 361, 74,  302, 284, 286, 288, 289, 287, 285, 284, 290, 360, 361,
        291, 285, 296, 310, 362, 297, 295, 363, 310, 312, 368, 369, 313, 311,
        312, 382, 370, 371, 383, 313, 376, 370, 382, 377, 315, 383, 350, 386,
        384, 351, 349, 385, 384, 386, 320, 385, 319, 321, 298, 300, 384, 385,
        301, 299, 300, 344, 342, 343, 345, 301, 322, 380, 378, 323, 315, 379,
        378, 380, 324, 379, 317, 325, 386, 380, 322, 387, 321, 323, 352, 380,
        386, 353, 351, 387, 388, 414, 402, 389, 401, 403, 400, 402, 404, 405,
        403, 401, 404, 406, 396, 405, 399, 397, 406, 408, 394, 407, 397, 395,
        408, 410, 392, 409, 395, 393, 392, 410, 412, 413, 411, 393, 410, 420,
        418, 419, 421, 411, 408, 422, 420, 421, 423, 409, 424, 422, 408, 425,
        407, 409, 426, 424, 406, 427, 405, 407, 428, 426, 404, 429, 403, 405,
        402, 414, 416, 417, 415, 403, 320, 444, 442, 321, 319, 443, 390, 412,
        444, 391, 321, 445, 310, 318, 442, 443, 319, 311, 382, 430, 414, 415,
        431, 383, 412, 418, 440, 441, 419, 413, 446, 444, 440, 447, 439, 441,
        434, 446, 438, 439, 447, 435, 448, 446, 434, 449, 433, 435, 448, 432,
        450, 449, 431, 451, 430, 450, 416, 431, 415, 417, 312, 448, 430, 431,
        449, 313, 442, 446, 448, 443, 313, 449, 416, 450, 452, 453, 451, 417,
        432, 462, 452, 433, 451, 453, 432, 434, 460, 461, 435, 433, 436, 458,
        460, 437, 435, 461, 438, 456, 458, 439, 437, 459, 438, 440, 454, 455,
        441, 439, 440, 418, 474, 475, 419, 441, 428, 416, 476, 477, 417, 429,
        426, 428, 464, 465, 429, 427, 424, 426, 466, 467, 427, 425, 424, 468,
        470, 425, 423, 471, 422, 470, 472, 423, 421, 473, 420, 472, 474, 421,
        419, 475, 456, 480, 478, 457, 459, 479, 480, 482, 484, 481, 479, 485,
        484, 482, 488, 489, 483, 485, 488, 490, 492, 489, 487, 493, 464, 476,
        486, 487, 477, 465, 484, 486, 476, 485, 453, 477, 462, 478, 484, 463,
        453, 485, 458, 478, 462, 463, 479, 459, 474, 480, 456, 475, 455, 457,
        472, 482, 480, 481, 483, 473, 488, 482, 472, 489, 471, 473, 490, 488,
        470, 491, 469, 471, 466, 492, 490, 491, 493, 467, 392, 390, 504, 505,
        391, 393, 394, 392, 502, 503, 393, 395, 394, 500, 498, 395, 397, 499,
        396, 498, 496, 397, 399, 497, 398, 496, 494, 399, 401, 495, 400, 494,
        506, 401, 389, 507, 502, 504, 506, 503, 495, 507, 494, 496, 500, 501,
        497, 495, 382, 388, 506, 383, 315, 507, 314, 506, 504, 505, 507, 315,
        320, 322, 504, 505, 323, 321};

    for (auto& i : indices)
    {
        i -= 1;
    }

    m_triangles.clear();
    int numTriangles = (int)indices.size() / 3;
    for (int i = 0; i < numTriangles; i++)
    {
        int ia = indices[i * 3 + 0];
        int ib = indices[i * 3 + 1];
        int ic = indices[i * 3 + 2];

        m_triangles.push_back(
            Triangle(vertices[ia], vertices[ib], vertices[ic]));
    }

    for (const auto& t : m_triangles)
    {
        m_boundingBox.expandBy(t.boundingBox());
    }
}

void DebugMesh::initRandomMesh()
{
    int num_tris = 1000;
    m_triangles.clear();
    for (int i = 0; i < num_tris; i++)
    {
        auto a = vector3f(randf() * 2.0f - 1.0f);
        auto b = vector3f(randf() * 2.0f - 1.0f);
        auto c = vector3f(randf() * 2.0f - 1.0f);

        m_triangles.push_back(Triangle(a, b, c));
    }

    for (const auto& t : m_triangles)
    {
        m_boundingBox.expandBy(t.boundingBox());
    }
}

bool DebugMesh::intersect(Ray& ray, HitInfo& hit_info) const
{
    hit_info.reset();

    // coarse grained culling using mesh bounding box
    if (false == m_boundingBox.overlap(ray))
    {
        return false;
    }

    HitInfo hit_info_child;  // reused for all contained triangles

    for (const auto& t : m_triangles)
    {
        bool hit = t.intersect(ray, hit_info_child);
        if (hit && hit_info_child < hit_info)
        {
            hit_info = hit_info_child;
        }
    }

    if (hit_info.m_t > ray.tmin && hit_info.m_t <= ray.tmax)
    {
        return true;
    }
    else
    {
        hit_info.reset();
        return false;
    }
}
