#include <sstream>
#include "consoledebug.h"
#include "vectors.h"

void debugPrint(const vec3& a)
{
    printf("vec3: %f, %f, %f\n", a.x, a.y, a.z);
}

std::string str(const vec3& v)
{
    std::stringstream ss;
    ss << "( " << v.x << ", " << v.y << ", " << v.z << " )";
    return ss.str();
}

std::string str(int i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

std::string str(float f)
{
    std::stringstream ss;
    ss << f;
    return ss.str();
}
