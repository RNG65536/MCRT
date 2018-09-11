#ifndef vectors_h__
#define vectors_h__

// vector and operators

class vec3
{
public:
    float x, y, z; /*char w[4]*/;  //::dummy member to help alignment

    friend vec3 operator*(float a, const vec3 &b);

    vec3();
    explicit vec3(float a);
    vec3(float _x, float _y, float _z);
    vec3 operator+(const float &v) const;
    vec3 operator+(const vec3 &v) const;
    vec3 operator-(const float &v) const;
    vec3 operator-(const vec3 &v) const;
    vec3& operator+=(const vec3 &v);
    vec3& operator-=(const vec3 &v);
    vec3 operator*(float f) const;
    vec3 operator*(const vec3& v) const;
    vec3 operator/(const vec3& v) const;
    vec3& operator*=(float f);
    vec3& operator*=(const vec3& v);
    vec3 operator-() const;
    vec3 operator/(float f) const;
    vec3 &operator/=(float f);
//     float operator % (const vec3 &b) const;
    float operator[](int i) const;
    float &operator[](int i);
    vec3 norm() const;
    vec3 normalize() const;
    float lengthSquared() const;
    float length() const;
    bool operator==(const vec3 &v) const;
    bool operator!=(const vec3 &v) const;
//     vec3 mult(const vec3 &v2) const;
//     vec3 operator ^ (const vec3 &b) const;
};

void debugPrint(const vec3& a);

vec3 operator + (const float a, const vec3& v);
vec3 operator - (const float a, const vec3& v);
vec3 mult(const vec3 &v1, const vec3 &v2);
vec3 div(const vec3 &v1, const vec3 &v2);
float dot(const vec3 &v1, const vec3 &v2);
vec3 cross(const vec3 &v1, const vec3 &v2);
float length(const vec3& v);
float lengthSquared(const vec3& v);
vec3 normalize(const vec3& v);

vec3 exp(vec3 a);
vec3 sqrt(vec3 a);
vec3 sq(vec3 a);

vec3 f_min(const vec3& a, const vec3& b);
vec3 f_max(const vec3& a, const vec3& b);

float maxComponent(const vec3& v);
float minComponent(const vec3& v);
bool isZero(const vec3& v);

class int3
{
public:
    int x, y, z;
    int3();
    int3(int a, int b, int c);
};

class vec3d
{
public:
    double x, y, z;

    vec3d();
    explicit vec3d(double a);
    vec3d(double _x, double _y, double _z);
    explicit vec3d(const vec3& v);
    vec3 toFloat() const;
    vec3d operator+(const vec3d &v) const;
    vec3d operator-(const vec3d &v) const;
    vec3d operator*(double f) const;
};

double length(const vec3d &v);
vec3d normalize(const vec3d &v);
double dot(const vec3d &v1, const vec3d &v2);

typedef vec3 vec2; // todo : implement vec2
typedef int3 int2;

#endif // vectors_h__
