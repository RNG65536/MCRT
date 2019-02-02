#include <cmath>
#include "bssrdf.h"
#include "constants.h"
#include "numeric.h"

static vec3 reflect(const vec3& wi, const vec3& nl)
{
    return wi - (2 * dot(wi, nl)) * nl;
}

static vec3 refract(const vec3& wi, const vec3& nl, float eta)
{
    float dotValue(dot(nl, wi));
    float k(1.0 - eta * eta * (1.0 - dotValue * dotValue));
    return (eta * wi - (eta * dotValue + std::sqrt(k)) * nl) * float(k >= 0);
}

static float fresnel_R(float cos_theta, float eta)
{
    float sin_theta_t_sqr = 1.0 / (eta * eta) * (1.0 - cos_theta * cos_theta);
    if (sin_theta_t_sqr >= 1.0)
    {
        return 1.0;
    }
    float cos_theta_t = std::sqrt(1.0 - sin_theta_t_sqr);
    float a = eta * cos_theta;
    float b = eta * cos_theta_t;
    float r_s = (cos_theta - b) / (cos_theta + b);
    float r_p = (a - cos_theta_t) / (a + cos_theta_t);
    return (r_s * r_s + r_p * r_p) * 0.5;
}

static float two_C1(float n)
{
    float r;
    if (n > 1.0)
    {
        r = -9.23372 + n * (22.2272 + n * (-20.9292 + n * (10.2291 + n * (-2.54396 + n * 0.254913))));
    }
    else
    {
        r = 0.919317 + n * (-3.4793 + n * (6.75335 + n * (-7.80989 + n * (4.98554 - n * 1.36881))));
    }
    return r;
}

static float three_C2(float n)
{
    float r;
    if (n > 1.0)
    {
        r = -1641.1 + n * (1213.67 + n * (-568.556 + n * (164.798 + n * (-27.0181 + n * 1.91826))));
        r += (((135.926 / n) - 656.175) / n + 1376.53) / n;
    }
    else
    {
        r = 0.828421 + n * (-2.62051 + n * (3.36231 + n * (-1.95284 + n * (0.236494 + n * 0.145787))));
    }
    return r;
}

float S_d_prime(const vec3& x, const vec3& w, float r, const vec3& n, float sigma_tr, float D, float Cp_norm, float Cp, float Ce)
{
    float s_tr_r = sigma_tr * r;
    float s_tr_r_one = 1.0 + s_tr_r;
    float x_dot_w = dot(x, w);
    float r_sqr = r * r;
    float t0 = Cp_norm * M_1_4PIPI * std::exp(-s_tr_r) / (r * r_sqr);
    float t1 = r_sqr / D + 3.0 * s_tr_r_one * x_dot_w;
    float t2 = 3.0 * D * s_tr_r_one * dot(w, n);
    float t3 = (s_tr_r_one + 3.0 * D * (3.0 * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * dot(x, n);
    return t0 * (Cp * t1 - Ce * (t2 - t3));
}
float dir_bssrdf(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no,
    float sigma_s, float sigma_a, float g, float eta, float Cp_norm, float Cp, float Ce, float A)
{
    float sigma_t = sigma_s + sigma_a;
    float sigma_t_p = sigma_s * (1.0 - g) + sigma_a;
    float D = 1.0 / (3.0 * sigma_t_p);
    float sigma_tr = std::sqrt(sigma_a / D);
    float de = 2.131 * D / std::sqrt(1.0 - sigma_a / sigma_t_p);
    vec3 x = xo - xi;
    float r_sqr = dot(x, x);
    vec3 wr = refract(-wi, ni, 1.0 / eta);
    vec3 ni_ast = r_sqr < 1.0e-12 ? ni : normalize(ni * r_sqr - x * dot(x, ni));
    vec3 wv = reflect(wr, ni_ast);
    float x_dot_wr = dot(x, wr);
    float numerator = r_sqr - x_dot_wr * x_dot_wr;
    float cos_beta = -std::sqrt(numerator / (r_sqr + de * de));
    float mu0 = -dot(no, wr);
    float dr = mu0 > 0.0 ? std::sqrt(D * mu0 * (D * mu0 - de * cos_beta * 2.0) + r_sqr)
        : std::sqrt(1.0 / (3.0 * sigma_t * 3.0 * sigma_t) + r_sqr);
    vec3 xoxv = xo - (xi + ni_ast * (2.0 * A * de));
    float dv = length(xoxv);
    return S_d_prime(x, wr, dr, no, sigma_tr, D, Cp_norm, Cp, Ce)
        - S_d_prime(xoxv, wv, dv, no, sigma_tr, D, Cp_norm, Cp, Ce);
}

//////////////////////////////////////////////////////////////////////////

float std_bssrdf(float r, float sigma_s, float sigma_a, float g, float A)
{
    float sigma_t_p = sigma_s * (1.0 - g) + sigma_a;
    float D = 1.0 / (3.0 * sigma_t_p);
    float sigma_tr = std::sqrt(sigma_a / D);
    float zr = 3.0 * D;
    float zv = zr + 4.0 * A * D;

    float dr = std::sqrt(zr * zr + r * r);
    float dv = std::sqrt(zv * zv + r * r);
    float real = zr * (1.0 + sigma_tr * dr)
        / (dr * dr * dr) * std::exp(-sigma_tr * dr);
    float virt = zv * (1.0 + sigma_tr * dv)
        / (dv * dv * dv) * std::exp(-sigma_tr * dv);

    float albedo = 1.0 - sigma_a * zr; // alpha_prime
    return albedo * M_1_4PIPI * (real + virt);
}

float phase_HG(float cos_theta, float g)
{
    float g_sqr = g * g;
    float demon = 1.0 + g_sqr - g * (2.0 * cos_theta);
    return (1.0 - g_sqr) / pow(demon, 1.5) * M_1_4PI;
}

float single_approx(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no,
    float sigma_s, float sigma_a, float g, float eta)
{
    float sigma_t = sigma_s + sigma_a;
    float sigma_s_p = sigma_s * (1.0 - g);
    float sigma_t_p = sigma_s_p + sigma_a;
    vec3 w12 = refract(-wi, ni, 1.0 / eta);
    float mu0 = std::abs(dot(no, w12));
    float d1 = mu0 / (3.0 * sigma_t_p); // 1.0/sigma_t_p; //
    vec3 xs = xi + w12 * d1;
    vec3 w21 = xo - xs;
    float d2 = length(w21);
    w21 /= d2;
    return sigma_s_p * d1 * phase_HG(dot(w12, w21), g) * std::exp(-sigma_t_p * d1 - sigma_t * d2) / (d2 * d2);
}

//////////////////////////////////////////////////////////////////////////

vec3 DirectionalDipoleBSSRDF::evaluate(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no, const vec3& wo)
{
    return vec3(
        evaluate_channel(xi, ni, wi, xo, no, wo, 0),
        evaluate_channel(xi, ni, wi, xo, no, wo, 1),
        evaluate_channel(xi, ni, wi, xo, no, wo, 2));
}

DirectionalDipoleBSSRDF::DirectionalDipoleBSSRDF(const vec3& sigma_s, const vec3& sigma_a, float g, const vec3& eta) :
sigma_s(sigma_s), sigma_a(sigma_a), g(g), eta(eta)
{
    for (int i = 0; i < 3; i++)
    {
        Cp_norm[i] = 1.0 / (1.0 - two_C1(1.0 / eta[i]));
        Cp[i] = (1.0 - two_C1(eta[i])) / 4.0;
        Ce[i] = (1.0 - three_C2(eta[i])) / 2.0;
        A[i] = (1.0 - Ce[i]) / (2.0 * Cp[i]);

        float _fdr = -1.440 / (eta[i] * eta[i]) + 0.710 / eta[i] + 0.668 + 0.0636 * eta[i];
        Fdt[i] = 1 - _fdr;
    }
}

float DirectionalDipoleBSSRDF::evaluate_channel(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no, const vec3& wo, int i)
{
    float sss = dir_bssrdf(xi, ni, wi, xo, no,
        sigma_s[i], sigma_a[i], g, eta[i], Cp_norm[i], Cp[i], Ce[i], A[i]);

    // these two results are amazingly close
    // 1)
    float T12 = 1.0 - fresnel_R(dot(wi, ni), eta[i]);
    sss *= T12;
    sss *= M_4PI * Cp[i]; // why this ??
    // 2)
//     sss *= Fdt[i] * M_PI; // there's no where to put this pi, but it helps to modulate the result nicely when F_dt is used

    sss = f_max(sss, 0);
    return sss;
}

//////////////////////////////////////////////////////////////////////////

vec3 StandardDipoleBSSRDF::evaluate(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no, const vec3& wo)
{
    return vec3(
        evaluate_channel(xi, ni, wi, xo, no, wo, 0),
        evaluate_channel(xi, ni, wi, xo, no, wo, 1),
        evaluate_channel(xi, ni, wi, xo, no, wo, 2));
}

StandardDipoleBSSRDF::StandardDipoleBSSRDF(const vec3& sigma_s, const vec3& sigma_a, float g, const vec3& eta) :
sigma_s(sigma_s), sigma_a(sigma_a), g(g), eta(eta)
{
    for (int i = 0; i < 3; i++)
    {
        Cp[i] = (1.0 - two_C1(eta[i])) / 4.0;
        float _fdr = -1.440 / (eta[i] * eta[i]) + 0.710 / eta[i] + 0.668 + 0.0636 * eta[i];
        A[i] = (1 + _fdr) / (1 - _fdr);
        Fdt[i] = 1 - _fdr;
    }
}

float StandardDipoleBSSRDF::evaluate_channel(const vec3& xi, const vec3& ni, const vec3& wi, const vec3& xo, const vec3& no, const vec3& wo, int i)
{
    float sss = std_bssrdf(length(xo - xi), sigma_s[i], sigma_a[i], g, A[i]);
    sss += single_approx(xi, ni, wi, xo, no, sigma_s[i], sigma_a[i], g, eta[i]);

    // 1)
    float T12 = 1.0 - fresnel_R(dot(wi, ni), eta[i]);
    sss *= T12;
    // 2)
//     sss *= Fdt[i];

    sss *= M_4PI * Cp[i]; // is it ??

    sss = f_max(sss, 0);
    return sss;
}

BSSRDF::~BSSRDF()
{

}
