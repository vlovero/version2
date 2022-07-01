#ifndef VCL_COMPLEX_TRIG_H
#define VCL_COMPLEX_TRIG_H

#include <vcl/vectorclass.h>
#include <vcl/vectormath_trig.h>
#include <vcl/vectormath_exp.h>
#include <vcl/vectormath_hyp.h>
#include <vcl/complexvec1.h>


template <typename T>
constexpr bool is_complex_vector()
{
    bool b = std::is_same<T, Complex1d>::value || std::is_same<T, Complex1f>::value || std::is_same<T, Complex2d>::value || std::is_same<T, Complex2f>::value || std::is_same<T, Complex4d>::value || std::is_same<T, Complex4f>::value || std::is_same<T, Complex8f>::value;
    return b;
}

// perform i * z = i * (re + i * im)
template <typename CVTYPE>
static inline CVTYPE mul_i(CVTYPE const z, const decltype(CVTYPE().to_vector()[0]) s = 1.0)
{
    // i (a + i b) = -b + a i
    using VTYPE = decltype(z.to_vector());
    const VTYPE alt = CVTYPE(-s, s).to_vector();
    VTYPE vec = z.to_vector();


    // swap real and imag parts
    // vec = b + i a (from a + i b)
    if constexpr (vec.size() == 2) {
        vec = permute2<1, 0>(vec);
    }
    else if constexpr (vec.size() == 4) {
        vec = permute4<1, 0, 3, 2>(vec);
    }
    else if constexpr (vec.size() == 8) {
        vec = permute8<1, 0, 3, 2, 5, 4, 7, 6>(vec);
    }
    else if constexpr (vec.size() == 16) {
        vec = permute16<1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 12, 13, 15, 14>(vec);
    }
    vec *= alt;
    return CVTYPE(vec);
}

// perform -i * z = -i * (re + i * im)
template <typename CVTYPE>
static inline CVTYPE nmul_i(CVTYPE const z)
{
    // -i (a + i b) = b - a i
    using VTYPE = decltype(z.to_vector());
    const VTYPE alt = CVTYPE(1.0, -1.0).to_vector();
    VTYPE vec = z.to_vector();


    // swap real and imag parts
    // vec = b + i a (from a + i b)
    if constexpr (vec.size() == 2) {
        vec = permute2<1, 0>(vec);
    }
    else if constexpr (vec.size() == 4) {
        vec = permute4<1, 0, 3, 2>(vec);
    }
    else if constexpr (vec.size() == 8) {
        vec = permute8<1, 0, 3, 2, 5, 4, 7, 6>(vec);
    }
    else if constexpr (vec.size() == 16) {
        vec = permute16<1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 12, 13, 15, 14>(vec);
    }
    vec *= alt;
    return CVTYPE(vec);
}

template <typename CVTYPE>
static inline CVTYPE sin(CVTYPE const z)
{
    static_assert(is_complex_vector<CVTYPE>());
    auto vec = z.to_vector(); // z as vector of interleaved real and imaginary parts
    using VTYPE = decltype(vec);
    // get real and imaginary parts in separate vectors
    VTYPE im, re, c, s, a, ai;

    if constexpr (vec.size() == 2) {
        re = permute2<0, 0>(vec);
        im = permute2<1, 1>(vec);
    }
    else if constexpr (vec.size() == 4) {
        re = permute4<0, 0, 2, 2>(vec);
        im = permute4<1, 1, 3, 3>(vec);
    }
    else if constexpr (vec.size() == 8) {
        re = permute8<0, 0, 2, 2, 4, 4, 6, 6>(vec);
        im = permute8<1, 1, 3, 3, 5, 5, 7, 7>(vec);
    }
    else if constexpr (vec.size() == 16) {
        re = permute16<0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14>(vec);
        im = permute16<1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15>(vec);
    }

    s = sincos(&c, re);
    a = exp(im);
    ai = (decltype(re[0]))(1) / a;
    re = s * (a + ai);
    im = c * (a - ai);

    if constexpr (vec.size() == 2) {
        vec = blend2<0, 3>(re, im);
    }
    else if constexpr (vec.size() == 4) {
        vec = blend4<0, 5, 2, 7>(re, im);
    }
    else if constexpr (vec.size() == 8) {
        vec = blend8<0, 9, 2, 11, 4, 13, 6, 15>(re, im);
    }
    else if constexpr (vec.size() == 16) {
        vec = blend16<0, 17, 2, 19, 4, 21, 6, 23, 8, 25, 10, 27, 12, 29, 14, 31>(re, im);
    }

    return CVTYPE(vec * 0.5);
}

template <typename CVTYPE>
static inline CVTYPE cos(CVTYPE const z)
{
    static_assert(is_complex_vector<CVTYPE>());
    auto vec = z.to_vector(); // z as vector of interleaved real and imaginary parts
    using VTYPE = decltype(vec);
    // get real and imaginary parts in separate vectors
    VTYPE im, re, c, s, a, ai;

    if constexpr (vec.size() == 2) {
        re = permute2<0, 0>(vec);
        im = permute2<1, 1>(vec);
    }
    else if constexpr (vec.size() == 4) {
        re = permute4<0, 0, 2, 2>(vec);
        im = permute4<1, 1, 3, 3>(vec);
    }
    else if constexpr (vec.size() == 8) {
        re = permute8<0, 0, 2, 2, 4, 4, 6, 6>(vec);
        im = permute8<1, 1, 3, 3, 5, 5, 7, 7>(vec);
    }
    else if constexpr (vec.size() == 16) {
        re = permute16<0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14>(vec);
        im = permute16<1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15>(vec);
    }

    s = sincos(&c, re);
    a = exp(im);
    ai = (decltype(re[0]))(1) / a;
    re = c * (a + ai);
    im = s * (ai - a);

    if constexpr (vec.size() == 2) {
        vec = blend2<0, 3>(re, im);
    }
    else if constexpr (vec.size() == 4) {
        vec = blend4<0, 5, 2, 7>(re, im);
    }
    else if constexpr (vec.size() == 8) {
        vec = blend8<0, 9, 2, 11, 4, 13, 6, 15>(re, im);
    }
    else if constexpr (vec.size() == 16) {
        vec = blend16<0, 17, 2, 19, 4, 21, 6, 23, 8, 25, 10, 27, 12, 29, 14, 31>(re, im);
    }

    return CVTYPE(vec * 0.5);
}

template <typename CVTYPE>
static inline CVTYPE tan(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());

    auto vec = z.to_vector(); // z as vector of interleaved real and imaginary parts
    using VTYPE = decltype(vec);

    CVTYPE r;
    const VTYPE one(1.0);
    const CVTYPE eye(0.0, 1.0);
    const VTYPE alt = CVTYPE(1.0, -1.0).to_vector();

    // swap real and imag parts
    // vec = b + i a (from a + i b)
    if constexpr (vec.size() == 2) {
        vec = permute2<1, 0>(vec);
    }
    else if constexpr (vec.size() == 4) {
        vec = permute4<1, 0, 3, 2>(vec);
    }
    else if constexpr (vec.size() == 8) {
        vec = permute8<1, 0, 3, 2, 5, 4, 7, 6>(vec);
    }
    else if constexpr (vec.size() == 16) {
        vec = permute16<1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 12, 13, 15, 14>(vec);
    }
    // vec = b - i a
    vec *= alt;
    z = cexp(CVTYPE(vec));
    r = 1.0 / z;
    return eye * (z - r) / (z + r);
}


template <typename CVTYPE>
static inline CVTYPE exp(CVTYPE const z)
{
    static_assert(is_complex_vector<CVTYPE>());
    return cexp(z);
}


template <typename CVTYPE>
static inline CVTYPE exp10(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    auto vec = z.to_vector();
    using VTYPE = decltype(vec);

    VTYPE im, re, c, s, a, ai;

    // get real and imaginary parts in separate vectors
    if constexpr (vec.size() == 2) {
        re = permute2<0, 0>(vec);
        im = permute2<1, 1>(vec);
    }
    else if constexpr (vec.size() == 4) {
        re = permute4<0, 0, 2, 2>(vec);
        im = permute4<1, 1, 3, 3>(vec);
    }
    else if constexpr (vec.size() == 8) {
        re = permute8<0, 0, 2, 2, 4, 4, 6, 6>(vec);
        im = permute8<1, 1, 3, 3, 5, 5, 7, 7>(vec);
    }
    else if constexpr (vec.size() == 16) {
        re = permute16<0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14>(vec);
        im = permute16<1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15>(vec);
    }

    im *= VM_LN10;
    // sincos(imag * ln(10))
    s = sincos(&c, im);
    a = exp10(re);

    if constexpr (vec.size() == 2) {
        vec = blend2<0, 3>(c, s);
    }
    else if constexpr (vec.size() == 4) {
        vec = blend4<0, 5, 2, 7>(c, s);
    }
    else if constexpr (vec.size() == 8) {
        vec = blend8<0, 9, 2, 11, 4, 13, 6, 15>(c, s);
    }
    else if constexpr (vec.size() == 16) {
        vec = blend16<0, 17, 2, 19, 4, 21, 6, 23, 8, 25, 10, 27, 12, 29, 14, 31>(c, s);
    }

    return CVTYPE(vec * a);
}

template <typename CVTYPE>
static inline CVTYPE exp2(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    auto vec = z.to_vector();
    using VTYPE = decltype(vec);

    VTYPE im, re, c, s, a, ai;

    // get real and imaginary parts in separate vectors
    if constexpr (vec.size() == 2) {
        re = permute2<0, 0>(vec);
        im = permute2<1, 1>(vec);
    }
    else if constexpr (vec.size() == 4) {
        re = permute4<0, 0, 2, 2>(vec);
        im = permute4<1, 1, 3, 3>(vec);
    }
    else if constexpr (vec.size() == 8) {
        re = permute8<0, 0, 2, 2, 4, 4, 6, 6>(vec);
        im = permute8<1, 1, 3, 3, 5, 5, 7, 7>(vec);
    }
    else if constexpr (vec.size() == 16) {
        re = permute16<0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14>(vec);
        im = permute16<1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15>(vec);
    }

    im *= VM_LN2;
    // sincos(imag * ln(10))
    s = sincos(&c, im);
    a = exp10(re);

    if constexpr (vec.size() == 2) {
        vec = blend2<0, 3>(c, s);
    }
    else if constexpr (vec.size() == 4) {
        vec = blend4<0, 5, 2, 7>(c, s);
    }
    else if constexpr (vec.size() == 8) {
        vec = blend8<0, 9, 2, 11, 4, 13, 6, 15>(c, s);
    }
    else if constexpr (vec.size() == 16) {
        vec = blend16<0, 17, 2, 19, 4, 21, 6, 23, 8, 25, 10, 27, 12, 29, 14, 31>(c, s);
    }

    return CVTYPE(vec * a);
}

template <typename CVTYPE>
static inline CVTYPE sqrt(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    return csqrt(z);
}

template <typename CVTYPE>
static inline CVTYPE sinh(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    using VTYPE = decltype(z.to_vector());
    CVTYPE a = cexp(z);
    CVTYPE ai = -0.5 / a;

    VTYPE b = a.to_vector(), bi = ai.to_vector(), half(0.5);
    return CVTYPE(mul_add(half, b, bi));
}

template <typename CVTYPE>
static inline CVTYPE cosh(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    using VTYPE = decltype(z.to_vector());
    CVTYPE a = cexp(z);
    CVTYPE ai = 0.5 / a;

    VTYPE b = a.to_vector(), bi = ai.to_vector(), half(0.5);
    return CVTYPE(mul_add(half, b, bi));
}

template <typename CVTYPE>
static inline CVTYPE tanh(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    CVTYPE a = cexp(z);
    CVTYPE ai = 1.0 / a;

    return (a - ai) / (a + ai);
}

template <typename CVTYPE>
static inline CVTYPE asin(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    CVTYPE eye(0.0, 1.0), one(1.0);
    return -eye * clog(eye * z + csqrt(one - z * z));
}

template <typename CVTYPE>
static inline CVTYPE acos(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    CVTYPE eye(0.0, 1.0), one(1.0);
    return VM_PI_2 + eye * clog(eye * z + csqrt(one - z * z));
}

template <typename CVTYPE>
static inline CVTYPE atan(CVTYPE z)
{
    // (i / 2) log((1 - i (a + i b)) / (1 + i (a + i b)))
    static_assert(is_complex_vector<CVTYPE>());
    using VTYPE = decltype(z.to_vector());

    VTYPE vec = z.to_vector();
    VTYPE alt1 = CVTYPE(-1.0, 1.0).to_vector();
    VTYPE alt2 = CVTYPE(1.0, -1.0).to_vector();

    CVTYPE eye(0.0, 1.0), one(1.0), num, den;

    // swap real and imag parts
    // vec = b + i a (from a + i b)
    if constexpr (vec.size() == 2) {
        vec = permute2<1, 0>(vec);
    }
    else if constexpr (vec.size() == 4) {
        vec = permute4<1, 0, 3, 2>(vec);
    }
    else if constexpr (vec.size() == 8) {
        vec = permute8<1, 0, 3, 2, 5, 4, 7, 6>(vec);
    }
    else if constexpr (vec.size() == 16) {
        vec = permute16<1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 12, 13, 15, 14>(vec);
    }

    // numerator -i (a + i b) = b - a i
    // denominator i (a + i b) = -b + a i
    // return (i / 2) log(num / den)
    num = CVTYPE(mul_add(alt2, vec, one.to_vector()));
    den = CVTYPE(mul_add(alt1, vec, one.to_vector()));
    return mul_i(clog(num / den), 0.5);
}

template <typename CVTYPE>
static inline CVTYPE log(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    return clog(z);
}

template <typename CVTYPE>
static inline CVTYPE asinh(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    const CVTYPE one(1);

    return clog(z + sqrt(one + z * z));
}

template <typename CVTYPE>
static inline CVTYPE acosh(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    const CVTYPE one(1);

    return clog(z + sqrt(z * z - one));
}

template <typename CVTYPE>
static inline CVTYPE atanh(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    CVTYPE one(1), num, den;

    num = one + z;
    den = one - z;
    return 0.5 * clog(num / den);
}

template <typename CVTYPE>
static inline CVTYPE log2(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    return VM_LOG2E * clog(z);
}

template <typename CVTYPE>
static inline CVTYPE log10(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    return VM_LOG10E * clog(z);
}

template <typename CVTYPE>
static inline CVTYPE log1p(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    return clog(z + 1.0);
}

template<typename VTYPE>
static inline VTYPE cosb_expm1_d(VTYPE const initial_x, VTYPE const cosb, VTYPE &exp_x)
{
    // Taylor coefficients, 1/n!
    // Not using minimax approximation because we prioritize precision close to x = 0
    const double p2 = 1. / 2.;
    const double p3 = 1. / 6.;
    const double p4 = 1. / 24.;
    const double p5 = 1. / 120.;
    const double p6 = 1. / 720.;
    const double p7 = 1. / 5040.;
    const double p8 = 1. / 40320.;
    const double p9 = 1. / 362880.;
    const double p10 = 1. / 3628800.;
    const double p11 = 1. / 39916800.;
    const double p12 = 1. / 479001600.;
    const double p13 = 1. / 6227020800.;

    // maximum abs(x), value depends on BA, defined below
    // The lower limit of x is slightly more restrictive than the upper limit.
    // We are specifying the lower limit, except for BA = 1 because it is not used for negative x
    double max_x;

    // data vectors
    VTYPE  x, r, z, n2;
    max_x = 708.39;
    const double ln2d_hi = 0.693145751953125;
    const double ln2d_lo = 1.42860682030941723212E-6;
    x = initial_x;
    r = round(initial_x * VM_LOG2E);
    // subtraction in two steps for higher precision
    x = nmul_add(r, ln2d_hi, x);             //  x -= r * ln2d_hi;
    x = nmul_add(r, ln2d_lo, x);             //  x -= r * ln2d_lo;

    z = polynomial_13m(x, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13);

    // multiply by power of 2
    n2 = vm_pow2n(r);
    exp_x = (z + 1.0) * n2;
    n2 *= cosb;
    z = mul_add(z, n2, n2 - 1.0);            // z = z * n2 + (n2 - 1.0);
#ifdef SIGNED_ZERO                           // pedantic preservation of signed zero
    z = select(initial_x == 0., initial_x, z);
#endif

    // check for overflow
    auto inrange = abs(initial_x) < max_x;      // boolean vector
    // check for INF and NAN
    inrange &= is_finite(initial_x);

    if (horizontal_and(inrange)) {
        // fast normal path
        return z;
    }
    else {
        // overflow, underflow and NAN
        r = select(sign_bit(initial_x), 0. - (1 & 1), infinite_vec<VTYPE>()); // value in case of +/- overflow or INF
        z = select(inrange, z, r);                         // +/- underflow
        z = select(is_nan(initial_x), initial_x, z);       // NAN goes through

        r = select(sign_bit(initial_x), 0. - (0 & 1), infinite_vec<VTYPE>()); // value in case of +/- overflow or INF
        exp_x = select(inrange, exp_x, r);                         // +/- underflow
        exp_x = select(is_nan(initial_x), initial_x, exp_x);       // NAN goes through
        return z;
    }
}

template<typename VTYPE>
static inline VTYPE cosb_expm1_f(VTYPE const initial_x, VTYPE const cosb, VTYPE &exp_x)
{

    // Taylor coefficients
    const float P0expf = 1.f / 2.f;
    const float P1expf = 1.f / 6.f;
    const float P2expf = 1.f / 24.f;
    const float P3expf = 1.f / 120.f;
    const float P4expf = 1.f / 720.f;
    const float P5expf = 1.f / 5040.f;

    VTYPE x, r, x2, z, n2;                      // data vectors

    // maximum abs(x), value depends on BA, defined below
    // The lower limit of x is slightly more restrictive than the upper limit.
    // We are specifying the lower limit, except for BA = 1 because it is not used for negative x
    float max_x;
    const float ln2f_hi = 0.693359375f;
    const float ln2f_lo = -2.12194440e-4f;
    max_x = 87.3f;

    x = initial_x;
    r = round(initial_x * float(VM_LOG2E));
    x = nmul_add(r, VTYPE(ln2f_hi), x);      //  x -= r * ln2f_hi;
    x = nmul_add(r, VTYPE(ln2f_lo), x);      //  x -= r * ln2f_lo;

    x2 = x * x;
    z = polynomial_5(x, P0expf, P1expf, P2expf, P3expf, P4expf, P5expf);
    z = mul_add(z, x2, x);                       // z *= x2;  z += x;

    // multiply by power of 2
    n2 = vm_pow2n(r);
    exp_x = (z + 1.0f) * n2;
    n2 *= cosb;
    z = mul_add(z, n2, n2 - 1.0f);           //  z = z * n2 + (n2 - 1.0f);
#ifdef SIGNED_ZERO                           // pedantic preservation of signed zero
    z = select(initial_x == 0.f, initial_x, z);
#endif

    // check for overflow
    auto inrange = abs(initial_x) < max_x;      // boolean vector
    // check for INF and NAN
    inrange &= is_finite(initial_x);

    if (horizontal_and(inrange)) {
        // fast normal path
        return z;
    }
    else {
        // overflow, underflow and NAN
        r = select(sign_bit(initial_x), 0.f - (1 & 1), infinite_vec<VTYPE>()); // value in case of +/- overflow or INF
        z = select(inrange, z, r);                         // +/- underflow
        z = select(is_nan(initial_x), initial_x, z);       // NAN goes through

        r = select(sign_bit(initial_x), 0.f - (0 & 1), infinite_vec<VTYPE>()); // value in case of +/- overflow or INF
        exp_x = select(inrange, exp_x, r);                         // +/- underflow
        exp_x = select(is_nan(initial_x), initial_x, exp_x);       // NAN goes through
        return z;
    }
}

template <typename CVTYPE>
static inline CVTYPE expm1(CVTYPE z)
{
    static_assert(is_complex_vector<CVTYPE>());
    using VTYPE = decltype(z.to_vector());
    using RTYPE = decltype(VTYPE()[0]);

    VTYPE re, im, vec, s, c;

    vec = z.to_vector();

    // get real and imaginary parts in separate vectors
    if constexpr (vec.size() == 2) {
        re = permute2<0, 0>(vec);
        im = permute2<1, 1>(vec);
    }
    else if constexpr (vec.size() == 4) {
        re = permute4<0, 0, 2, 2>(vec);
        im = permute4<1, 1, 3, 3>(vec);
    }
    else if constexpr (vec.size() == 8) {
        re = permute8<0, 0, 2, 2, 4, 4, 6, 6>(vec);
        im = permute8<1, 1, 3, 3, 5, 5, 7, 7>(vec);
    }
    else if constexpr (vec.size() == 16) {
        re = permute16<0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12, 14, 14>(vec);
        im = permute16<1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15>(vec);
    }

    s = sincos(&c, im);
    if constexpr (std::is_same<RTYPE, double>::value) {
        // im = exp(a), vec = exp(a) cos(b) - 1
        vec = cosb_expm1_d(re, c, im);
    }
    else if constexpr (std::is_same<RTYPE, float>::value) {
        // im = exp(a), vec = exp(a) cos(b) - 1
        vec = cosb_expm1_f(re, c, im);
    }
    else {
        return 0;
    }
    im *= s;

    if constexpr (vec.size() == 2) {
        vec = blend2<0, 3>(vec, im);
    }
    else if constexpr (vec.size() == 4) {
        vec = blend4<0, 5, 2, 7>(vec, im);
    }
    else if constexpr (vec.size() == 8) {
        vec = blend8<0, 9, 2, 11, 4, 13, 6, 15>(vec, im);
    }
    else if constexpr (vec.size() == 16) {
        vec = blend16<0, 17, 2, 19, 4, 21, 6, 23, 8, 25, 10, 27, 12, 29, 14, 31>(vec, im);
    }

    return CVTYPE(vec);
}

#endif // VCL_COMPLEX_TRIG_H