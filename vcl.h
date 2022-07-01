#ifndef VCL_VCL_H
#define VCL_VCL_H

#include <vcl/vectorclass.h>
#include <vcl/vectormath_trig.h>
#include <vcl/vectormath_exp.h>
#include <vcl/vectormath_hyp.h>
#include <vcl/complexvec_trig.h>
#include <ostream>


#define DECL_VEC_PRINT(VTYPE) \
    std::ostream &operator<<(std::ostream &o, const VTYPE &v)       \
    {                                                               \
        o << '[';                                                   \
        for (size_t i = 0; i < VTYPE::size(); i++) {                \
            o << v[i] << (i != (VTYPE::size() - 1) ? ", " : "]\n"); \
        }                                                           \
        return o;                                                   \
    } 

#define DECL_VCL_FMOD(VTYPE) \
    static inline VTYPE fmod(VTYPE x, VTYPE y) { return x - (y * truncate(x / y)); }

#define DECL_VCL_HYPOT(VTYPE) \
static inline VTYPE hypot(VTYPE x, VTYPE y) \
{                                           \
    VTYPE m;                                \
    x = abs(x);                             \
    y = abs(y);                             \
    m = max(x, y);                          \
    x = min(x, y) / m;                      \
    return m * sqrt(mul_add(x, x, 1.0));    \
}

static inline Vec2d rint(Vec2d x)
{
    return _mm_round_pd(x, 4);
}

static inline Vec4f rint(Vec4f x)
{
    return _mm_round_ps(x, 4);
}

#if MAX_VECTOR_SIZE >= 256

static inline Vec4d rint(Vec4d x)
{
    return _mm256_round_pd(x, 4);
}

static inline Vec8f rint(Vec8f x)
{
    return _mm256_round_ps(x, 4);
}

#if MAX_VECTOR_SIZE >= 512

static inline Vec8d rint(Vec8d x)
{
#if INSTRSET >= 8
    return Vec8d(rint(x.get_low()), rint(x.get_high()));
#else
    return _mm512_roundscale_pd((__m512d)x, 4);
#endif // INSTRSET >= 8
}

static inline Vec16f rint(Vec16f x)
{
#if INSTRSET >= 8
    return Vec16f(rint(x.get_low()), rint(x.get_high()));
#else
    return _mm512_roundscale_ps((__m512d)x, 4);
#endif // INSTRSET >= 8
}

#endif // MAX_VECTOR_SIZE >= 512

#endif // MAX_VECTOR_SIZE >= 256

template <typename T>
static inline T rint(T z)
{
    static_assert(is_complex_vector<T>());
    return rint(z.to_vector());
}

template <typename T>
static inline T round(T z)
{
    static_assert(is_complex_vector<T>());
    return round(z.to_vector());
}

template <typename T>
static inline T floor(T z)
{
    static_assert(is_complex_vector<T>());
    return floor(z.to_vector());
}

template <typename T>
static inline T ceil(T z)
{
    static_assert(is_complex_vector<T>());
    return ceil(z.to_vector());
}

template <typename T>
static inline T truncate(T z)
{
    static_assert(is_complex_vector<T>());
    return truncate(z.to_vector());
}


DECL_VEC_PRINT(Vec2d)
DECL_VEC_PRINT(Vec4d)
DECL_VEC_PRINT(Vec8d)
DECL_VEC_PRINT(Vec4f)
DECL_VEC_PRINT(Vec8f)
DECL_VEC_PRINT(Vec16f)

DECL_VEC_PRINT(Vec2db)
DECL_VEC_PRINT(Vec4db)
DECL_VEC_PRINT(Vec8db)
DECL_VEC_PRINT(Vec4fb)
DECL_VEC_PRINT(Vec8fb)
DECL_VEC_PRINT(Vec16fb)

DECL_VEC_PRINT(Vec16c)
DECL_VEC_PRINT(Vec16cb)
DECL_VEC_PRINT(Vec16uc)

DECL_VEC_PRINT(Vec8s)
DECL_VEC_PRINT(Vec8sb)
DECL_VEC_PRINT(Vec8us)

DECL_VEC_PRINT(Vec4i)
DECL_VEC_PRINT(Vec4ib)
DECL_VEC_PRINT(Vec4ui)

DECL_VEC_PRINT(Vec2q)
DECL_VEC_PRINT(Vec2qb)
DECL_VEC_PRINT(Vec2uq)

DECL_VEC_PRINT(Vec32c)
DECL_VEC_PRINT(Vec32cb)
DECL_VEC_PRINT(Vec32uc)

DECL_VEC_PRINT(Vec16s)
DECL_VEC_PRINT(Vec16sb)
DECL_VEC_PRINT(Vec16us)

DECL_VEC_PRINT(Vec8i)
DECL_VEC_PRINT(Vec8ib)
DECL_VEC_PRINT(Vec8ui)

DECL_VEC_PRINT(Vec4q)
DECL_VEC_PRINT(Vec4qb)
DECL_VEC_PRINT(Vec4uq)

DECL_VEC_PRINT(Vec64c)
DECL_VEC_PRINT(Vec64cb)
DECL_VEC_PRINT(Vec64uc)

DECL_VEC_PRINT(Vec32s)
DECL_VEC_PRINT(Vec32sb)
DECL_VEC_PRINT(Vec32us)

DECL_VEC_PRINT(Vec16i)
DECL_VEC_PRINT(Vec16ib)
DECL_VEC_PRINT(Vec16ui)

DECL_VEC_PRINT(Vec8q)
DECL_VEC_PRINT(Vec8qb)
DECL_VEC_PRINT(Vec8uq)

DECL_VCL_FMOD(Vec2d)
DECL_VCL_FMOD(Vec4d)
DECL_VCL_FMOD(Vec8d)
DECL_VCL_FMOD(Vec4f)
DECL_VCL_FMOD(Vec8f)
DECL_VCL_FMOD(Vec16f)

DECL_VCL_HYPOT(Vec2d)
DECL_VCL_HYPOT(Vec4d)
DECL_VCL_HYPOT(Vec8d)
DECL_VCL_HYPOT(Vec4f)
DECL_VCL_HYPOT(Vec8f)
DECL_VCL_HYPOT(Vec16f)

#endif // VCL_VCL_H