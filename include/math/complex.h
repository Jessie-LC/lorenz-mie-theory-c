#if !defined _COMPLEX_MATH
#define _COMPLEX_MATH

#include <immintrin.h>
#include <math.h>

typedef struct ComplexDouble {
    alignas(16) double val[2];
} fcomplex64_t;

static inline fcomplex64_t SetComplex64(double real, double imag) {
    fcomplex64_t complex = { real, imag };
    return complex;
}

static inline fcomplex64_t Complex64_AddC(fcomplex64_t lhs, fcomplex64_t rhs) {
    __m128d op = _mm_add_pd(
        _mm_load_pd(&lhs.val[0]), 
        _mm_load_pd(&rhs.val[0]) 
    );
    fcomplex64_t result;
    _mm_store_pd(&result.val[0], op);
    return result;
}
static inline fcomplex64_t Complex64_AddS(fcomplex64_t lhs, double rhs) {
    fcomplex64_t result;
    result.val[0] = lhs.val[0] + rhs;
    result.val[1] = lhs.val[1];
    return result;
}

static inline fcomplex64_t Complex64_SubC(fcomplex64_t lhs, fcomplex64_t rhs) {
    __m128d op = _mm_sub_pd(
        _mm_load_pd(&lhs.val[0]), 
        _mm_load_pd(&rhs.val[0]) 
    );
    fcomplex64_t result;
    _mm_store_pd(&result.val[0], op);
    return result;
}
static inline fcomplex64_t Complex64_SubS(fcomplex64_t lhs, double rhs) {
    fcomplex64_t result;
    result.val[0] = lhs.val[0] - rhs;
    result.val[1] = lhs.val[1];
    return result;
}

static inline fcomplex64_t Complex64_MulC(fcomplex64_t lhs, fcomplex64_t rhs) {
    __m128d opLhs = _mm_mul_pd(
        _mm_load_pd(&lhs.val[0]), 
        _mm_load_pd1(&rhs.val[0]) 
    );
    __m128d opRhs = _mm_mul_pd(
        _mm_set_pd(lhs.val[0], lhs.val[1]), 
        _mm_load_pd1(&rhs.val[1]) 
    );
    __m128d op = _mm_addsub_pd(opLhs, opRhs);

    fcomplex64_t result = { 0.0, 0.0 };
    _mm_store_pd(&result.val[0], op);
    return result;
}
static inline fcomplex64_t Complex64_MulS(fcomplex64_t lhs, double rhs) {
    fcomplex64_t result = {
        lhs.val[0] * rhs, 
        lhs.val[1] * rhs
    };
    return result;
}

static inline fcomplex64_t Complex64_DivC(fcomplex64_t lhs, fcomplex64_t rhs) {
    double div = rhs.val[0]*rhs.val[0] + rhs.val[1]*rhs.val[1];
    fcomplex64_t result = {
        (lhs.val[0]*rhs.val[0] + lhs.val[1]*rhs.val[1]) / div,
        (lhs.val[1]*rhs.val[0] - lhs.val[0]*rhs.val[1]) / div
    };
    return result;
}
static inline fcomplex64_t Complex64_DivS(fcomplex64_t lhs, double rhs) {
    fcomplex64_t result = {
        lhs.val[0] / rhs,
        lhs.val[1] / rhs
    };
    return result;
}

static inline fcomplex64_t Complex64_Conj(fcomplex64_t z) {
    fcomplex64_t result = {
         z.val[0],
        -z.val[1]
    };
    return result;
}

static inline double Complex64_Abs(fcomplex64_t z) {
    return sqrt(z.val[0]*z.val[0] + z.val[1]*z.val[1]);
}

static inline fcomplex64_t Complex64_Exp(fcomplex64_t z) {
    fcomplex64_t v = { cos(z.val[1]), sin(z.val[1]) };
    return Complex64_MulS(v, exp(z.val[0]));
}
static inline fcomplex64_t Complex64_Log(fcomplex64_t z) {
    fcomplex64_t result = {
        0.5 * log(z.val[0]*z.val[0] + z.val[1]*z.val[1]), 
        atan2(z.val[1], z.val[0]) 
    };
    return result;
}

#endif
