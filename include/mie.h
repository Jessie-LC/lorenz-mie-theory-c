#if !defined _MIE_THEORY
#define _MIE_THEORY

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#include "./math/math.h"

void CalculateLorenzMieTheory(
    double theta, 
    double lambda, 
    double radius, 
    fcomplex64_t nHost, 
    fcomplex64_t nParticle, 

    double* extinction, 
    double* scattering, 

    double* s_polarized, 
    double* p_polarized, 
    double* unpolarized 
) {
    fcomplex64_t k = Complex64_DivS(Complex64_MulS(nHost, 2.0 * M_PI), lambda);
    fcomplex64_t x = Complex64_DivS(Complex64_MulS(nHost, 2.0 * M_PI * radius), lambda);
    fcomplex64_t y = Complex64_DivS(Complex64_MulS(nParticle, 2.0 * M_PI * radius), lambda);
    double    size = Complex64_Abs(x);

    uint32_t M = (uint32_t)(size + 4.3 * pow(size, 1.0/3.0) + 1.0);

    double*              piN = (double*)malloc(sizeof(double[M]));
    double*             tauN = (double*)malloc(sizeof(double[M]));
    fcomplex64_t* AparticleN = (fcomplex64_t*)malloc(sizeof(fcomplex64_t[M+2u]));
    fcomplex64_t*     AhostN = (fcomplex64_t*)malloc(sizeof(fcomplex64_t[M+2u]));
    fcomplex64_t*         aN = (fcomplex64_t*)malloc(sizeof(fcomplex64_t[M+1u]));
    fcomplex64_t*         bN = (fcomplex64_t*)malloc(sizeof(fcomplex64_t[M+1u]));
    
    for(uint32_t n = 0; n <= M+1u; ++n) {
        AparticleN[n] = SetComplex64(0.0, 0.0);
        AhostN[n]     = SetComplex64(0.0, 0.0);
    }

    for(uint32_t n = 1u; n < M; ++n) {
        legendrepolynomial_t legendre = CalculatePnAndDerivatives(n, cos(theta));
        piN[n]  = legendre.Pn_d1;
        tauN[n] = cos(theta) * legendre.Pn_d1 - sin(theta)*sin(theta) * legendre.Pn_d2;
    }

    AparticleN[M+1u] = SetComplex64(1.0, 0.0);
    AhostN[M+1u]     = SetComplex64(1.0, 0.0);
    for(uint32_t n = M; n <= M; --n) {
        fcomplex64_t tmpH = Complex64_DivC(SetComplex64((double)n + 1.0, 0.0), x);
        fcomplex64_t tmpP = Complex64_DivC(SetComplex64((double)n + 1.0, 0.0), y);
        AhostN[n]     = Complex64_SubC(tmpH, Complex64_DivC(SetComplex64(1.0, 0.0), Complex64_AddC(tmpH, AhostN[n+1u])));
        AparticleN[n] = Complex64_SubC(tmpP, Complex64_DivC(SetComplex64(1.0, 0.0), Complex64_AddC(tmpP, AparticleN[n+1u])));
    }

    fcomplex64_t           B = { 0.0, 1.0 };
    fcomplex64_t     psiZeta = Complex64_MulC(SetComplex64(0.5, 0.0), Complex64_SubC(SetComplex64(1.0, 0.0), Complex64_Exp(SetComplex64(-2.0 * x.val[1],  2.0 * x.val[0]))));
    fcomplex64_t psiOverZeta = Complex64_MulC(SetComplex64(0.5, 0.0), Complex64_SubC(SetComplex64(1.0, 0.0), Complex64_Exp(SetComplex64( 2.0 * x.val[1], -2.0 * x.val[0]))));
    for(uint32_t n = 1u; n < M; ++n) {
        fcomplex64_t nZ = Complex64_DivC(SetComplex64((double)n, 0.0), x);
        psiZeta     = Complex64_MulC(psiZeta, Complex64_MulC(Complex64_SubC(nZ, AhostN[n-1u]), Complex64_SubC(nZ, B)));
        B           = Complex64_AddC(AhostN[n], Complex64_DivC(SetComplex64(0.0, 1.0), psiZeta));
        psiOverZeta = Complex64_MulC(psiOverZeta, Complex64_DivC(Complex64_AddC(B, nZ), Complex64_AddC(AhostN[n], nZ)));

        aN[n] = Complex64_MulC(psiOverZeta, Complex64_DivC(
            Complex64_SubC(Complex64_MulC(nHost, AparticleN[n]), Complex64_MulC(nParticle, AhostN[n])), 
            Complex64_SubC(Complex64_MulC(nHost, AparticleN[n]), Complex64_MulC(nParticle, B)) 
        ));
        bN[n] = Complex64_MulC(psiOverZeta, Complex64_DivC(
            Complex64_SubC(Complex64_MulC(nParticle, AparticleN[n]), Complex64_MulC(nHost, AhostN[n])), 
            Complex64_SubC(Complex64_MulC(nParticle, AparticleN[n]), Complex64_MulC(nHost, B)) 
        ));
    }

    fcomplex64_t S1 = { 0.0, 0.0 };
    fcomplex64_t S2 = { 0.0, 0.0 };
    fcomplex64_t  S = { 0.0, 0.0 };
    double cscattering = 0.0;
    for(uint32_t n = 1u; n < M; ++n) {
        S = Complex64_AddC(S, Complex64_MulS(Complex64_AddC(aN[n], bN[n]), (double)(2 * n + 1)));
        double term = (double)(2 * n + 1) / (double)(n * n + 1);
        S1 = Complex64_AddC(S1, Complex64_MulS(Complex64_AddC(Complex64_MulS(aN[n], piN[n]), Complex64_MulS(bN[n], tauN[n])), term));
        S2 = Complex64_AddC(S2, Complex64_MulS(Complex64_AddC(Complex64_MulS(aN[n], tauN[n]), Complex64_MulS(bN[n], piN[n])), term));
        cscattering += (double)(2 * n + 1) * ((Complex64_Abs(aN[n])*Complex64_Abs(aN[n])) + (Complex64_Abs(bN[n])*Complex64_Abs(bN[n])));
    }
    S = Complex64_MulS(S, 0.5);
    fcomplex64_t tmpExt = Complex64_DivC(S, Complex64_MulC(k, k));
    double cextinction = 4.0 * M_PI * tmpExt.val[0];

    double nkhSquared = Complex64_Abs(nHost)*Complex64_Abs(nHost);
    double      alpha = 4.0 * M_PI * radius * nHost.val[1] / lambda;
    double      gamma = alpha <= 1e-6 ? 1.0 : (2.0 * (1.0 + (alpha - 1.0) * exp(alpha))) / (alpha*alpha);
    cscattering = (((lambda*lambda) * exp(-alpha)) / (2.0 * M_PI * gamma * nkhSquared)) * cscattering;

    double absS1 = Complex64_Abs(S1);
    double absS2 = Complex64_Abs(S2);
    double phase_unpolarized = (absS1*absS1 + absS2*absS2) / (2.0 * (Complex64_Abs(k)*Complex64_Abs(k)) * cscattering);
    double phase_s_polarized = (absS1*absS1) / ((Complex64_Abs(k)*Complex64_Abs(k)) * cscattering);
    double phase_p_polarized = (absS2*absS2) / ((Complex64_Abs(k)*Complex64_Abs(k)) * cscattering);

    if (extinction != NULL) {
        *extinction = cextinction;
    }
    if (scattering != NULL) {
        *scattering = cscattering;
    }
    if (s_polarized != NULL) {
        *s_polarized = phase_s_polarized;
    }
    if (p_polarized != NULL) {
        *p_polarized = phase_p_polarized;
    }
    if (unpolarized != NULL) {
        *unpolarized = phase_unpolarized;
    }

    free(piN);
    free(tauN);
    free(AparticleN);
    free(AhostN);
    free(aN);
    free(bN);
}

#endif
