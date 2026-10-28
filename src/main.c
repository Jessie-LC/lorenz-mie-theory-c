#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mie.h"

static inline double NormalDistribution(double x, double sigma, double mean) {
    return (1.0 / sqrt(2.0 * M_PI * (sigma*sigma))) * exp(-((pow(x - mean, 2.0)) / (2.0 * sigma*sigma)));
}
static inline double LogNormalDistribution(double x, double sigma, double mean) {
    return (1.0 / (x * sigma * sqrt(2.0 * M_PI))) * exp(-((pow(log(x) - mean, 2.0)) / (2.0 * sigma*sigma)));
}
static inline double WendlandKernelC6(const double dist, const double radius) {
    const double normalization = (78.0 / (7.0 * M_PI)) / pow(radius, 2.0); 
    const double q = dist / radius;
    if (q > 1.0) {
        return 0.0;
    }
    if (q <= 0.0) {
        return normalization;
    }
    const double  oneMinus_q = 1.0 - q;
    const double qPolynomial = 32.0 * pow(q, 3.0) + 25.0 * pow(q, 2.0) + 8.0 * q + 1.0;
    const double      kernel = fmax(pow(oneMinus_q, 8.0), 0.0) * qPolynomial;
    return normalization * kernel;
}

typedef struct ParticleDistribution {
    double* radii;
    double* numberDensity;

    fcomplex64_t refractiveIndex;

    uint32_t binCount;
} pdistribution_t;

static inline pdistribution_t CreateNormalDistribution(
    double minimumRadius, 
    double maximumRadius, 
    double radiusStep, 

    double Vdesired, 

    double sigma, 
    double mean 
) {
    pdistribution_t dist;
    dist.binCount      = (uint32_t)((abs((maximumRadius * 1e6) - (minimumRadius * 1e6)) / (radiusStep * 1e6)));
    dist.radii         = (double*)malloc(sizeof(double[dist.binCount]));
    dist.numberDensity = (double*)malloc(sizeof(double[dist.binCount]));
    double V = 0.0;
    for(uint32_t i = 0u; i < dist.binCount; ++i) {
        double       radius = ((double)i / (double)(dist.binCount-1u)) * maximumRadius + minimumRadius;
        double distribution = NormalDistribution(radius, sigma/1e6, mean/1e6);
        dist.radii[i]         = radius;
        dist.numberDensity[i] = distribution;
        V += pow(radius, 3.0) * distribution * radiusStep;
    }
    V = V * ((4.0 * M_PI) / 3.0);
    for(uint32_t i = 0u; i < dist.binCount; ++i) {
        dist.numberDensity[i] = dist.numberDensity[i] * (Vdesired / V);
    }

    return dist;
}
static inline pdistribution_t CreateLogNormalDistribution(
    double minimumRadius, 
    double maximumRadius, 
    double radiusStep, 

    double Vdesired, 

    double sigma, 
    double mean 
) {
    pdistribution_t dist;
    dist.binCount      = (uint32_t)((abs((maximumRadius * 1e6) - (minimumRadius * 1e6)) / (radiusStep * 1e6)));
    dist.radii         = (double*)malloc(sizeof(double[dist.binCount]));
    dist.numberDensity = (double*)malloc(sizeof(double[dist.binCount]));
    double V = 0.0;
    for(uint32_t i = 0u; i < dist.binCount; ++i) {
        double       radius = ((double)i / (double)(dist.binCount-1u)) * maximumRadius + minimumRadius;
        double distribution = LogNormalDistribution(radius*1e6, sigma, mean);
        dist.radii[i]         = radius;
        dist.numberDensity[i] = distribution;
        V += pow(radius, 3.0) * distribution * radiusStep;
    }
    V = V * ((4.0 * M_PI) / 3.0);
    for(uint32_t i = 0u; i < dist.binCount; ++i) {
        dist.numberDensity[i] = dist.numberDensity[i] * (Vdesired / V);
    }

    return dist;
}
static inline pdistribution_t CreateWendlandDistribution(
    double minimumRadius, 
    double maximumRadius, 
    double radiusStep, 

    double Vdesired, 

    double sigma, 
    double mean 
) {
    pdistribution_t dist;
    dist.binCount      = (uint32_t)((abs((maximumRadius * 1e6) - (minimumRadius * 1e6)) / (radiusStep * 1e6)));
    dist.radii         = (double*)malloc(sizeof(double[dist.binCount]));
    dist.numberDensity = (double*)malloc(sizeof(double[dist.binCount]));
    double V = 0.0;
    for(uint32_t i = 0u; i < dist.binCount; ++i) {
        double       radius = ((double)i / (double)(dist.binCount-1u)) * maximumRadius + minimumRadius;
        double distribution = WendlandKernelC6(sqrt(pow((radius*1e6) - mean, 2.0)), sigma);
        dist.radii[i]         = radius;
        dist.numberDensity[i] = distribution;
        V += pow(radius, 3.0) * distribution * radiusStep;
    }
    V = V * ((4.0 * M_PI) / 3.0);
    for(uint32_t i = 0u; i < dist.binCount; ++i) {
        dist.numberDensity[i] = dist.numberDensity[i] * (Vdesired / V);
    }

    return dist;
}
static inline void DeleteDistribution(pdistribution_t* dist) {
    free(dist->radii);
    free(dist->numberDensity);
}

int main(int argv, char** argc) {
    fcomplex64_t particle = { 1.3333, 1e-4 };
    fcomplex64_t     host = { 1.00028, 0.0 };
    double         lambda = 500e-9;
    double  minimumRadius = 1e-7;
    double  maximumRadius = 1e-5;
    double     radiusStep = 1e-7;
    double   waterDensity = 1000.0;
    double    waterWeight = 0.0001;
    double    waterVolume = waterWeight < 1e-12 ? 0.0 : waterWeight / waterDensity;
    double      airVolume = 1.0 - waterVolume;
    pdistribution_t  dist = CreateWendlandDistribution(
        minimumRadius, 
        maximumRadius, 
        radiusStep, 

        (waterVolume / airVolume), 

        4.0, 
        0.8 
    );

    const uint32_t      angles = 720u;
    const uint32_t wavelengths = 471u/4u;
    FILE*            phaseFile = fopen("phase.dat", "w+b");
    double*      phaseFunction = (double*)malloc(sizeof(double)*angles*wavelengths);
    for(uint32_t angle = 0u; angle < angles; ++angle) {
        double dtheta = M_PI / (double)(angles - 1u);
        double theta = (double)angle * dtheta;

        for(uint32_t wavelength = 0u; wavelength < wavelengths; ++wavelength) {
            double             lambda = (((double)wavelength / (double)wavelengths) * 471.0 + 360.0) * 1e-9;
            double ensembleScattering = 0.0;
            double ensembleExtinction = 0.0;
            double      ensemblePhase = 0.0;
            for(uint32_t i = 0u; i < dist.binCount; ++i) {
                if (dist.numberDensity[i] == 0.0) {
                    // Literally no particles exist with this radius
                    continue;
                }
                double scattering;
                double extinction;
                double s_polarized;
                double p_polarized;
                double unpolarized;
                CalculateLorenzMieTheory(
                    theta, 
                    lambda, 
                    dist.radii[i], 
                    host, 
                    particle, 

                    &scattering, 
                    &extinction, 

                    &s_polarized, 
                    &p_polarized, 
                    &unpolarized 
                );

                ensembleScattering += dist.numberDensity[i] * scattering * radiusStep;
                ensemblePhase += dist.numberDensity[i] * scattering * unpolarized * radiusStep;
            }
            phaseFunction[wavelength+wavelengths*angle] = ensemblePhase;
            char*          value = (char*)malloc(sizeof(float));
            float fEnsemblePhase = (float)ensemblePhase;
            memcpy(&value[0], &fEnsemblePhase, sizeof(float));
            fwrite(&value[0], sizeof(float), 1u, phaseFile);
            printf(" %f\t\r", ((double)angle/(double)angles)*100.0);
        }
    }
    DeleteDistribution(&dist);

    // fwrite((void*)phaseFunction, sizeof(double), angles*wavelengths, phaseFile);
    fclose(phaseFile);
    free(phaseFunction);

    return 0;
}
