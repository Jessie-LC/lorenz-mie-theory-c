#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mie.h"
#include "refractive-indices.h"

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

typedef fcomplex64_t(*iorfunction_t)(double lambda);

typedef struct ParticleDistribution {
    double* radii;
    double* numberDensity;

    iorfunction_t refractiveIndex;

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
    dist.refractiveIndex = NULL;
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
    dist.refractiveIndex = NULL;
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
    dist.refractiveIndex = NULL;
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
    double  minimumRadius = 1e-7;
    double  maximumRadius = 8e-6;
    double     radiusStep = 1e-7;
    double   waterDensity = 1000.0;
    double    waterWeight = 0.0001;
    double    waterVolume = waterWeight < 1e-12 ? 0.0 : waterWeight / waterDensity;
    double      airVolume = 1.0 - waterVolume;
    pdistribution_t  dist = CreateLogNormalDistribution(
        minimumRadius, 
        maximumRadius, 
        radiusStep, 

        (waterVolume / airVolume), 

        4.0, 
        0.8 
    );
    // dist.refractiveIndex = &Water_Complex;
    const uint32_t      angles = 720u;
    const uint32_t wavelengths = 471u;
    FILE*            phaseFile = fopen("phase.dat", "w+b");
    FILE*              cdfFile = fopen("cdf.dat", "w+b");
    double*      phaseFunction = (double*)malloc(sizeof(double)*angles*wavelengths);
    double*     anglesFunction = (double*)malloc(sizeof(double)*angles);
    for(uint32_t angle = 0u; angle < angles; ++angle) {
        double dtheta = M_PI / (double)(angles - 1u);
        double  theta = (double)angle * dtheta;
        anglesFunction[angle] = cos(theta);

        for(uint32_t wavelength = 0u; wavelength < wavelengths; ++wavelength) {
            double             lambda = (((double)wavelength / (double)wavelengths) * 471.0 + 360.0) * 1e-9;
            double ensembleScattering = 0.0;
            double ensembleExtinction = 0.0;
            double      ensemblePhase = 0.0;
            fcomplex64_t          air = { N2_Peck0C(lambda * 1e6), 0.0 };
            host = air;
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
                    dist.refractiveIndex == NULL ? particle : dist.refractiveIndex(lambda * 1e6), 

                    &scattering, 
                    &extinction, 

                    &s_polarized, 
                    &p_polarized, 
                    &unpolarized 
                );

                ensembleScattering += dist.numberDensity[i] * scattering * radiusStep;
                ensemblePhase += dist.numberDensity[i] * scattering * unpolarized * radiusStep;
            }
            phaseFunction[wavelength+wavelengths*angle] = (1.0 / ensembleScattering) * ensemblePhase;
            // char*          value = (char*)malloc(sizeof(float));
            // float fEnsemblePhase = (float)ensemblePhase;
            // memcpy(&value[0], &fEnsemblePhase, sizeof(float));
            // fwrite(&value[0], sizeof(float), 1u, phaseFile);
            // free(value);
            printf(" %f\t\r", ((double)angle/(double)angles)*100.0);
        }
    }
    DeleteDistribution(&dist);

    // CDF calculation from https://github.com/BelmuTM/Lorenz-Mie-Solver/blob/master/scripts/phase_cdf.py
    double* cdf = (double*)malloc(sizeof(double)*angles*wavelengths);
    for(uint32_t wavelength = 0u; wavelength < wavelengths; ++wavelength) {
        for(uint32_t angle = 1u; angle < angles; ++angle) {
            cdf[wavelength + angle * wavelengths] = 0.0f;
        }
    }
    for(uint32_t wavelength = 0u; wavelength < wavelengths; ++wavelength) {
        for(uint32_t angle = 1u; angle < angles; ++angle) {
            double p1 = phaseFunction[wavelength+wavelengths*(angle - 1)] + phaseFunction[wavelength+wavelengths*angle];
            double p2 = anglesFunction[angle-1] - anglesFunction[angle];
            cdf[wavelength + angle * wavelengths] = cdf[wavelength + (angle-1) * wavelengths] + 2.0 * M_PI * p1 * 0.5 * p2;
        }
    }
    for(uint32_t wavelength = 0u; wavelength < wavelengths; ++wavelength) {
        double maxCDF = cdf[wavelength+wavelengths*(angles-1)];
        for(uint32_t angle = 0u; angle < angles; ++angle) {
            phaseFunction[wavelength+wavelengths*angle] /= maxCDF;
            cdf[wavelength + angle * wavelengths]       /= maxCDF;
        }
    }
    for(uint32_t angle = 0u; angle < angles; ++angle) {
        for(uint32_t wavelength = 0u; wavelength < wavelengths; ++wavelength) {
            char* cdfValue = (char*)malloc(sizeof(float));
            float     fCdf = (float)cdf[wavelength + angle * wavelengths];
            memcpy(&cdfValue[0], &fCdf, sizeof(float));
            fwrite(&cdfValue[0], sizeof(float), 1u, cdfFile);
            free(cdfValue);

            char* phaseValue = (char*)malloc(sizeof(float));
            float     fPhase = (float)phaseFunction[wavelength + wavelengths * angle];
            memcpy(&phaseValue[0], &fPhase, sizeof(float));
            fwrite(&phaseValue[0], sizeof(float), 1u, phaseFile);
            free(phaseValue);
        }
    }
    fclose(phaseFile);
    fclose(cdfFile);
    free(phaseFunction);

    /*

            cdf[w][a] = (
                 cdf[w][a - 1]
               + 2.0 * np.pi * (phase[w][a - 1] + phase[w][a]) * 0.5
               * (angles[a - 1] - angles[a])
            )

    # Normalize both the phase and the CDF
    for w in range(0, wavelengths_n):

        max_cdf = cdf[w][-1]

        for a in range(0, angles_n):

            phase[w][a] /= max_cdf
            cdf[w][a]   /= max_cdf

    return phase, cdf

    */

    return 0;
}
