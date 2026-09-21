#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include "mie.h"

int main(int argv, char** argc) {
    fcomplex64_t particle = { 1.33333, 1e-4 };
    fcomplex64_t     host = { 1.00028, 0.0 };
    double         lambda = 550e-9;
    double         radius = 1e-6;
    for(uint32_t angle = 0u; angle < 180u; ++angle) {
        double dtheta = M_PI / (double)(180u - 1u);
        double theta = (double)angle * dtheta;

        double scattering;
        double extinction;
        double s_polarized;
        double p_polarized;
        double unpolarized;
        CalculateLorenzMieTheory(
            theta, 
            lambda, 
            radius, 
            host, 
            particle, 

            &scattering, 
            &extinction, 

            &s_polarized, 
            &p_polarized, 
            &unpolarized 
        );

        printf("%f;%f\n", theta, log(unpolarized));
    }

    return 0;
}
