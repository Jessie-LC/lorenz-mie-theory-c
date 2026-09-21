#if !defined _MATH_LEGENDRE
#define _MATH_LEGENDRE

typedef struct LegendrePolynomial {
    double Pn;
    double Pn_d1;
    double Pn_d2;
} legendrepolynomial_t;

static inline legendrepolynomial_t CalculatePnAndDerivatives(uint32_t n, double x) {
    /*

    I will explain the derivation of the first derivative of the Legendre polynomial with respect to x. 
    The second derivative follows the same basic steps, but takes way longer to write, so I have decided not to include the derivation for it.


    Split this into two halves, one before the minus sign and one after
    (2 * n - 1) * (x * f(n-1,x)) - (n - 1) * f(n-2,x)

    f1(x) = f(n-1,x)
    f2(x) = f(n-2,x)

    First half : (2 * n - 1) * (x * f1(x))
    Second half : (n - 1) * f2(x)

    Find the derivative of the first half using the product rule. The constant rule applies here, (2 * n - 1) is a constant, so this is the special case 
    where (a*f1(x))' = a*f1'(x). As such we can ignore the first half of the product rule calculate, and only do the second half.
    But let's assume for now that a = 1, simply because that makes explaining the derivation simpler.
    We can expand f1'(x) into (x * f1(x))', we then substitute x with function u(x) = x, making f1'(x) = (u(x) * f1(x))', which we then expand into 
    the final algorithm using the product rule.
    (u(x) * f1(x))' = u'(x) * f1(x) + u(x) * f1'(x)
    We now expand u'(x)
    u'(x) * f(n-1,x) = 0 = (1 * x^(1-1)) * f(n-1,x) = 0
    Replace u(x) with x.
    The derivative of the first half with respect to x, once we factor a back in and replace u(x) with x, is: (2 * n - 1) * (x * f1'(x))

    Now we find the derivative of the second half using the product rule. The constant rule also applies here.
    I was going to write a detailed explanation of how I got to the final algorithm for the second half like I did for the first half, 
    but honestly it is literally just a*f2'(x). It's not like the first half where f1'(x) expands into some new algorithm. 
    Finding the derivative of (n - 1) * f2(x), once the constant rule is factored in, is as simple as saying 
    "Oh, this is just the case of (a * f(x))' = a * f'(x)".


    We then use the difference/sum rule to get our first derivative.
    (2 * n - 1) * (x * f1'(x)) - (n - 1) * f2'(x)

    */
    if (n == 0) {
        legendrepolynomial_t p0 = {1.0, 0.0, 0.0};
        return p0;
    } else if (n == 1) {
        legendrepolynomial_t p1 = {x, 1.0, 0.0};
        return p1;
    }

    double Pn_minus1 = x;
    double Pn_minus2 = 1.0;
    double PD1n_minus1 = 1.0;
    double PD1n_minus2 = 0.0;
    double PD2n_minus1 = 0.0;
    double PD2n_minus2 = 0.0;
    for(uint32_t i = 2; i <= n; ++i) {
        double Pn_current = ((double)(2 * i - 1) * x * Pn_minus1 - ((double)i - 1) * Pn_minus2) / (double)i;
        double PD1n_current = ((double)(2 * i - 1) * (Pn_minus1 + x * PD1n_minus1) - ((double)i - 1) * PD1n_minus2) / (double)i;
        double PD2n_current = ((double)(2 * i - 1) * (PD1n_minus1 + (PD1n_minus1 + x * PD2n_minus1)) - ((double)i - 1) * PD2n_minus2) / (double)i;

        Pn_minus2 = Pn_minus1;
        Pn_minus1 = Pn_current;
        PD1n_minus2 = PD1n_minus1;
        PD1n_minus1 = PD1n_current;
        PD2n_minus2 = PD2n_minus1;
        PD2n_minus1 = PD2n_current;
    }

    legendrepolynomial_t pn = {Pn_minus1, PD1n_minus1, PD2n_minus1};

    return pn;
}

#endif
