//
// Created by EYx on 2023/4/2.
//

#ifndef FRNS_BOOTSTRAP_MYCHEBYSHEVFUNC_H
#define FRNS_BOOTSTRAP_MYCHEBYSHEVFUNC_H

#include "Common.h"
#include "Ciphertext.h"
#include "Scheme.h"


/**
 * Computes the quotient and remainder of the long division of two polynomials in the Chebyshev series basis
 *
 * @param &f the vector of coefficients of the dividend.
 * @param &g the vector of coefficients of the divisor.
 * @return a struct with the coefficients for the quotient and remainder.
 */
//std::shared_ptr<longDiv> LongDivisionChebyshev(const std::vector<double>& f, const std::vector<double>& g);

//implement function evaluation relative procedure
void EvalChebyshevSeries(Ciphertext& result, Ciphertext x, double* coefficients,
                               double a, double b, uint32_t poly_degree,
                               Scheme scheme);
void EvalChebyshevSeriesPS(Ciphertext& result, Ciphertext x, double* coefficients,
                                 double a, double b, uint32_t coefficients_len,
                                 Scheme scheme);

/**
 * Gets the degree of a polynomial specified by its coefficients.
 *
 * @param &coefficients vector of coefficients of a polynomial.
 * @return the integer degree of the polynomial.
 */
//uint32_t Degree(const std::vector<double>& coefficients);
uint32_t Degree(double *coefficients, uint32_t poly_degree);

/* Compute positive integers k,m such that n < k(2^m-1) and k close to sqrt(n/2) */
std::vector<uint32_t> ComputeDegreesPS(const uint32_t n);

void CheckAndAdjusLevel(Ciphertext& rct1, Ciphertext& rct2, Ciphertext& ct1, Ciphertext& ct2, Scheme scheme);

/**
 * Computes the quotient and remainder of the long division of two polynomials in the Chebyshev series basis
 *
 * @param &f the vector of coefficients of the dividend.
 * @param &g the vector of coefficients of the divisor.
 * @return a struct with the coefficients for the quotient and remainder.
 */
void LongDivisionChebyshev(double*& quotient,uint32_t & quotient_len,
                           double*& remainder,uint32_t & remainder_len,
                           double* f, uint32_t f_len,
                           double* g, uint32_t g_len
);

inline bool IsNotEqualOne(double val);

void EvalLinearWSumMutable(Ciphertext & wsum,
                           Ciphertext * ciphertexts, uint32_t ciphertexts_num, double * constants,
                           Scheme scheme);

void InnerEvalChebyshevPS(Ciphertext& result, const Ciphertext x,
                                double* coefficients, uint32_t coefficients_len,
                                uint32_t k, uint32_t m,
                                Ciphertext* T,Ciphertext*T2,
                                Scheme scheme);

#endif //FRNS_BOOTSTRAP_MYCHEBYSHEVFUNC_H
