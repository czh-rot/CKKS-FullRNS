//
// Created by EYx on 2023/4/29.
//

#ifndef FRNS_BOOTSTRAP_MYUTILS_H
#define FRNS_BOOTSTRAP_MYUTILS_H

#include "Ciphertext.h"
#include "Scheme.h"

void SetZero(uint64_t * a, long Start, long End);

void CopyValue(uint64_t * res, uint64_t * input, long Start, long End);

/**
* Converts the value to an double.
*
* @return double representation of the value.
*/
double ConvertToDouble(uint64_t input);

//------------------------------------------------------------------------------
// Running Approximate Mod Reduction
//------------------------------------------------------------------------------

// Double-angle iterations are applied in the case of OPTIMIZED/uniform secrets
void ApplyDoubleAngleIterations(Ciphertext& ciphertext,Scheme scheme);


#endif //FRNS_BOOTSTRAP_MYUTILS_H
