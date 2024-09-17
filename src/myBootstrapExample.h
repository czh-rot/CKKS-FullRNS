//
// Created by EYx on 2023/5/5.
//

#ifndef FRNS_BOOTSTRAP_MYBOOTSTRAPEXAMPLE_H
#define FRNS_BOOTSTRAP_MYBOOTSTRAPEXAMPLE_H

#include "Common.h"
#include "Scheme.h"

#include "constants.h"

void BootstrapExampleClean(uint32_t n=64, uint32_t slots=8, uint32_t levelsRemaining=2, SecretKeyDist secretKeyDist=UNIFORM_TERNARY);
uint32_t GetBootstrapDepth(uint32_t approxModDepth, uint32_t* levelBudget,
                           SecretKeyDist secretKeyDist=UNIFORM_TERNARY);

void FindBootstrapRotationIndices(Scheme scheme, int32_t* fullIndexList, long fullIndexList_len, uint32_t slots, uint32_t M);

void AdjustCiphertext(Scheme scheme, Ciphertext& ciphertext, double correction);

void SwitchModulus(uint64_t* res, uint64_t newModulus, uint64_t* input, uint64_t oldModulus, uint64_t N);

/**
* Sets all parameters for the linear method for the FFT-like method
*
* @param levelBudget - vector of budgets for the amount of levels in encoding
* and decoding
* @param dim1 - vector of inner dimension in the baby-step giant-step routine
* for encoding and decoding
* @param slots - number of slots to be bootstrapped
* @param correctionFactor - value to rescale message by to improve precision. If set to 0, we use the default logic. This value is only used when NATIVE_SIZE=64.
*/
void EvalBootstrapSetup(Scheme& scheme, uint32_t* levelBudget, uint32_t* dim1, uint32_t slots, uint32_t correctionFactor=0);


/**
* Virtual function to define the generation of all automorphism keys for EvalBT (with FFT evaluation).
* EvalBTKeyGen uses the baby-step/giant-step strategy.
*
* @param privateKey private key.
* @param slots - number of slots to be bootstrapped
* @return the dictionary of evaluation key indices.
*/
void EvalBootstrapKeyGen(Scheme& scheme, SecretKey secretKey, uint32_t slots);

/**
* Defines the bootstrapping evaluation of ciphertext
*
* The flavor of bootstrapping that uses the numIterations and precision parameters is described
* in the Meta-BTS paper.
* Source: Bae Y., Cheon J., Cho W., Kim J., and Kim T. META-BTS: Bootstrapping Precision
* Beyond the Limit. Cryptology ePrint Archive, Report
* 2022/1167. (https://eprint.iacr.org/2022/1167.pdf)
*
* @param ciphertext the input ciphertext.
* @param numIterations number of iterations to run iterative bootstrapping (Meta-BTS). Increasing the iterations increases the precision of bootstrapping.
* @param precision precision of initial bootstrapping algorithm. This value is
* determined by the user experimentally by first running EvalBootstrap with numIterations = 1 and precision = 0 (unused).
* @return the refreshed ciphertext.
*/
void EvalBootstrap(Scheme scheme, Ciphertext& result, Ciphertext ciphertext, uint32_t numIterations=1,uint32_t precision=0);

#endif //FRNS_BOOTSTRAP_MYBOOTSTRAPEXAMPLE_H
