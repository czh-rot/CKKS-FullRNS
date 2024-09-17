//
// Created by EYx on 2023/5/27.
//

#ifndef FRNS_BOOTSTRAP_CKKSRNS_LEVELEDSHE_H
#define FRNS_BOOTSTRAP_CKKSRNS_LEVELEDSHE_H

#include "Scheme.h"
#include "Ciphertext.h"


typedef __int128 int128_t_;

//Note: 需要提前给factors分配ciphertext.limbs 大小的空间。
void GetElementForEvalMult(uint64_t *factors, Ciphertext ciphertext, double constant, Scheme scheme);

void EvalMultCoreInPlace(Ciphertext &ciphertext, double constant, Scheme scheme);

void EvalMultInPlace(Ciphertext &ciphertext, double constant, Scheme scheme);

//void EvalMult(Ciphertext &res, const Ciphertext input, double constant, Scheme scheme);
Ciphertext EvalMult( const Ciphertext input, double constant, Scheme scheme);

void DropLastElementAndScale(uint64_t*& a, long l, uint64_t* QlQlInvModqlDivqlModq, Scheme scheme);

void ModReduceInternalInPlace(Ciphertext& ciphertext, size_t levels, Scheme scheme);

#endif //FRNS_BOOTSTRAP_CKKSRNS_LEVELEDSHE_H
