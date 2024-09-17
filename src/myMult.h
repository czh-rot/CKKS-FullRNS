//
// Created by EYx on 2023/4/20.
//

#ifndef FRNS_BOOTSTRAP_MYMULT_H
#define FRNS_BOOTSTRAP_MYMULT_H
#include "Common.h"
#include "Context.h"
#include "Scheme.h"

//my Mult procedure
void myMult(Ciphertext& result,Ciphertext cipher1, Ciphertext cipher2,Scheme scheme);
void myMultCore(Ciphertext &cipher1, Ciphertext &cipher2,
                uint64_t* axbx1, uint64_t* axbx2, uint64_t* axax, uint64_t* bxbx, Scheme scheme);
void myMultAndEqual(Ciphertext& cipher1, Ciphertext cipher2, Scheme scheme);
void mySquare(Ciphertext& result, Ciphertext cipher, Scheme scheme);
void mySquareAndEqual(Ciphertext& cipher, Scheme scheme);

#endif //FRNS_BOOTSTRAP_MYMULT_H
