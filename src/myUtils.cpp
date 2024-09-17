//
// Created by EYx on 2023/4/29.
//

#include "Ciphertext.h"

#include "ckksrns-fhe.h"
#include "ckksrns-leveledshe.h"

#include "myUtils.h"
#include "myMult.h"

void SetZero(uint64_t * a, long Start, long End)
{
    for (int i = Start; i < End; ++i) {
        a[i]=0;
    }
}

void CopyValue(uint64_t *res, uint64_t *input, long Start, long End) {
    for (int i = Start; i < End; ++i) {
        res[i]=input[i];
    }
}

double ConvertToDouble(uint64_t input) {
    return static_cast<double>(input);
}

void ApplyDoubleAngleIterations(Ciphertext& ciphertext,Scheme scheme) {
//    auto cc = ciphertext->GetCryptoContext();

    int32_t r = R;
    for (int32_t j = 1; j < r + 1; j++) {
        mySquareAndEqual(ciphertext,scheme);//cc->EvalSquareInPlace(ciphertext);
        scheme.addAndEqual(ciphertext,ciphertext);//ciphertext    = cc->EvalAdd(ciphertext, ciphertext);

        //fixme: openfhe 中这个 reScale 是放在循环体的最后做的，但是这里如果和现在的addconstandequal 交换就会算不对
        //因为openfhe的addconst和比这里实现的高级
//        scheme.reScaleByAndEqual(ciphertext,1);//cc->ModReduceInPlace(ciphertext);
        ModReduceInternalInPlace(ciphertext,1,scheme);

        double scalar = -1.0 / std::pow((2.0 * M_PI), std::pow(2.0, j - r));
        scheme.addConstAndEqual(ciphertext,scalar);//cc->EvalAddInPlace(ciphertext, scalar);

    }
}