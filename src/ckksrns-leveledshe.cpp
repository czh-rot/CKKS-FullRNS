//
// Created by EYx on 2023/5/27.
//

#include "ckksrns-leveledshe.h"

#include "Ciphertext.h"
#include "Scheme.h"

#include "myBootstrapExample.h"
#include "constants.h"

//std::vector<DCRTPoly::Integer> LeveledSHECKKSRNS::GetElementForEvalMult(ConstCiphertext<DCRTPoly> ciphertext,
//                                                                        double constant)
//Note: 需要提前给factors分配ciphertext.limbs 大小的空间。
void GetElementForEvalMult(uint64_t *factors, Ciphertext ciphertext, double constant, Scheme scheme) {
    uint32_t numTowers = ciphertext.curr_limbs;
    uint64_t p = scheme.context.p;
    uint64_t *qVec = scheme.context.qVec;
//    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

//    const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
//    uint32_t numTowers              = cv[0].GetNumOfElements();
//    std::vector<DCRTPoly::Integer> moduli(numTowers);
//    for (usint i = 0; i < numTowers; i++) {
//        moduli[i] = cv[0].GetElementAtIndex(i).GetModulus();
//    }

//    double scFactor = cryptoParams->GetScalingFactorReal(ciphertext->GetLevel());
    double scFactor = p;

//#if defined(HAVE_INT128)
    typedef int128_t_ DoubleInteger;
    int32_t MAX_BITS_IN_WORD = 126;
//#else
//    typedef int64_t DoubleInteger;
//    int32_t MAX_BITS_IN_WORD = LargeScalingFactorConstants::MAX_BITS_IN_WORD;
//#endif

    // Compute approxFactor, a value to scale down by, in case the value exceeds a 64-bit integer.
    int32_t logSF = static_cast<int32_t>(ceil(log2(fabs(scFactor))));
    int32_t logValid = (logSF <= MAX_BITS_IN_WORD) ? logSF : MAX_BITS_IN_WORD;
    int32_t logApprox = logSF - logValid;
    double approxFactor = pow(2, logApprox);

    DoubleInteger large = static_cast<DoubleInteger>(constant / approxFactor * scFactor + 0.5);
    DoubleInteger large_abs = (large < 0 ? -large : large);
    DoubleInteger bound = (uint64_t) 1 << 63;

//    std::vector<DCRTPoly::Integer> factors(numTowers);

    if (large_abs > bound) {
        for (uint32_t i = 0; i < numTowers; i++) {
            DoubleInteger reduced = large % qVec[i];

            factors[i] = (reduced < 0) ? static_cast<uint64_t>(reduced + qVec[i]) :
                         static_cast<uint64_t>(reduced);
        }
    } else {
        int64_t scConstant = static_cast<int64_t>(large);
        for (uint32_t i = 0; i < numTowers; i++) {
            int64_t reduced = scConstant % static_cast<int64_t>(qVec[i]);

            factors[i] = (reduced < 0) ? reduced + qVec[i] : reduced;
        }
    }

    // Scale back up by approxFactor within the CRT multiplications.
    if (logApprox > 0) { //Note: default to be 0 in FIXEDMANUAL mode
//        int32_t logStep = (logApprox <= LargeScalingFactorConstants::MAX_LOG_STEP) ?
//                          logApprox :
//                          LargeScalingFactorConstants::MAX_LOG_STEP;
//        DCRTPoly::Integer intStep = uint64_t(1) << logStep;
//        std::vector<DCRTPoly::Integer> crtApprox(numTowers, intStep);
//        logApprox -= logStep;
//
//        while (logApprox > 0) {
//            int32_t logStep = (logApprox <= LargeScalingFactorConstants::MAX_LOG_STEP) ?
//                              logApprox :
//                              LargeScalingFactorConstants::MAX_LOG_STEP;
//            DCRTPoly::Integer intStep = uint64_t(1) << logStep;
//            std::vector<DCRTPoly::Integer> crtSF(numTowers, intStep);
//            crtApprox = CKKSPackedEncoding::CRTMult(crtApprox, crtSF, moduli);
//            logApprox -= logStep;
//        }
//        factors = CKKSPackedEncoding::CRTMult(factors, crtApprox, moduli);
    }

//    return factors;
}

//void LeveledSHECKKSRNS::EvalMultCoreInPlace(Ciphertext<DCRTPoly>& ciphertext, double constant) const {
void EvalMultCoreInPlace(Ciphertext &ciphertext, double constant, Scheme scheme) {
//    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

//    std::vector<DCRTPoly::Integer> factors = GetElementForEvalMult(ciphertext, constant);
    uint64_t logN = scheme.context.logN;
    uint64_t curr_limbs = ciphertext.curr_limbs;
    uint64_t *factors = new uint64_t[curr_limbs];
    GetElementForEvalMult(factors, ciphertext, constant, scheme);

//    std::vector<DCRTPoly>& cv              = ciphertext->GetElements();
//    for (usint i = 0; i < cv.size(); ++i) {
//        cv[i] = cv[i] * factors;
//    }
    for (int i = 0; i < curr_limbs; ++i) {
        uint64_t *axj = ciphertext.ax + (i << logN);
        uint64_t *bxj = ciphertext.bx + (i << logN);
        scheme.context.qiMulConstAndEqual(axj, factors[i], i);
        scheme.context.qiMulConstAndEqual(bxj, factors[i], i);
    }
//    ciphertext->SetNoiseScaleDeg(ciphertext->GetNoiseScaleDeg() + 1);
//
//    double scFactor = cryptoParams->GetScalingFactorReal(ciphertext->GetLevel());
//    ciphertext->SetScalingFactor(ciphertext->GetScalingFactor() * scFactor);
}

//void LeveledSHECKKSRNS::EvalMultInPlace(Ciphertext<DCRTPoly>& ciphertext, double constant) const {
void EvalMultInPlace(Ciphertext &ciphertext, double constant, Scheme scheme) {
//    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

//    if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
//        if (ciphertext->GetNoiseScaleDeg() == 2) {
//            ModReduceInternalInPlace(ciphertext, BASE_NUM_LEVELS_TO_DROP);
//        }
//    }

    EvalMultCoreInPlace(ciphertext, constant, scheme);
}

Ciphertext EvalMult( const Ciphertext input, double constant, Scheme scheme)
{

    Ciphertext res=input;
    EvalMultCoreInPlace(res, constant, scheme);
    return res;
}

void DropLastElementAndScale(uint64_t*& a, long sizeQl, uint64_t* QlQlInvModqlDivqlModq, Scheme scheme)
{
    uint64_t logN=scheme.context.logN;
    uint64_t N=scheme.context.N;
    uint64_t* qVec=scheme.context.qVec;
    uint64_t* qrVec=scheme.context.qrVec;
    long* qTwok=scheme.context.qTwok;
    uint64_t** qInvModq=scheme.context.qInvModq;

    uint64_t *ra = new uint64_t[(sizeQl - 1) << logN]();
    uint64_t *lastPoly = a + ((sizeQl - 1) << logN); //copy the last poly
    uint64_t *extra = ra;

    //写法1：同openfhe的计算流程
    scheme.context.qiINTTAndEqual(lastPoly, sizeQl - 1); //lastPoly.SetFormat(Format::COEFFICIENT);

    for (long i = 0; i < sizeQl - 1; ++i) { //openfhe: l-1=extra.m_vectors.size()
        uint64_t *extrai = extra + (i << logN);
        copy(lastPoly, lastPoly + N, extrai); //auto temp = lastPoly;
        SwitchModulus(extrai,qVec[i],lastPoly,qVec[sizeQl-1],N);//temp.SwitchModulus()
        for (long n = 0; n < N; ++n) {
            mulModBarrett(extrai[n],extrai[n],QlQlInvModqlDivqlModq[i],
                          qVec[i], qrVec[i], qTwok[i]);
        }
    }

    scheme.context.NTTAndEqual(extra, sizeQl - 1);

    for (int i = 0; i < sizeQl - 1; ++i) {
        uint64_t *rai = ra + (i << logN);
        uint64_t *extrai = extra + (i << logN);
        uint64_t *ai = a + (i << logN);
        for (long n = 0; n < N; ++n) {
            mulModBarrett(ai[n], ai[n], qInvModq[sizeQl - 1][i], qVec[i], qrVec[i], qTwok[i]);
            addMod(rai[n], ai[n], extrai[n], qVec[i]);
        }
    }
    delete[] a;
    a = ra;

//    //写法2：用openfhe的算法，在FRNS原来的rescale上改的
//    scheme.context.qiINTTAndEqual(lastPoly, sizeQl - 1); //lastPoly.SetFormat(Format::COEFFICIENT);
//
//    for (long i = 0; i < sizeQl - 1; ++i) { //openfhe: l-1=extra.m_vectors.size()
//        uint64_t *rai = ra + (i << logN);
//        uint64_t *extrai = extra + (i << logN);
//        uint64_t *ai = a + (i << logN);
//
//        copy(lastPoly, lastPoly + N, extrai); //auto temp = lastPoly;
//        SwitchModulus(extrai,qVec[i],lastPoly,qVec[sizeQl-1],N);//openfhe: temp.SwitchModulus()
//        for (long n = 0; n < N; ++n) {
//            mulModBarrett(extrai[n],extrai[n],QlQlInvModqlDivqlModq[i],
//                          qVec[i], qrVec[i], qTwok[i]);
//        }
//        scheme.context.qiNTTAndEqual(extrai, i);
//
//        for (long n = 0; n < N; ++n) {
//            mulModBarrett(ai[n], ai[n], qInvModq[sizeQl - 1][i], qVec[i], qrVec[i], qTwok[i]);
//            addMod(rai[n], ai[n], extrai[n], qVec[i]);
//        }
//    }
//    delete[] a;
//    a = ra;
}

//void LeveledSHECKKSRNS::ModReduceInternalInPlace(Ciphertext<DCRTPoly>& ciphertext, size_t levels) const {
void ModReduceInternalInPlace(Ciphertext& ciphertext, size_t levels, Scheme scheme)
{
    size_t sizeQ  = scheme.context.L;
    size_t sizeQl = ciphertext.curr_limbs;
    size_t diffQl = sizeQ - sizeQl;

    uint64_t **QlQlInvModqlDivqlModq=scheme.context.QlQlInvModqlDivqlModq;

    for (size_t l = 0; l < levels; ++l) {
        uint64_t curr_limbs=ciphertext.curr_limbs;
        DropLastElementAndScale(ciphertext.bx, curr_limbs, QlQlInvModqlDivqlModq[diffQl+l], scheme);
        DropLastElementAndScale(ciphertext.ax, curr_limbs, QlQlInvModqlDivqlModq[diffQl+l], scheme);
        ciphertext.curr_limbs-=1;
    }
}

////std::vector<DCRTPoly::Integer> LeveledSHECKKSRNS::GetElementForEvalAddOrSub(ConstCiphertext<DCRTPoly> ciphertext,double constant)
////Note: 需要提前给 crtConstant 分配ciphertext.limbs 大小的空间。
//void GetElementForEvalAddOrSub(uint64_t *crtConstant, Ciphertext ciphertext, double constant, Scheme scheme)
//{
////    const std::vector<DCRTPoly>& cv = ciphertext->GetElements();
//    uint32_t L=scheme.context.L; //usint sizeQl = cv[0].GetNumOfElements();
//    uint64_t p=scheme.context.p; //usint sizeQl = cv[0].GetNumOfElements();
//    uint32_t sizeQl=ciphertext.curr_limbs; //usint sizeQl = cv[0].GetNumOfElements();
//    uint64_t *moduli = scheme.context.qVec;
//    ScalingTechnique rescaleTech=scheme.rescaleTech;
////    std::vector<DCRTPoly::Integer> moduli(sizeQl);
////    for (usint i = 0; i < sizeQl; i++) {
////        moduli[i] = cv[0].GetElementAtIndex(i).GetModulus();
////    }
//
////    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());
//
//    double scFactor = 0;
////    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT && ciphertext->GetLevel() == 0)
//    if (rescaleTech == FLEXIBLEAUTOEXT && ciphertext.curr_limbs==L) {
////        scFactor = cryptoParams->GetScalingFactorRealBig(ciphertext->GetLevel());
//    }
//    else {
//        if (rescaleTech == FIXEDMANUAL){
//            scFactor = p;
//        }
//        else{
////        scFactor = cryptoParams->GetScalingFactorReal(ciphertext->GetLevel());
//        }
//    }
//
//    // Compute approxFactor, a value to scale down by, in case the value exceeds a 64-bit integer.
//    int32_t logSF    = static_cast<int32_t>(ceil(log2(fabs(scFactor))));
//    int32_t logValid = (logSF <= LargeScalingFactorConstants::MAX_BITS_IN_WORD) ?
//                       logSF :
//                       LargeScalingFactorConstants::MAX_BITS_IN_WORD;
//    int32_t logApprox   = logSF - logValid;
//    double approxFactor = pow(2, logApprox);
//
//    uint64_t scConstant = static_cast<uint64_t>(constant * scFactor / approxFactor + 0.5);
////    std::vector<DCRTPoly::Integer> crtConstant(sizeQl, scConstant);
//    for (int i = 0; i < sizeQl; ++i) {
//        crtConstant[i]=scConstant;
//    }
//
//    // Scale back up by approxFactor within the CRT multiplications.
//    if (logApprox > 0) {
////        int32_t logStep = (logApprox <= LargeScalingFactorConstants::MAX_LOG_STEP) ?
////                          logApprox :
////                          LargeScalingFactorConstants::MAX_LOG_STEP;
////        DCRTPoly::Integer intStep = uint64_t(1) << logStep;
////        std::vector<DCRTPoly::Integer> crtApprox(sizeQl, intStep);
////        logApprox -= logStep;
////
////        while (logApprox > 0) {
////            int32_t logStep = (logApprox <= LargeScalingFactorConstants::MAX_LOG_STEP) ?
////                              logApprox :
////                              LargeScalingFactorConstants::MAX_LOG_STEP;
////            DCRTPoly::Integer intStep = uint64_t(1) << logStep;
////            std::vector<DCRTPoly::Integer> crtSF(sizeQl, intStep);
////            crtApprox = CKKSPackedEncoding::CRTMult(crtApprox, crtSF, moduli);
////            logApprox -= logStep;
////        }
////        crtConstant = CKKSPackedEncoding::CRTMult(crtConstant, crtApprox, moduli);
//    }
//
//    // In FLEXIBLEAUTOEXT mode at level 0, we don't use the depth to calculate the scaling factor,
//    // so we return the value before taking the depth into account.
////    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT && ciphertext->GetLevel() == 0) {
//    if (rescaleTech == FLEXIBLEAUTOEXT && ciphertext.curr_limbs==L) {
////        return crtConstant;
//    }
//
//    uint64_t intScFactor = static_cast<uint64_t>(scFactor + 0.5);
//    std::vector<DCRTPoly::Integer> crtScFactor(sizeQl, intScFactor);
//
//    for (usint i = 1; i < ciphertext->GetNoiseScaleDeg(); i++) {
//        //fixme: 这个函数的意义就在这个循环，但是需要给密文引入NoiseScaleDeg这个参数。然后再看下面这个crtmult能不能搞成mulmodBarrett之类的
//        crtConstant = CKKSPackedEncoding::CRTMult(crtConstant, crtScFactor, moduli);
//    }
//
//    return crtConstant;
//}