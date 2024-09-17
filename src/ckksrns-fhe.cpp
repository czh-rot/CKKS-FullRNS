//
// Created by EYx on 2023/5/4.
//

#include "ckksrns-fhe.h"

#include "Scheme.h"

//------------------------------------------------------------------------------
// Precomputations for CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------

std::vector<Plaintext> EvalLinearTransformPrecompute(
        Scheme scheme, const std::vector<std::vector<std::complex<double>>>& A,
        double scale, uint32_t L)
{
    if (A[0].size() != A.size()) {
        cout<< "The matrix passed to EvalLTPrecompute is not square"<<endl;
    }

    uint32_t slots = A.size();
    CKKSBootstrapPrecom precom=scheme.precom;
    uint32_t M=scheme.context.M;
//    auto pair = m_bootPrecomMap.find(slots);
//    if (pair == m_bootPrecomMap.end()) {
//    std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
//                         std::string(" slots were not generated") +
//                         std::string(" Need to call EvalBootstrapSetup to proceed"));
//    OPENFHE_THROW(type_error, errorMsg);
//    }
//    const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;
//    uint32_t M = cc.GetCyclotomicOrder();

    // Computing the baby-step bStep and the giant-step gStep.
    int bStep = (precom.m_dim1 == 0) ? ceil(sqrt(slots)) : precom.m_dim1;
    int gStep = ceil(static_cast<double>(slots) / bStep);

    //fixme: // make sure the plaintext is created only with the necessary amount of moduli

//    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

//    ILDCRTParams<DCRTPoly::Integer> elementParams = *(cryptoParams->GetElementParams());


//
//    uint32_t towersToDrop = 0;
//    if (L != 0) {
//        towersToDrop = elementParams.GetParams().size() - L - 1;
//    }
//
//    for (uint32_t i = 0; i < towersToDrop; i++) {
//        elementParams.PopLastParam();
//    }
//
//    auto paramsQ = elementParams.GetParams();
//    usint sizeQ  = paramsQ.size();
//    auto paramsP = cryptoParams->GetParamsP()->GetParams();
//    usint sizeP  = paramsP.size();
//
//    std::vector<NativeInteger> moduli(sizeQ + sizeP);
//    std::vector<NativeInteger> roots(sizeQ + sizeP);
//
//    for (size_t i = 0; i < sizeQ; i++) {
//        moduli[i] = paramsQ[i]->GetModulus();
//        roots[i]  = paramsQ[i]->GetRootOfUnity();
//    }
//
//    for (size_t i = 0; i < sizeP; i++) {
//        moduli[sizeQ + i] = paramsP[i]->GetModulus();
//        roots[sizeQ + i]  = paramsP[i]->GetRootOfUnity();
//    }
//
//    auto elementParamsPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);
//    //  auto elementParamsPtr2 = std::dynamic_pointer_cast<typename DCRTPoly::Params>(elementParamsPtr);
//
    std::vector<Plaintext> result(slots);
//// parallelizing the loop (below) with OMP causes a segfault on MinGW
//// see https://github.com/openfheorg/openfhe-development/issues/176
//#if !defined(__MINGW32__) && !defined(__MINGW64__)
////    #pragma omp parallel for
//#endif
//    for (int j = 0; j < gStep; j++) {
//        int offset = -bStep * j;
//        for (int i = 0; i < bStep; i++) {
//            if (bStep * j + i < static_cast<int>(slots)) {
//            auto diag = ExtractShiftedDiagonal(A, bStep * j + i);
//            for (uint32_t k = 0; k < diag.size(); k++)
//                diag[k] *= scale;
//
//                result[bStep * j + i] =
//                    MakeAuxPlaintext(cc, elementParamsPtr, Rotate(diag, offset), 1, towersToDrop, diag.size());
//            }
//        }
//    }
    return result;
}