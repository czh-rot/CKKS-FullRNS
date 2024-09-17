//
// Created by EYx on 2023/5/4.
//

#include "Ciphertext.h"
#include "Scheme.h"
#include "myLT.h"
#include "myRotation.h"
#include "myKeySwitch.h"
#include "myUtils.h"
#include "ckksrns-utils.h"
#include "ckksrns-leveledshe.h"

//------------------------------------------------------------------------------
// EVALUATION: CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------


void EvalLinearTransform(Ciphertext &result,
                         Plaintext *A, uint32_t A_len,
                         const Ciphertext &ct,
                         Scheme scheme) {
    uint32_t slots = A_len;

    CKKSBootstrapPrecom precom = scheme.precom;
    // Computing the baby-step bStep and the giant-step gStep.
    uint32_t bStep = (precom.m_dim1 == 0) ? ceil(sqrt(slots)) : precom.m_dim1;
    uint32_t gStep = ceil(static_cast<double>(slots) / bStep);

    uint32_t logN = scheme.context.logN;
    uint32_t N = scheme.context.N;
    uint32_t M = scheme.context.M;
    uint32_t K = scheme.context.K;

    //Note: 本函数不改变密文的模数个数，所以统一的在这里声明一些参数
    long curr_limbs = ct.curr_limbs;
    long special_limbs=K;
    long limbsExt = curr_limbs + special_limbs;
    long lenExt = (limbsExt<<logN);
    long len = (curr_limbs<<logN);

    uint32_t alpha = scheme.context.alpha;
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    // computes the NTTs for each CRT limb (for the hoisted automorphisms used
    // later on)
    //auto digits = cc->EvalFastRotationPrecompute(ct);
    long digits_len = beta * (limbsExt << logN);
    uint64_t *digits = new uint64_t[digits_len];//allocate at where it calls
    EvalFastRotationPrecompute(digits, ct, scheme);

    Ciphertext *fastRotationExt = new Ciphertext[bStep - 1];

    // hoisted automorphisms
    //#pragma omp parallel for
    for (uint32_t j = 1; j < bStep; j++) {
        fastRotationExt[j - 1].N=ct.N;
        fastRotationExt[j - 1].slots=ct.slots;
        fastRotationExt[j - 1].curr_limbs=ct.curr_limbs;
        fastRotationExt[j - 1].special_limbs=special_limbs;
        fastRotationExt[j - 1].ax = new uint64_t[lenExt]();
        fastRotationExt[j - 1].bx = new uint64_t[lenExt]();

        //fastRotation[j - 1] = cc->EvalFastRotationExt(ct, j, digits, true);
        FastRotate_KeyGen(scheme.secretKey, j, scheme);//todo: 将密钥生成挪到预计算中去
        EvalFastRotationExt(fastRotationExt[j - 1], ct, j, digits, true, scheme);
    }
    delete[]digits;

    uint64_t *resultExt_ax = new uint64_t[lenExt]();
    uint64_t *resultExt_bx = new uint64_t[lenExt]();
    Ciphertext resultExt(resultExt_ax, resultExt_bx, ct.N, ct.slots,
                         curr_limbs,special_limbs);

    uint64_t *first = new uint64_t[len];

    for (uint32_t j = 0; j < gStep; j++) {
        //Ciphertext<DCRTPoly> inner = EvalMultExt(cc->KeySwitchExt(ct, true), A[bStep * j]);
        // <=>
        // 1) Ciphertext<DCRTPoly> ctExt = cc->KeySwitchExt(ct, true),
        // 2) Ciphertext<DCRTPoly> inner = EvalMultExt(ctExt, A[bStep * j]);

        uint64_t *inner_ax = new uint64_t[lenExt]();
        uint64_t *inner_bx = new uint64_t[lenExt]();
        Ciphertext innerExt(inner_ax, inner_bx, ct.N, ct.slots, curr_limbs, special_limbs);

        //todo: innerExt 的Extent部分可能是零值？这里是否可以化简成不Ext的Mult
        KeySwitchExt(innerExt, ct, NORMAL_CIPHER_SIZE, true, scheme);
        scheme.EvalMultExtInPlace(innerExt, A[bStep * j]);

        for (uint32_t i = 1; i < bStep; i++) {
            if (bStep * j + i < slots) {
                //EvalAddExtInPlace(inner, EvalMultExt(fastRotation[i - 1], A[bStep * j + i]));
                // <=>
                // 1) Ciphertext<DCRTPoly> tmp = EvalMultExt(fastRotation[i - 1], A[bStep * j + i])
                // 2) EvalAddExtInPlace(inner, tmp)
                uint64_t *tmp_ax = new uint64_t[lenExt]();
                uint64_t *tmp_bx = new uint64_t[lenExt]();
                Ciphertext tmpExt(tmp_ax, tmp_bx, ct.N, ct.slots, curr_limbs, special_limbs);
                scheme.EvalMultExt(tmpExt, fastRotationExt[i - 1], A[bStep * j + i]);//todo: check correctness
                scheme.EvalAddExtInPlace(innerExt, tmpExt);
            }
        }

        if (j == 0) {
            //first = cc->KeySwitchDownFirstElement(inner);
            KeySwitchDownFirstElement(first, innerExt.bx, curr_limbs, scheme);
            //auto elements = inner->GetElements();
            //elements[0].SetValuesToZero();
            //inner->SetElements(elements);
            SetZero(innerExt.bx, 0, lenExt);

            resultExt=innerExt;
        }
        else {
            uint64_t *innerKSDown_ax = new uint64_t[len]();
            uint64_t *innerKSDown_bx = new uint64_t[len]();
            Ciphertext innerKSDown(innerKSDown_ax, innerKSDown_bx, ct.N, ct.slots, curr_limbs);

            // inner = cc->KeySwitchDown(inner);
            // => innerKSDown = cc->KeySwitchDown(inner);
            KeySwitchDown(innerKSDown.ax, innerKSDown.bx, innerExt.ax, innerExt.bx, curr_limbs, scheme);

            // Find the automorphism index that corresponds to rotation index.
            long autoIndex = FindAutomorphismIndex2nComplex(bStep * j,M);
            //std::vector<usint> map(N);
            uint32_t map_len = N;
            uint32_t *map = new uint32_t[map_len]();
            PrecomputeAutoMap(N, autoIndex, map);//PrecomputeAutoMap(N, autoIndex, &map);

            //DCRTPoly firstCurrent = innerKSDown->GetElements()[0].AutomorphismTransform(autoIndex, map);
            uint64_t *firstCurrent = new uint64_t[curr_limbs << logN];
            AutomorphismTransform(firstCurrent, innerKSDown.bx, curr_limbs, N, autoIndex, map);
            scheme.context.addAndEqual(first, firstCurrent, curr_limbs);//first += firstCurrent;

            //auto innerDigits = cc->EvalFastRotationPrecompute(innerKSDown);
            uint64_t *innerDigits = new uint64_t[digits_len];//allocate at where it calls
            EvalFastRotationPrecompute(innerDigits, innerKSDown,scheme);

            uint64_t *innerKSDownExt_ax = new uint64_t[lenExt]();
            uint64_t *innerKSDownExt_bx = new uint64_t[lenExt]();
            Ciphertext innerKSDownExt(innerKSDownExt_ax, innerKSDownExt_bx, ct.N, ct.slots,
                                      curr_limbs, special_limbs);
            FastRotate_KeyGen(scheme.secretKey, bStep * j, scheme);//todo: 将密钥生成挪到预计算中去

            //Ciphertext<DCRTPoly> inner_Ext=cc->EvalFastRotationExt(inner, bStep * j, innerDigits, false);
            EvalFastRotationExt(innerKSDownExt, innerKSDown, bStep * j, innerDigits, false, scheme);
            scheme.EvalAddExtInPlace(resultExt, innerKSDownExt);//EvalAddExtInPlace(result, inner_Ext);
        }
    }

    //result = cc->KeySwitchDown(resultExt); <= result = cc->KeySwitchDown(result);
    KeySwitchDown(result.ax, result.bx, resultExt.ax, resultExt.bx, curr_limbs, scheme);
    scheme.context.addAndEqual(result.bx, first, curr_limbs);//    elements[0] += first;

    delete[] first;
//    return result;
}


//DongYilin add
void EvalSlotsToCoeffs(Ciphertext &result,
                       Plaintext** A,uint32_t A_len,
                       Ciphertext ctxt,
                       Scheme scheme)
{
    uint32_t slots = A_len;//uint32_t slots = ctxt->GetSlots();

    // auto pair = m_bootPrecomMap.find(slots);
    // if (pair == m_bootPrecomMap.end()) {
    //     std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
    //                          std::string(" slots were not generated") +
    //                          std::string(" Need to call EvalBootstrapSetup and EvalBootstrapKeyGen to proceed"));
    //     OPENFHE_THROW(type_error, errorMsg);
    // }

    // const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;

    // auto cc = ctxt->GetCryptoContext();

    uint32_t N = scheme.context.N;//uint32_t N = cc->GetRingDimension();
    uint32_t M = scheme.context.M;//uint32_t M = cc->GetCyclotomicOrder();
    uint32_t K = scheme.context.K;
    uint32_t logN = scheme.context.logN;
    uint32_t special_limbs=K;

    CKKSBootstrapPrecom& precom=scheme.precom;
    int32_t levelBudget     = precom.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];
    int32_t layersCollapse  = precom.m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = precom.m_paramsDec[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = precom.m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = precom.m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = precom.m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = precom.m_paramsDec[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = precom.m_paramsDec[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = precom.m_paramsDec[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    //auto algo = cc->GetScheme();

    int32_t flagRem = 0;

    if (remCollapse != 0) {
        flagRem = 1;
    }

    // precompute the inner and outer rotations

    uint32_t** rot_in=new uint32_t*[levelBudget];//std::vector<std::vector<int32_t>> rot_in(levelBudget);

    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        if (flagRem == 1 && i == uint32_t(levelBudget - 1)) {
            // remainder corresponds to index 0 in encoding and to last index in decoding
            rot_in[i]=new uint32_t[numRotationsRem+1];//rot_in[i] = std::vector<int32_t>(numRotationsRem + 1);
        }
        else {
            rot_in[i]=new uint32_t[numRotations+1];//rot_in[i] = std::vector<int32_t>(numRotations + 1);
        }
    }

    uint32_t** rot_out=new uint32_t*[levelBudget];//std::vector<std::vector<int32_t>> rot_out(levelBudget);

    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        rot_out[i]=new uint32_t[b+bRem];//rot_out[i] = std::vector<int32_t>(b + bRem);
    }

    for (int32_t s = 0; s < levelBudget - flagRem; s++) {
        for (int32_t j = 0; j < g; j++) {
            rot_in[s][j] =
                    ReduceRotation((j - int32_t((numRotations + 1) / 2) + 1) * (1 << (s * layersCollapse)), M / 4);
        }

        for (int32_t i = 0; i < b; i++) {
            rot_out[s][i] = ReduceRotation((g * i) * (1 << (s * layersCollapse)), M / 4);
        }
    }

    if (flagRem) {
        int32_t s = levelBudget - flagRem;
        for (int32_t j = 0; j < gRem; j++) {
            rot_in[s][j] =
                    ReduceRotation((j - int32_t((numRotationsRem + 1) / 2) + 1) * (1 << (s * layersCollapse)), M / 4);
        }

        for (int32_t i = 0; i < bRem; i++) {
            rot_out[s][i] = ReduceRotation((gRem * i) * (1 << (s * layersCollapse)), M / 4);
        }
    }

    //  No need for Encrypted Bit Reverse
    //Ciphertext<DCRTPoly> result = ctxt->Clone();
    result=ctxt;//fixme: 是否需要改掉C2S的 tmp_result

    // hoisted automorphisms
    for (int32_t s = 0; s < levelBudget - flagRem; s++) {
        if (s != 0) {
            ModReduceInternalInPlace(result,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP);
        }

        long curr_limbs = result.curr_limbs;
        long limbsExt = curr_limbs + special_limbs;
        long lenExt = (limbsExt<<logN);
        long len = (curr_limbs<<logN);
        uint32_t alpha = scheme.context.alpha;
        long beta = std::ceil((curr_limbs * 1.0 / alpha));

        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        //auto digits = cc->EvalFastRotationPrecompute(result);
        
        long digits_len = beta * lenExt;
        uint64_t* digits = new uint64_t[digits_len];
        EvalFastRotationPrecompute(digits,result,scheme);


        Ciphertext* fastRotationExt=new Ciphertext[g];//std::vector<Ciphertext<DCRTPoly>> fastRotation(g);
//#pragma omp parallel for
        for (int32_t j = 0; j < g; j++) {
            fastRotationExt[j].N=N;
            fastRotationExt[j].slots=slots;
            fastRotationExt[j].curr_limbs=curr_limbs;
            fastRotationExt[j].special_limbs=special_limbs;
            fastRotationExt[j].ax = new uint64_t[lenExt]();
            fastRotationExt[j].bx = new uint64_t[lenExt]();
            if (rot_in[s][j] != 0) {
                //fastRotation[j] = cc->EvalFastRotationExt(result, rot_in[s][j], digits, true);

                FastRotate_KeyGen(scheme.secretKey,rot_in[s][j],scheme);//先生成密钥
                EvalFastRotationExt(fastRotationExt[j], result, rot_in[s][j], digits, true, scheme);
            }
            else {
                KeySwitchExt(fastRotationExt[j],result, NORMAL_CIPHER_SIZE,true,scheme);//fastRotation[j] = cc->KeySwitchExt(result, true);
            }
        }
        delete[]digits;

        //Ciphertext<DCRTPoly> outer;
        uint64_t* outer_ax = new uint64_t[lenExt]();
        uint64_t* outer_bx = new uint64_t[lenExt]();
        Ciphertext outer(outer_ax,outer_bx,N,slots,curr_limbs,special_limbs);

        uint64_t* first= new uint64_t [len];//DCRTPoly first;

        for (int32_t i = 0; i < b; i++) {
            // for the first iteration with j=0:
            int32_t G = g * i;
            //Ciphertext<DCRTPoly> inner= EvalMultExt(fastRotation[0], A[s][G]);

            uint64_t *inner_ax = new uint64_t[lenExt]();
            uint64_t *inner_bx = new uint64_t[lenExt]();
            Ciphertext inner(inner_ax, inner_bx, ctxt.N, ctxt.slots, curr_limbs, special_limbs);
            scheme.EvalMultExt(inner, fastRotationExt[0], A[s][G]);

            // continue the loop
            for (int32_t j = 1; j < g; j++) {
                if ((G + j) != int32_t(numRotations)) {
                    //EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[s][G + j]));
                    uint64_t *tmp_ax = new uint64_t[lenExt]();
                    uint64_t *tmp_bx = new uint64_t[lenExt]();
                    Ciphertext tmpExt(tmp_ax, tmp_bx, ctxt.N, ctxt.slots, curr_limbs, special_limbs);
                    scheme.EvalMultExt(tmpExt, fastRotationExt[j], A[s][G + j]);//todo: check correctness
                    scheme.EvalAddExtInPlace(inner, tmpExt);
                }
            }

            if (i == 0) {
                //first         = cc->KeySwitchDownFirstElement(inner);
                KeySwitchDownFirstElement(first,inner.bx,curr_limbs,scheme);

                // auto elements = inner->GetElements();
                // elements[0].SetValuesToZero();
                // inner->SetElements(elements);
                SetZero(inner.bx,0,lenExt);

                outer = inner;
                // CopyValue(outer.ax,inner.ax,0,inner_len);
                // CopyValue(outer.bx,inner.bx,0,inner_len);
            }
            else {
                if (rot_out[s][i] != 0) {
                    //inner = cc->KeySwitchDown(inner);
                    uint64_t* innerKSDown_ax = new uint64_t[len]();
                    uint64_t* innerKSDown_bx = new uint64_t[len]();
                    Ciphertext innerKSDown(innerKSDown_ax, innerKSDown_bx, N, slots, curr_limbs);
                    KeySwitchDown(innerKSDown.ax, innerKSDown.bx, inner.ax, inner.bx, curr_limbs, scheme);

                    // Find the automorphism index that corresponds to rotation index.
                    long autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);
                    uint32_t map_len=N;
                    uint32_t * map  = new uint32_t [map_len]();
                    PrecomputeAutoMap(N, autoIndex, map);

                    //first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    uint64_t * firstCurrent= new uint64_t[len];
                    AutomorphismTransform(firstCurrent, innerKSDown.bx, curr_limbs, N, autoIndex, map);
                    scheme.context.addAndEqual(first,firstCurrent,curr_limbs);//first += firstCurrent;

                    //auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    uint64_t* innerDigits = new uint64_t[digits_len];
                    EvalFastRotationPrecompute(innerDigits, innerKSDown, scheme);

                    //EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[s][i], innerDigits, false));
                    uint64_t *innerKSDownExt_ax = new uint64_t[lenExt]();
                    uint64_t *innerKSDownExt_bx = new uint64_t[lenExt]();
                    Ciphertext innerKSDownExt(innerKSDownExt_ax, innerKSDownExt_bx, ctxt.N, ctxt.slots,curr_limbs, special_limbs);
                    FastRotate_KeyGen(scheme.secretKey, rot_out[s][i], scheme);//todo: 将密钥生成挪到预计算中去

                    EvalFastRotationExt(innerKSDownExt, innerKSDown, rot_out[s][i], innerDigits, false, scheme);
                    scheme.EvalAddExtInPlace(outer, innerKSDownExt);
                }
                else {
                    //first += cc->KeySwitchDownFirstElement(inner);
                    uint64_t * tmp_first= new uint64_t[len];
                    KeySwitchDownFirstElement(tmp_first,inner.bx,curr_limbs,scheme);
                    scheme.context.addAndEqual(first,tmp_first,curr_limbs);//first += tmp_first;

                    // auto elements = inner->GetElements();
                    // elements[0].SetValuesToZero();
                    // inner->SetElements(elements);
                    SetZero(inner.bx,0,lenExt);

                    scheme.EvalAddExtInPlace(outer, inner);//EvalAddExtInPlace(outer, inner);
                }
            }
        }

        //result                          = cc->KeySwitchDown(outer);
        KeySwitchDown(result.ax, result.bx, outer.ax, outer.bx, curr_limbs, scheme);

        // std::vector<DCRTPoly>& elements = result->GetElements();
        // elements[0] += first;
        scheme.context.addAndEqual(result.bx,first,curr_limbs);
    }
    
    if (flagRem) {
        ModReduceInternalInPlace(result,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP);
        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        //auto digits = cc->EvalFastRotationPrecompute(result);

        long curr_limbs = result.curr_limbs;
        long limbsExt = curr_limbs + special_limbs;
        long lenExt = (limbsExt<<logN);
        long len = (curr_limbs<<logN);
        uint32_t alpha = scheme.context.alpha;
        long beta = std::ceil((curr_limbs * 1.0 / alpha));
        
        long digits_len = beta * lenExt;
        uint64_t* digits = new uint64_t[digits_len];
        EvalFastRotationPrecompute(digits,result,scheme);

        Ciphertext* fastRotationExt=new Ciphertext[gRem];//std::vector<Ciphertext<DCRTPoly>> fastRotation(gRem);

        int32_t s = levelBudget - flagRem;
//#pragma omp parallel for
        for (int32_t j = 0; j < gRem; j++) {
            fastRotationExt[j].N=N;
            fastRotationExt[j].slots=slots;
            fastRotationExt[j].curr_limbs=curr_limbs;
            fastRotationExt[j].special_limbs=special_limbs;
            fastRotationExt[j].ax = new uint64_t[lenExt]();
            fastRotationExt[j].bx = new uint64_t[lenExt]();
            if (rot_in[s][j] != 0) {
                //fastRotation[j] = cc->EvalFastRotationExt(result, rot_in[s][j], digits, true);
                FastRotate_KeyGen(scheme.secretKey,rot_in[s][j],scheme);//先生成密钥
                EvalFastRotationExt(fastRotationExt[j], result, rot_in[s][j], digits, true, scheme);
            }
            else {
                KeySwitchExt(fastRotationExt[j],result,NORMAL_CIPHER_SIZE,true,scheme);//fastRotation[j] = cc->KeySwitchExt(result, true);
            }
        }
        delete[]digits;


        Ciphertext outer;//Ciphertext<DCRTPoly> outer;
        uint64_t* first= new uint64_t [len];//DCRTPoly first;

        for (int32_t i = 0; i < bRem; i++) {
            //Ciphertext<DCRTPoly> inner;
            // for the first iteration with j=0:
            int32_t GRem = gRem * i;

            //inner        = EvalMultExt(fastRotation[0], A[s][GRem]);
            uint64_t *inner_ax = new uint64_t[lenExt]();
            uint64_t *inner_bx = new uint64_t[lenExt]();
            Ciphertext inner(inner_ax, inner_bx, N, slots, curr_limbs, special_limbs);
            scheme.EvalMultExt(inner, fastRotationExt[0], A[s][GRem]);

            // continue the loop
            for (int32_t j = 1; j < gRem; j++) {
                if ((GRem + j) != int32_t(numRotationsRem)){
                    //EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[s][GRem + j]));
                    uint64_t *tmp_ax = new uint64_t[lenExt]();
                    uint64_t *tmp_bx = new uint64_t[lenExt]();
                    Ciphertext tmpExt(tmp_ax, tmp_bx, ctxt.N, ctxt.slots, curr_limbs, special_limbs);
                    scheme.EvalMultExt(tmpExt, fastRotationExt[j], A[s][GRem + j]);//todo: check correctness
                    scheme.EvalAddExtInPlace(inner, tmpExt);
                }
            }

            if (i == 0) {
                //first         = cc->KeySwitchDownFirstElement(inner);
                KeySwitchDownFirstElement(first,inner.bx,curr_limbs,scheme);

                // auto elements = inner->GetElements();
                // elements[0].SetValuesToZero();
                // inner->SetElements(elements);
                SetZero(inner.bx,0,lenExt);

                outer = inner;
                // CopyValue(outer.ax,inner.ax,0,inner_len);
                // CopyValue(outer.bx,inner.bx,0,inner_len);
            }
            else {
                if (rot_out[s][i] != 0) {
                    //inner = cc->KeySwitchDown(inner);
                    uint64_t* innerKSDown_ax = new uint64_t[len]();
                    uint64_t* innerKSDown_bx = new uint64_t[len]();
                    Ciphertext innerKSDown(innerKSDown_ax, innerKSDown_bx, N, slots, curr_limbs);
                    KeySwitchDown(innerKSDown.ax, innerKSDown.bx, inner.ax, inner.bx, curr_limbs, scheme);

                    // Find the automorphism index that corresponds to rotation index.
                    long autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);//usint autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);

                    //std::vector<usint> map(N);
                    uint32_t map_len=N;
                    uint32_t * map  = new uint32_t [map_len]();
                    PrecomputeAutoMap(N, autoIndex, map);//PrecomputeAutoMap(N, autoIndex, &map);

                    //first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    uint64_t * firstCurrent= new uint64_t[len];
                    AutomorphismTransform(firstCurrent, innerKSDown.bx, curr_limbs, N, autoIndex, map);
                    scheme.context.addAndEqual(first,firstCurrent,curr_limbs);//first += firstCurrent;

                    //auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    uint64_t* innerDigits = new uint64_t[digits_len];
                    EvalFastRotationPrecompute(innerDigits, innerKSDown, scheme);

                    //EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[s][i], innerDigits, false));
                    uint64_t *innerKSDownExt_ax = new uint64_t[lenExt]();
                    uint64_t *innerKSDownExt_bx = new uint64_t[lenExt]();
                    Ciphertext innerKSDownExt(innerKSDownExt_ax, innerKSDownExt_bx, N, slots, curr_limbs, special_limbs);
                    FastRotate_KeyGen(scheme.secretKey, rot_out[s][i], scheme);//todo: 将密钥生成挪到预计算中去

                    EvalFastRotationExt(innerKSDownExt, innerKSDown, rot_out[s][i], innerDigits, false, scheme);
                    scheme.EvalAddExtInPlace(outer, innerKSDownExt);
                }
                else {
                    //first += cc->KeySwitchDownFirstElement(inner);
                    uint64_t* tmp_first= new uint64_t [len];
                    KeySwitchDownFirstElement(tmp_first,inner.bx,curr_limbs,scheme);
                    scheme.context.addAndEqual(first,tmp_first,curr_limbs);//first += tmp_first;

                    // auto elements = inner->GetElements();
                    // elements[0].SetValuesToZero();
                    // inner->SetElements(elements);
                    SetZero(inner.bx,0,lenExt);

                    scheme.EvalAddExtInPlace(outer, inner);//EvalAddExtInPlace(outer, inner);
                }
            }
        }

        //result                          = cc->KeySwitchDown(outer);
        KeySwitchDown(result.ax, result.bx, outer.ax, outer.bx, curr_limbs, scheme);

        // std::vector<DCRTPoly>& elements = result->GetElements();
        // elements[0] += first;
        scheme.context.addAndEqual(result.bx,first,curr_limbs);
    }

    //return result;
}

// //DongYilin add
void EvalCoeffsToSlots(Ciphertext &result,
                       Plaintext** A,uint32_t A_len,
                       Ciphertext ctxt,
                       Scheme scheme)
{
    uint32_t slots = A_len;//uint32_t slots = ctxt->GetSlots();

    // auto pair = m_bootPrecomMap.find(slots);
    // if (pair == m_bootPrecomMap.end()) {
    //     std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
    //                          std::string(" slots were not generated") +
    //                          std::string(" Need to call EvalBootstrapSetup and EvalBootstrapKeyGen to proceed"));
    //     OPENFHE_THROW(type_error, errorMsg);
    // }
    // const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;

    // auto cc    = ctxt->GetCryptoContext();


    uint32_t N = scheme.context.N;//uint32_t N = cc->GetRingDimension();
    uint32_t M = scheme.context.M;//uint32_t M = cc->GetCyclotomicOrder();
    uint32_t K = scheme.context.K;
    uint32_t logN = scheme.context.logN;
    uint32_t special_limbs=K;
    CKKSBootstrapPrecom& precom=scheme.precom;


    int32_t levelBudget     = precom.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET];
    int32_t layersCollapse  = precom.m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_COLL];
    int32_t remCollapse     = precom.m_paramsEnc[CKKS_BOOT_PARAMS::LAYERS_REM];
    int32_t numRotations    = precom.m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS];
    int32_t b               = precom.m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP];
    int32_t g               = precom.m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP];
    int32_t numRotationsRem = precom.m_paramsEnc[CKKS_BOOT_PARAMS::NUM_ROTATIONS_REM];
    int32_t bRem            = precom.m_paramsEnc[CKKS_BOOT_PARAMS::BABY_STEP_REM];
    int32_t gRem            = precom.m_paramsEnc[CKKS_BOOT_PARAMS::GIANT_STEP_REM];

    int32_t stop    = -1;
    int32_t flagRem = 0;

    //auto algo = cc->GetScheme();

    if (remCollapse != 0) {
        stop    = 0;
        flagRem = 1;
    }

    // precompute the inner and outer rotations
    //std::vector<std::vector<int32_t>> rot_in(levelBudget);
    uint32_t** rot_in=new uint32_t*[levelBudget];

    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        if (flagRem == 1 && i == 0) {
            // remainder corresponds to index 0 in encoding and to last index in decoding
            rot_in[i]=new uint32_t[numRotationsRem+1];//rot_in[i] = std::vector<int32_t>(numRotationsRem + 1);
        }
        else {
            rot_in[i]=new uint32_t[numRotations+1];//rot_in[i] = std::vector<int32_t>(numRotations + 1);
        }
    }

    //std::vector<std::vector<int32_t>> rot_out(levelBudget);
    uint32_t** rot_out=new uint32_t*[levelBudget];

    for (uint32_t i = 0; i < uint32_t(levelBudget); i++) {
        rot_out[i] = new uint32_t[b+bRem];//rot_out[i] = std::vector<int32_t>(b + bRem);
    }

    //ReduceRotation在ckksrns-util.cpp中实现
    for (int32_t s = levelBudget - 1; s > stop; s--) {
        for (int32_t j = 0; j < g; j++) {
            rot_in[s][j] = ReduceRotation(
                    (j - int32_t((numRotations + 1) / 2) + 1) * (1 << ((s - flagRem) * layersCollapse + remCollapse)),
                    slots);
        }

        for (int32_t i = 0; i < b; i++) {
            rot_out[s][i] = ReduceRotation((g * i) * (1 << ((s - flagRem) * layersCollapse + remCollapse)), M / 4);
        }
    }

    if (flagRem) {
        for (int32_t j = 0; j < gRem; j++) {
            rot_in[stop][j] = ReduceRotation((j - int32_t((numRotationsRem + 1) / 2) + 1), slots);
        }

        for (int32_t i = 0; i < bRem; i++) {
            rot_out[stop][i] = ReduceRotation((gRem * i), M / 4);
        }
    }

    //Ciphertext<DCRTPoly> result = ctxt->Clone();
//    Ciphertext tmp_result=ctxt;
    result=ctxt;

    // hoisted automorphisms
    for (int32_t s = levelBudget - 1; s > stop; s--) {
        if (s != levelBudget - 1) {
            ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP, scheme);//algo->ModReduceInternalInPlace(resultExt, BASE_NUM_LEVELS_TO_DROP);
        }

        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        //auto digits = cc->EvalFastRotationPrecompute(resultExt);
        long curr_limbs = result.curr_limbs;
        long limbsExt = curr_limbs + special_limbs;
        long lenExt = (limbsExt<<logN);
        long len = (curr_limbs<<logN);
        uint32_t alpha = scheme.context.alpha;
        long beta = std::ceil((curr_limbs * 1.0 / alpha));
        
        long digits_len = beta * lenExt;
        uint64_t* digits = new uint64_t[digits_len];
        EvalFastRotationPrecompute(digits, result, scheme);

        Ciphertext* fastRotationExt=new Ciphertext[g];//std::vector<ciphertext> fastRotation(g);

//#pragma omp parallel for
        for (int32_t j = 0; j < g; j++) {
            fastRotationExt[j].N=N;
            fastRotationExt[j].slots=slots;
            fastRotationExt[j].curr_limbs=curr_limbs;
            fastRotationExt[j].special_limbs=special_limbs;
            fastRotationExt[j].ax = new uint64_t[lenExt]();
            fastRotationExt[j].bx = new uint64_t[lenExt]();
            if (rot_in[s][j] != 0) {
                //fastRotation[j] = cc->EvalFastRotationExt(resultExt, rot_in[s][j], digits, true);
                FastRotate_KeyGen(scheme.secretKey,rot_in[s][j],scheme);//todo: 将密钥生成挪到预计算中去
                EvalFastRotationExt(fastRotationExt[j], result, rot_in[s][j], digits, true, scheme);
            }
            else {
                KeySwitchExt(fastRotationExt[j], result, NORMAL_CIPHER_SIZE, true, scheme);//fastRotation[j] = cc->KeySwitchExt(resultExt, true);
            }
        }
        delete[]digits;


        //Ciphertext<DCRTPoly> outer;
        uint64_t* outer_ax = new uint64_t[lenExt]();
        uint64_t* outer_bx = new uint64_t[lenExt]();
        Ciphertext outer(outer_ax,outer_bx,N,slots,curr_limbs,special_limbs);
        //DCRTPoly first;
        uint64_t* first= new uint64_t [len];

        for (int32_t i = 0; i < b; i++) {
            // for the first iteration with j=0:
            int32_t G                  = g * i;
            //Ciphertext<DCRTPoly> inner= EvalMultExt(fastRotation[0], A[s][G]);
            uint64_t *inner_ax = new uint64_t[lenExt]();
            uint64_t *inner_bx = new uint64_t[lenExt]();
            Ciphertext inner(inner_ax, inner_bx, N, slots, curr_limbs, special_limbs);
            scheme.EvalMultExt(inner, fastRotationExt[0], A[s][G]);

            // continue the loop
            for (int32_t j = 1; j < g; j++) {
                if ((G + j) != int32_t(numRotations)) {
                    uint64_t *tmp_ax = new uint64_t[lenExt]();
                    uint64_t *tmp_bx = new uint64_t[lenExt]();
                    Ciphertext tmpExt(tmp_ax, tmp_bx, N, slots, curr_limbs, special_limbs);
                    //EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[s][G + j]));
                    scheme.EvalMultExt(tmpExt, fastRotationExt[j], A[s][G + j]);
                    scheme.EvalAddExtInPlace(inner, tmpExt);
                }
            }

            if (i == 0) {
                //first         = cc->KeySwitchDownFirstElement(inner);
                KeySwitchDownFirstElement(first,inner.bx,curr_limbs,scheme);
                // auto elements = inner->GetElements();
                // elements[0].SetValuesToZero();
                // inner->SetElements(elements);
                SetZero(inner.bx,0,lenExt);

                outer = inner;
            }
            else {
                if (rot_out[s][i] != 0) {
                    // inner = cc->KeySwitchDown(inner);
                    // => innerKSDown = cc->KeySwitchDown(inner);
                    uint32_t innerKSDown_len= (curr_limbs << logN);;
                    uint64_t* innerKSDown_ax = new uint64_t[innerKSDown_len]();
                    uint64_t* innerKSDown_bx = new uint64_t[innerKSDown_len]();
                    Ciphertext innerKSDown(innerKSDown_ax, innerKSDown_bx, ctxt.N, ctxt.slots, curr_limbs);
                    KeySwitchDown(innerKSDown.ax, innerKSDown.bx, inner.ax, inner.bx, curr_limbs, scheme);

                    // Find the automorphism index that corresponds to rotation index.
                    long autoIndex = FindAutomorphismIndex2nComplex(rot_out[s][i], M);
                    uint32_t map_len=N;
                    uint32_t * map  = new uint32_t [map_len]();
                    PrecomputeAutoMap(N, autoIndex, map);

                    //first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    uint64_t * firstCurrent= new uint64_t[len];
                    AutomorphismTransform(firstCurrent, innerKSDown.bx, curr_limbs, N, autoIndex, map);
                    scheme.context.addAndEqual(first,firstCurrent,curr_limbs);//first += firstCurrent;

                    //auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    uint64_t* innerDigits = new uint64_t[digits_len];
                    EvalFastRotationPrecompute(innerDigits, innerKSDown, scheme);

                    //EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[s][i], innerDigits, false));
                    uint64_t *innerKSDownExt_ax = new uint64_t[lenExt]();
                    uint64_t *innerKSDownExt_bx = new uint64_t[lenExt]();
                    Ciphertext innerKSDownExt(innerKSDownExt_ax, innerKSDownExt_bx, ctxt.N, ctxt.slots,curr_limbs, special_limbs);
                    FastRotate_KeyGen(scheme.secretKey, rot_out[s][i], scheme);//todo: 将密钥生成挪到预计算中去

                    EvalFastRotationExt(innerKSDownExt, innerKSDown, rot_out[s][i], innerDigits, false, scheme);
                    scheme.EvalAddExtInPlace(outer, innerKSDownExt);
                }
                else {
                    //first += cc->KeySwitchDownFirstElement(inner);
                    uint64_t * tmp_first= new uint64_t[len];
                    KeySwitchDownFirstElement(tmp_first,inner.bx,curr_limbs,scheme);
                    scheme.context.addAndEqual(first,tmp_first,curr_limbs);//first += tmp_first;

                    // auto elements = inner->GetElements();
                    // elements[0].SetValuesToZero();
                    // inner->SetElements(elements);
                    SetZero(inner.bx,0,lenExt);

                    scheme.EvalAddExtInPlace(outer, inner);//EvalAddExtInPlace(outer, inner);
                }
            }
        }
        //result                          = cc->KeySwitchDown(outer);
        KeySwitchDown(result.ax, result.bx, outer.ax, outer.bx, curr_limbs, scheme);

        // std::vector<DCRTPoly>& elements = result->GetElements();
        // elements[0] += first;
        scheme.context.addAndEqual(result.bx, first, curr_limbs);
    }
    
    if (flagRem) {
        ModReduceInternalInPlace(result, BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(resultExt, BASE_NUM_LEVELS_TO_DROP);

        // computes the NTTs for each CRT limb (for the hoisted automorphisms used later on)
        //auto digits = cc->EvalFastRotationPrecompute(result);
        long curr_limbs = result.curr_limbs;
        long limbsExt = curr_limbs + special_limbs;
        long lenExt = (limbsExt<<logN);
        long len = (curr_limbs<<logN);
        uint32_t alpha = scheme.context.alpha;
        long beta = std::ceil((curr_limbs * 1.0 / alpha));

        long digits_len = beta * lenExt;
        uint64_t* digits = new uint64_t[digits_len];
        EvalFastRotationPrecompute(digits, result, scheme);

        Ciphertext* fastRotationExt=new Ciphertext[gRem];//std::vector<Ciphertext<DCRTPoly>> fastRotation(gRem);

//#pragma omp parallel for
        for (int32_t j = 0; j < gRem; j++) {
            fastRotationExt[j].N=N;
            fastRotationExt[j].slots=slots;
            fastRotationExt[j].curr_limbs=curr_limbs;
            fastRotationExt[j].special_limbs=special_limbs;
            fastRotationExt[j].ax = new uint64_t[lenExt]();
            fastRotationExt[j].bx = new uint64_t[lenExt]();
            if (rot_in[stop][j] != 0) {
                //fastRotation[j] = cc->EvalFastRotationExt(resultExt, rot_in[stop][j], digits, true);
                FastRotate_KeyGen(scheme.secretKey,rot_in[stop][j],scheme);//先生成密钥
                EvalFastRotationExt(fastRotationExt[j], result, rot_in[stop][j], digits, true, scheme);
            }
            else {
                KeySwitchExt(fastRotationExt[j], result, NORMAL_CIPHER_SIZE, true, scheme);//fastRotation[j] = cc->KeySwitchExt(resultExt, true);
            }
        }
        delete[]digits;


        Ciphertext outer;//Ciphertext<DCRTPoly> outer;
        uint64_t* first= new uint64_t[len];//DCRTPoly first;

        for (int32_t i = 0; i < bRem; i++) {

            // for the first iteration with j=0:
            int32_t GRem = gRem * i;
            //Ciphertext<DCRTPoly> inner;
            uint64_t *inner_ax = new uint64_t[lenExt]();
            uint64_t *inner_bx = new uint64_t[lenExt]();
            Ciphertext inner(inner_ax, inner_bx, N, slots, curr_limbs, special_limbs);
            //inner        = EvalMultExt(fastRotation[0], A[stop][GRem]);
            scheme.EvalMultExt(inner, fastRotationExt[0], A[stop][GRem]);

            // continue the loop
            for (int32_t j = 1; j < gRem; j++) {
                if ((GRem + j) != int32_t(numRotationsRem)) {
                    //EvalAddExtInPlace(inner, EvalMultExt(fastRotation[j], A[stop][GRem + j]));
                   
                    uint64_t *tmp_ax = new uint64_t[lenExt]();
                    uint64_t *tmp_bx = new uint64_t[lenExt]();
                    Ciphertext tmpExt(tmp_ax, tmp_bx, ctxt.N, ctxt.slots, curr_limbs, special_limbs);
                    scheme.EvalMultExt(tmpExt, fastRotationExt[j], A[stop][GRem + j]);//todo: check correctness
                    scheme.EvalAddExtInPlace(inner, tmpExt);
                }
            }

            if (i == 0) {
                //first         = cc->KeySwitchDownFirstElement(inner);
                KeySwitchDownFirstElement(first,inner.bx,curr_limbs,scheme);

                // auto elements = inner->GetElements();
                // elements[0].SetValuesToZero();
                // inner->SetElements(elements);
                SetZero(inner.bx,0,lenExt);

                outer = inner;
            }
            else {
                if (rot_out[stop][i] != 0) {
                    //inner = cc->KeySwitchDown(inner);
                    uint32_t innerKSDown_len= (curr_limbs << logN);;
                    uint64_t* innerKSDown_ax = new uint64_t[innerKSDown_len]();
                    uint64_t* innerKSDown_bx = new uint64_t[innerKSDown_len]();
                    Ciphertext innerKSDown(innerKSDown_ax, innerKSDown_bx, ctxt.N, ctxt.slots, curr_limbs);
                    KeySwitchDown(innerKSDown.ax, innerKSDown.bx, inner.ax, inner.bx, curr_limbs, scheme);

                    // Find the automorphism index that corresponds to rotation index.
                    long autoIndex = FindAutomorphismIndex2nComplex(rot_out[stop][i], M);
                    uint32_t map_len=N;
                    uint32_t * map  = new uint32_t[map_len]();
                    PrecomputeAutoMap(N, autoIndex, map);

                    //first += inner->GetElements()[0].AutomorphismTransform(autoIndex, map);
                    uint64_t * firstCurrent= new uint64_t[len];
                    AutomorphismTransform(firstCurrent, innerKSDown.bx, curr_limbs, N, autoIndex, map);
                    scheme.context.addAndEqual(first,firstCurrent,curr_limbs);//first += firstCurrent;

                    //auto innerDigits = cc->EvalFastRotationPrecompute(inner);
                    uint64_t* innerDigits = new uint64_t[digits_len];
                    EvalFastRotationPrecompute(innerDigits, innerKSDown, scheme);

                    uint64_t *innerKSDownExt_ax = new uint64_t[lenExt]();
                    uint64_t *innerKSDownExt_bx = new uint64_t[lenExt]();
                    Ciphertext innerKSDownExt(innerKSDownExt_ax, innerKSDownExt_bx, N, slots, curr_limbs, special_limbs);
                    FastRotate_KeyGen(scheme.secretKey, rot_out[stop][i], scheme);//todo: 将密钥生成挪到预计算中去
                    //EvalAddExtInPlace(outer, cc->EvalFastRotationExt(inner, rot_out[stop][i], innerDigits, false));
                    EvalFastRotationExt(innerKSDownExt, innerKSDown, rot_out[stop][i], innerDigits, false, scheme);
                    scheme.EvalAddExtInPlace(outer, innerKSDownExt);
                }
                else {
                    //first += cc->KeySwitchDownFirstElement(inner);
                    uint64_t* tmp_first= new uint64_t [len];
                    KeySwitchDownFirstElement(tmp_first,inner.bx,curr_limbs,scheme);
                    scheme.context.addAndEqual(first,tmp_first,curr_limbs);//first += tmp_first;

                    // auto elements = inner->GetElements();
                    // elements[0].SetValuesToZero();
                    // inner->SetElements(elements);
                    SetZero(inner.bx,0,lenExt);

                    scheme.EvalAddExtInPlace(outer, inner);//EvalAddExtInPlace(outer, inner);
                }
            }
        }

        //result                          = cc->KeySwitchDown(outer);
        KeySwitchDown(result.ax, result.bx, outer.ax, outer.bx, curr_limbs, scheme);

        // std::vector<DCRTPoly>& elements = result->GetElements();
        // elements[0] += first;
        scheme.context.addAndEqual(result.bx, first, curr_limbs);
    }

    //return result;
}

