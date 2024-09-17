//
// Created by EYx on 2023/5/5.
//

//Choose one for LT
//#include "../data/N64slots2isLTBootstrap1isSparse1lEnc16lDec2.h"
//#include "../data/N64slots8isLTBootstrap1isSparse1lEnc16lDec2.h"
//#include "../data/N64slots16isLTBootstrap1isSparse1lEnc16lDec2.h"
#include "../data/N64slots32isLTBootstrap1isSparse0lEnc16lDec2.h"

//Choose one for C2S, S2C
#include "../data/N64slots32isLTBootstrap0isSparse0lEnc19lDec2.h"


#include "../data/BootstrapKeyGen_data.h"
#include "../data/BSConstValue.h"

#include "Common.h"
#include "Context.h"
#include "Scheme.h"
#include "Ciphertext.h"
#include "Plaintext.h"
#include "StringUtils.h"

#include "myBootstrapExample.h"
#include "myUtils.h"
#include "myRotation.h"
#include "myLT.h"
#include "myChebyshevFunc.h"
#include "constants.h"
#include "ckksrns-leveledshe.h"//使用其中的EvalMultInPlace
#include "ckksrns-fhe.h"

#include <sstream>// 引入字符串流库

//#define DEBUG_EvalBootstrap

//void BootstrapExampleClean(uint32_t n, uint32_t slots, uint32_t levelsRemaining, SecretKeyDist secretKeyDist)
void BootstrapExampleClean(uint32_t n, uint32_t slots, uint32_t levelsRemaining, SecretKeyDist secretKeyDist)
{
    // giant step for baby-step-giant-step algorithm in linear transforms for encoding and decoding, respectively
    // Choose this a power of 2 preferably, otherwise an exact divisor of the number of elements in the sum
    std::uint32_t dim1[2] = {0, 0};

    // budget in levels for FFT for encoding and decoding, respectively
    // Choose a number smaller than ceil(log2(slots))
    std::uint32_t levelBudget[2] = {4, 4};
//    std::uint32_t levelBudget[2] = {1,1};

    // computes how many levels are needed for
    uint32_t depth = levelsRemaining + GetBootstrapDepth(9, levelBudget, secretKeyDist);
    long L0 = depth + 1;
    long dnum = 1; //Note: 满足dnum=1的模数pi全体，可以用作任意dnum>1时的计算
    long K = L0 / dnum;
    long L = L0;
    long logN = 6;
    // uint32_t firstMod = 60; /*firstMod*/ //Note: 对应本工程的 Q0_BIT_SIZE
    long logp = 59;//openfhe：uint32_t dcrtBits = 59;


    Context context("Debug",logN, logp, L, K, dnum);
    SecretKey secretKey(context);
    ScalingTechnique rescaleTech = FIXEDMANUAL;
    Scheme scheme(secretKey, context, rescaleTech, secretKeyDist, L0);//NOTE: 此处应按照L0个limbs这一峰值生成secretKey

    std::complex<double> *mvec= new complex<double>[slots];
    for (int i = 0; i < slots; ++i) {
//        mvec[i]=0.111111*i;//correct
//        mvec[i]=-0.11111*i;//correct
        mvec[i]=-1.5;//correct
    }
    //NOTE:  bootstrap 第一步 AdjustCiphertext 中需要做一次明密文乘，所以这里生成Limbs=2的密文
    Ciphertext cipher = scheme.fullencrypt(mvec, slots, 2);
    StringUtils::FullDecAndShow(scheme, cipher, "Before Bootstrap");

//    //Debug Bootstrap
//    Context context("Debug", logN, logp, L, K, dnum);//写死原根和模数
//    SecretKey secretKey(context, "Bootstrap");
//    SET_SWK=1;
//    ScalingTechnique rescaleTech = FIXEDMANUAL;
//    Scheme scheme(secretKey, context, rescaleTech, secretKeyDist, L0);
//    std::complex<double> *mvec= new complex<double>[slots];
//    for (int i = 0; i < slots; ++i) {
////        mvec[i]=0.111111*i;//correct
////        mvec[i]=-0.11111*i;//correct
//        mvec[i]=-1.5;//correct
//    }
//    Ciphertext cipher(new uint64_t[(2<<logN)],new uint64_t[(2<<logN)],(1<< logN),slots,2);
//    //使用Openfhe输出的密文多项式
//    //SET 密文初始值
//    for (int i = 0; i < 2; ++i) {
//        uint64_t *axj = cipher.ax + (i << logN);
//        uint64_t *bxj = cipher.bx + (i << logN);
//        for (int j = 0; j < context.N; ++j) {
//            bxj[j] = BS_NegCT[0][i][j];
//            axj[j] = BS_NegCT[1][i][j];
//        }
//    }
//    StringUtils::FullDecAndShow(scheme, cipher, "Before Bootstrap");

    EvalBootstrapSetup(scheme, levelBudget, dim1, slots);

    EvalBootstrapKeyGen(scheme, secretKey, slots);

    Ciphertext result = cipher;
    EvalBootstrap(scheme, result, cipher);

    complex<double> *dvec = scheme.fulldecrypt(secretKey, result);
    StringUtils::showcompare(mvec, dvec, slots, "Bootstrap_demo");

}

uint32_t GetBootstrapDepth(uint32_t approxModDepth, uint32_t* levelBudget,
                           SecretKeyDist secretKeyDist)
{
    if (secretKeyDist == UNIFORM_TERNARY) {
        approxModDepth += R - 1;
    }

    return approxModDepth + levelBudget[0] + levelBudget[1];
}

void FindBootstrapRotationIndices(Scheme scheme, int32_t* fullIndexList, long fullIndexList_len, uint32_t slots, uint32_t M) {

    //fixme: change to on the fly compute
    for (int i = 0; i < fullIndexList_len; ++i) {
        fullIndexList[i]=fullIndexList_N64slots8isLT1[i];
    }

//    auto pair = m_bootPrecomMap.find(slots);
//    if (pair == m_bootPrecomMap.end()) {
//        std::string errorMsg(std::string("Precomputations for ") + std::to_string(slots) +
//                             std::string(" slots were not generated") +
//                             std::string(" Need to call EvalBootstrapSetup to proceed"));
//        OPENFHE_THROW(type_error, errorMsg);
//    }
//    const std::shared_ptr<CKKSBootstrapPrecom> precom = pair->second;
//
//    std::vector<int32_t> fullIndexList;
//
//    bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
//                         (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);
//
//    if (isLTBootstrap) {
//        std::vector<int32_t> indexList = FindLinearTransformRotationIndices(slots, M);
//        fullIndexList.insert(fullIndexList.end(), indexList.begin(), indexList.end());
//    }
//    else {
//        std::vector<int32_t> indexListCtS = FindCoeffsToSlotsRotationIndices(slots, M);
//        std::vector<int32_t> indexListStC = FindSlotsToCoeffsRotationIndices(slots, M);
//
//        fullIndexList.insert(fullIndexList.end(), indexListCtS.begin(), indexListCtS.end());
//        fullIndexList.insert(fullIndexList.end(), indexListStC.begin(), indexListStC.end());
//    }
//
//    // Remove possible duplicates
//    sort(fullIndexList.begin(), fullIndexList.end());
//    fullIndexList.erase(unique(fullIndexList.begin(), fullIndexList.end()), fullIndexList.end());
//
//    // remove automorphisms corresponding to 0
//    fullIndexList.erase(std::remove(fullIndexList.begin(), fullIndexList.end(), 0), fullIndexList.end());
//    fullIndexList.erase(std::remove(fullIndexList.begin(), fullIndexList.end(), M / 4), fullIndexList.end());
//
//    return fullIndexList;
}

void AdjustCiphertext(Scheme scheme, Ciphertext& ciphertext, double correction)
{
    ScalingTechnique rescaleTech=scheme.rescaleTech;

//    if (cryptoParams->GetScalingTechnique() == FLEXIBLEAUTO || cryptoParams->GetScalingTechnique() == FLEXIBLEAUTOEXT)
    if (rescaleTech == FLEXIBLEAUTO || rescaleTech == FLEXIBLEAUTOEXT){
        //TODO: to be implemented
    }
    else{
        //#if NATIVEINT != 128
        // Scaling down the message by a correction factor to emulate using a larger q0.
        // This step is needed so we could use a scaling factor of up to 2^59 with q9 ~= 2^60.
        //cc->EvalMultInPlace(ciphertext, std::pow(2, -correction));
        //algo->ModReduceInternalInPlace(ciphertext, BASE_NUM_LEVELS_TO_DROP);
        auto cnst=std::pow(2, -correction);

        EvalMultInPlace(ciphertext,cnst,scheme);

#ifdef DEBUG_EvalBootstrap
        StringUtils::CompareAndSet(2,64,ciphertext,
                                   reinterpret_cast<uint64_t *>(after_mult),
                                   "after_mult ");
#endif

        ModReduceInternalInPlace(ciphertext,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(ciphertext, BASE_NUM_LEVELS_TO_DROP);

#ifdef DEBUG_EvalBootstrap
        StringUtils::CompareAndSet(1,64,ciphertext,
                                   reinterpret_cast<uint64_t *>(after_ModReduce),
                                   "after_ModReduce ");
#endif
    }
}

void SwitchModulus(uint64_t* res, uint64_t newModulus, uint64_t* input, uint64_t oldModulus, uint64_t N)
{
    uint64_t oldModulusByTwo=oldModulus>>1;
    uint64_t diff=(oldModulus > newModulus) ? (oldModulus - newModulus) : (newModulus - oldModulus);

    if (newModulus > oldModulus) {
        for (uint64_t i = 0; i < N; i++) {
            uint64_t n = input[i];//IntegerType n = this->m_data[i];
            if (n > oldModulusByTwo) {
                res[i]= input[i] + diff;//this->m_data[i] += diff;
            }
        }
    }
    else {  // newModulus <= oldModulus
        for (uint64_t i = 0; i < N; i++) {
            uint64_t n = input[i];//IntegerType n = this->m_data[i];
            uint64_t sub_diff = (n > oldModulusByTwo) ? diff : 0;
            subMod(res[i],n,sub_diff,newModulus); //this->m_data[i]      = n.ModSub(sub_diff, newModulus);
        }
    }

}

void EvalBootstrapSetup(Scheme &scheme, uint32_t* levelBudget, uint32_t* dim1, uint32_t numslots, uint32_t correctionFactor)
{
    uint32_t M=scheme.context.M;
    uint32_t N=scheme.context.N;
    uint32_t slots = (numslots==0)? M/4: numslots;
    ScalingTechnique rescaleTech=scheme.rescaleTech;
    SecretKeyDist secretKeyDist=scheme.secretKeyDist;
    CKKSBootstrapPrecom& precom=scheme.precom;
    uint64_t* qVec=scheme.context.qVec;
    uint32_t L0 = scheme.L0;//L0 is the total level after ModRaise

    // Set correction factor by default, if it is not already set.
    if (correctionFactor == 0) {
        if ( rescaleTech== FLEXIBLEAUTO ||
                rescaleTech == FLEXIBLEAUTOEXT) {
            // The default correction factors chosen yielded the best precision in our experiments.
            // We chose the best fit line from our experiments by running ckks-bootstrapping-precision.cpp.
            // The spreadsheet with our experiments is here:
            // https://docs.google.com/spreadsheets/d/1WqmwBUMNGlX6Uvs9qLXt5yeddtCyWPP55BbJPu5iPAM/edit?usp=sharing
            auto tmp = std::round(-0.265 * (2 * std::log2(M / 2) + std::log2(slots)) + 19.1);
            if (tmp < 7)
                m_correctionFactor = 7;
            else if (tmp > 13)
                m_correctionFactor = 13;
            else
                m_correctionFactor = static_cast<uint32_t>(tmp);
        }
        else {
            m_correctionFactor = 9;
        }
    }
    else {
        m_correctionFactor = correctionFactor;
    }

    precom.m_slots=slots;
    precom.m_dim1=dim1[0];

    double logSlots = std::log2(slots);  // TODO (dsuponit): can logSlots be cast to uint32_t on this line?
    // Perform some checks on the level budget and compute parameters
    uint32_t newBudget[2] = {levelBudget[0],levelBudget[1]};

    if (levelBudget[0] > logSlots) {
        std::cerr << "\nWarning, the level budget for encoding cannot be this large. The budget was changed to "
                  << uint32_t(logSlots) << std::endl;
        newBudget[0] = uint32_t(logSlots);
    }
    if (levelBudget[0] < 1) {
        std::cerr << "\nWarning, the level budget for encoding has to be at least 1. The budget was changed to " << 1
                  << std::endl;
        newBudget[0] = 1;
    }

    if (levelBudget[1] > logSlots) {
        std::cerr << "\nWarning, the level budget for decoding cannot be this large. The budget was changed to "
                  << uint32_t(logSlots) << std::endl;
        newBudget[1] = uint32_t(logSlots);
    }
    if (levelBudget[1] < 1) {
        std::cerr << "\nWarning, the level budget for decoding has to be at least 1. The budget was changed to " << 1
                  << std::endl;
        newBudget[1] = 1;
    }

    precom.m_paramsEnc = GetCollapsedFFTParams(slots, newBudget[0], dim1[0]);
    precom.m_paramsDec = GetCollapsedFFTParams(slots, newBudget[1], dim1[1]);


    //fixme: change to on-the-fly compute
    precom.m_U0Pre=new Plaintext [LTMatrix_Row];
    precom.m_U0hatTPre=new Plaintext [LTMatrix_Row];
    for (int i = 0; i < LTMatrix_Row; ++i) {
        //precom.m_U0hatTPre
        int m_U0hatTPre_len=LTMatrix_mx_len*m_U0hatTPre_limbs;
        uint64_t *m_U0hatTPre=new uint64_t[m_U0hatTPre_len];
        for (int j = 0; j < m_U0hatTPre_len; ++j) {
            m_U0hatTPre[j]=m_U0hatTPre_mx[i][j];
        }
        precom.m_U0hatTPre[i] = Plaintext(m_U0hatTPre,LTMatrix_mx_len,LTMatrix_Column,m_U0hatTPre_limbs);

        //precom.m_U0Pre
        int m_U0Pre_len=LTMatrix_mx_len*m_U0Pre_limbs;
        uint64_t *m_U0Pre=new uint64_t[m_U0Pre_len];
        for (int j = 0; j < m_U0Pre_len; ++j) {
            m_U0Pre[j]=m_U0Pre_mx[i][j];
        }
        precom.m_U0Pre[i] = Plaintext(m_U0Pre,LTMatrix_mx_len,LTMatrix_Column,m_U0Pre_limbs);
    }

    int RHScnt=0;
    precom.m_U0hatTPreFFT=new Plaintext* [m_U0hatTPreFFT_dim1];
    for(int i=0;i<m_U0hatTPreFFT_dim1;i++){
        long j_len=m_U0hatTPreFFT_dim2[i];
        precom.m_U0hatTPreFFT[i]=new Plaintext [j_len];
        int limbs=m_U0hatTPreFFT_limbs[i];
        int m_U0hatTPreFFT_len=mx_len*limbs;
        for (int j = 0; j < j_len; ++j) {
            uint64_t *m_U0hatTPreFFT=new uint64_t[m_U0hatTPreFFT_len];
            int LHScnt=0;
            for (int k = 0; k < limbs; ++k) {
                for (int l = 0; l < mx_len; ++l) {
                    m_U0hatTPreFFT[LHScnt]=m_U0hatTPreFFT_mx[RHScnt];
                    LHScnt++;
                    RHScnt++;
                }
            }
            precom.m_U0hatTPreFFT[i][j] =
            Plaintext(m_U0hatTPreFFT,mx_len,mx_slots,limbs);
        }
    }


    RHScnt=0;
    precom.m_U0PreFFT=new Plaintext* [m_U0PreFFT_dim1];
    for(int i=0;i<m_U0PreFFT_dim1;i++){
        long j_len=m_U0PreFFT_dim2[i];
        precom.m_U0PreFFT[i]=new Plaintext [j_len];
        int limbs=m_U0PreFFT_limbs[i];
        int m_U0PreFFT_len=mx_len*limbs;
        for (int j = 0; j < j_len; ++j) {
            uint64_t *m_U0PreFFT=new uint64_t[m_U0PreFFT_len];
            int LHScnt=0;
            for (int k = 0; k < limbs; ++k) {
                for (int l = 0; l < mx_len; ++l) {
                    m_U0PreFFT[LHScnt]=m_U0PreFFT_mx[RHScnt];
                    LHScnt++;
                    RHScnt++;
                }
            }
            precom.m_U0PreFFT[i][j] =
                    Plaintext(m_U0PreFFT,mx_len,mx_slots,limbs);
        }
    }

//    uint32_t m    = 4 * slots;
//    bool isSparse = (M != m) ? true : false;
//
//    // computes indices for all primitive roots of unity
//    std::vector<uint32_t> rotGroup(slots);
//    uint32_t fivePows = 1;
//    for (uint32_t i = 0; i < slots; ++i) {
//        rotGroup[i] = fivePows;
//        fivePows *= 5;
//        fivePows %= m;
//    }
//
//    // computes all powers of a primitive root of unity exp(2 * M_PI/m)
//    std::vector<std::complex<double>> ksiPows(m + 1);
//    for (uint32_t j = 0; j < m; ++j) {
//        double angle = 2.0 * M_PI * j / m;
//        ksiPows[j].real(cos(angle));
//        ksiPows[j].imag(sin(angle));
//    }
//    ksiPows[m] = ksiPows[0];
//
//    // Extract the modulus prior to bootstrapping
//    uint64_t q=qVec[0];//NativeInteger q = cryptoParams->GetElementParams()->GetParams()[0]->GetModulus().ConvertToInt();
//    double qDouble  = ConvertToDouble(q);
//
//    unsigned __int128 factor = ((unsigned __int128)1 << ((uint32_t)std::round(std::log2(qDouble))));
//    double pre               = qDouble / factor;
//    double k                 = (secretKeyDist == SPARSE_TERNARY) ? K_SPARSE : 1.0;
//    double scaleEnc          = pre / k;
//    double scaleDec          = 1 / pre;
//
//    uint32_t approxModDepth = 8;
//    if (secretKeyDist == UNIFORM_TERNARY) {
//        approxModDepth += R - 1;
//    }
//
//    uint32_t depthBT = approxModDepth + 1 + precom.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] +
//                       precom.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET];
//
//    // compute # of levels to remain when encoding the coefficients
//    //uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();//L0 is the total level after ModRaise
//    // for FLEXIBLEAUTOEXT we do not need extra modulus in auxiliary plaintexts
//    if (rescaleTech == FLEXIBLEAUTOEXT)
//        L0 -= 1;
//    uint32_t lEnc = L0 - precom.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] - 1;
//    uint32_t lDec = L0 - depthBT;
//
//    bool isLTBootstrap = (precom.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
//                         (precom.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);
//
//    if (isLTBootstrap) {
//        // allocate all vectors
//        std::vector<std::vector<std::complex<double>>> U0(slots, std::vector<std::complex<double>>(slots));
//        std::vector<std::vector<std::complex<double>>> U1(slots, std::vector<std::complex<double>>(slots));
//        std::vector<std::vector<std::complex<double>>> U0hatT(slots, std::vector<std::complex<double>>(slots));
//        std::vector<std::vector<std::complex<double>>> U1hatT(slots, std::vector<std::complex<double>>(slots));
//
//        for (size_t i = 0; i < slots; i++) {
//            for (size_t j = 0; j < slots; j++) {
//                U0[i][j]     = ksiPows[(j * rotGroup[i]) % m];
//                U0hatT[j][i] = std::conj(U0[i][j]);
//                U1[i][j]     = std::complex<double>(0, 1) * U0[i][j];
//                U1hatT[j][i] = std::conj(U1[i][j]);
//            }
//        }
//
////        if (!isSparse) {
////            precom.m_U0hatTPre = EvalLinearTransformPrecompute(scheme, U0hatT, scaleEnc, lEnc);
////            precom.m_U0Pre     = EvalLinearTransformPrecompute(scheme, U0, scaleDec, lDec);
////        }
////        else {
////            precom.m_U0hatTPre = EvalLinearTransformPrecompute(scheme, U0hatT, U1hatT, 0, scaleEnc, lEnc);
////            precom.m_U0Pre     = EvalLinearTransformPrecompute(scheme, U0, U1, 1, scaleDec, lDec);
////        }
//    }
//    else {
////        precom.m_U0hatTPreFFT = EvalCoeffsToSlotsPrecompute(scheme, ksiPows, rotGroup, false, scaleEnc, lEnc);
////        precom.m_U0PreFFT     = EvalSlotsToCoeffsPrecompute(scheme, ksiPows, rotGroup, false, scaleDec, lDec);
//    }
}

void EvalBootstrapKeyGen(Scheme &scheme, SecretKey secretKey, uint32_t slots)
{
    uint32_t M = scheme.context.M;

    if (slots == 0)
        slots = M / 4;

    //fixme: on-the-fly computation
    long fullIndexList_len=6;
    int32_t* fullIndexList=new int32_t [fullIndexList_len];
    FindBootstrapRotationIndices(scheme, fullIndexList,fullIndexList_len,slots,M);

    FastRotate_KeyGen_List(secretKey,fullIndexList,fullIndexList_len,scheme);

    scheme.addConjKey(secretKey);
}

//void EvalBootstrap(Scheme scheme, Ciphertext ciphertext,uint32_t numIterations=1,uint32_t precision);
void EvalBootstrap(Scheme scheme, Ciphertext& result, Ciphertext ciphertext,uint32_t numIterations, uint32_t precision)
{
    uint32_t M = scheme.context.M;
    uint32_t N = scheme.context.N;
    uint32_t logN = scheme.context.logN;
    uint32_t L0 = scheme.L0;
    uint32_t slots = ciphertext.slots;
    CKKSBootstrapPrecom precom = scheme.precom;
    ScalingTechnique rescaleTech = scheme.rescaleTech;
    SecretKeyDist secretKeyDist = scheme.secretKeyDist;
    uint64_t *moduliQ = scheme.context.qVec;
//    uint64_t *rootsQ=scheme.context.qRoots;
//    uint32_t sizeQ=scheme.L0;
//    uint32_t K=scheme.context.K;

    if (numIterations > 1) {
        //TODO: to be implemented
    }

    if (rescaleTech == FLEXIBLEAUTOEXT) {
        //TODO: to be implemented
        //elementParamsRaised.PopLastParam();
    }

    uint64_t q = moduliQ[0];
    double qDouble = ConvertToDouble(q);

    auto p = scheme.context.logp;//!!: 对应OpenFHE中的dcrbits //const auto p = cryptoParams->GetPlaintextModulus();
    double powP = pow(2, p);

    int32_t deg = std::round(std::log2(qDouble / powP));

    uint32_t correction = m_correctionFactor - deg;
    double post = std::pow(2, static_cast<double>(deg));

    double pre = 1. / post;
    uint64_t scalar = std::llround(post);

    //------------------------------------------------------------------------------
    // RAISING THE MODULUS
    //------------------------------------------------------------------------------

    // In FLEXIBLEAUTO, raising the ciphertext to a larger number
    // of towers is a bit more complex, because we need to adjust
    // it's scaling factor to the one that corresponds to the level
    // it's being raised to.
    // Increasing the modulus

    //Ciphertext<DCRTPoly> raised = ciphertext->Clone();
    //openFHE，在Bootstrap开始之前，处理未做的rescale操作//Note: 对于FIXEDMANUAL来说，下面这步不执行
//    algo->ModReduceInternalInPlace(raised, raised->GetNoiseScaleDeg() - 1);
    Ciphertext tmp = ciphertext;

    AdjustCiphertext(scheme, tmp, correction);

    uint64_t *raised_ax = new uint64_t[L0 << logN];
    uint64_t *raised_bx = new uint64_t[L0 << logN];
    Ciphertext raised(raised_ax, raised_bx, ciphertext.N, ciphertext.slots, L0);
    //Copy the first RNS poly of `Ciphertext tmp`, after `AdjustCiphertext`
    for (int i = 0; i < N; ++i) {
        raised.ax[i] = tmp.ax[i];
        raised.bx[i] = tmp.bx[i];
    }

    // We only use the level 0 ciphertext here. All other towers are automatically ignored to make
    // CKKS bootstrapping faster.

    //deal with cipher.bx
    scheme.context.INTTAndEqual(raised.bx, 1); //ctxtDCRT[i].SetFormat(COEFFICIENT);
    for (int i = 1; i < L0; ++i) {
        //        ctxtDCRT[i].SetFormat(COEFFICIENT);
        //        temp = ctxtDCRT[i].GetElementAtIndex(0);
        uint64_t *raised_bxj = raised.bx + (i << logN);
        SwitchModulus(raised_bxj, moduliQ[i], raised.bx, moduliQ[0], N);
    }
    scheme.context.NTTAndEqual(raised.bx, L0);     //        temp.SetFormat(EVALUATION);

    //deal with cipher.ax
    scheme.context.INTTAndEqual(raised.ax, 1);
    for (int i = 1; i < L0; ++i) {
        uint64_t *raised_axj = raised.ax + (i << logN);
        SwitchModulus(raised_axj, moduliQ[i], raised.ax, moduliQ[0], N);
    }
    scheme.context.NTTAndEqual(raised.ax, L0);

#ifdef DEBUG_EvalBootstrap
    StringUtils::DecAndShow(scheme,raised,"Dec after ModRaise");
    StringUtils::FullDecAndShow(scheme,raised,"Full Dec after ModRaise");
    StringUtils::CompareAndSet(raised.curr_limbs,N,raised,
                               reinterpret_cast<uint64_t *>(After_MultByCnst[0]),
                               "after multByConstAndEqual");
#endif

    //------------------------------------------------------------------------------
    // SETTING PARAMETERS FOR APPROXIMATE MODULAR REDUCTION
    //------------------------------------------------------------------------------

    // Coefficients of the Chebyshev series interpolating 1/(2 Pi) Sin(2 Pi K x)
    long coefficients_len;
    double *coefficients;
    double k = 0;

    if (secretKeyDist == SPARSE_TERNARY) {//if (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) {
        //        coefficients = g_coefficientsSparse;
        coefficients_len=135;
        coefficients= new double [coefficients_len];
        for (int i = 0; i < coefficients_len; ++i) {
            coefficients[i]=g_coefficientsSparse[i];
        }
        // k = K_SPARSE;
        k = 1.0;  // do not divide by k as we already did it during precomputation
    }
    else {
        //        coefficients = g_coefficientsUniform;
        coefficients_len = 89;
        coefficients = new double[coefficients_len];
        for (int i = 0; i < coefficients_len; ++i) {
            coefficients[i] = g_coefficientsUniform[i];
        }
        k = K_UNIFORM;
    }

    double constantEvalMult = pre * (1.0 / (k * N));

    scheme.multByConstAndEqual(raised, constantEvalMult);//cc->EvalMultInPlace(raised, constantEvalMult);


    // no linear transformations are needed for Chebyshev series as the range has been normalized to [-1,1]
    double coeffLowerBound = -1;
    double coeffUpperBound = 1;

    Ciphertext ctxtDec;//Ciphertext<DCRTPoly> ctxtDec;

    //    bool isLTBootstrap = (precom->m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
    //                         (precom->m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);
    bool isLTBootstrap = (precom.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) &&
                         (precom.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);


    if (slots == M / 4) {
        //------------------------------------------------------------------------------
        // FULLY PACKED CASE
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        // Running CoeffToSlot
        //------------------------------------------------------------------------------

        // need to call internal modular reduction so it also works for FLEXIBLEAUTO
        ModReduceInternalInPlace(raised,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(raised, BASE_NUM_LEVELS_TO_DROP);

        // only one linear transform is needed as the other one can be derived
//        auto ctxtEnc = (isLTBootstrap) ? EvalLinearTransform(precom->m_U0hatTPre, raised) :
//                       EvalCoeffsToSlots(precom->m_U0hatTPreFFT, raised);
        Ciphertext ctxtEnc;
        uint64_t *ctxtEnc_ax;
        uint64_t *ctxtEnc_bx;
        if (isLTBootstrap) {
            long ctxtEnc_limbs = raised.curr_limbs;
            ctxtEnc_ax = new uint64_t[ctxtEnc_limbs << logN];
            ctxtEnc_bx = new uint64_t[ctxtEnc_limbs << logN];
            ctxtEnc = Ciphertext(ctxtEnc_ax, ctxtEnc_bx, raised.N, raised.slots, ctxtEnc_limbs);
            EvalLinearTransform(ctxtEnc, precom.m_U0hatTPre, LTMatrix_Row, raised, scheme);
        } else {
//            long ctxtEnc_limbs = raised.curr_limbs-precom.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET]+1;
//            ctxtEnc_ax = new uint64_t[ctxtEnc_limbs << logN];
//            ctxtEnc_bx = new uint64_t[ctxtEnc_limbs << logN];
//            ctxtEnc = Ciphertext(ctxtEnc_ax, ctxtEnc_bx, raised.N, raised.slots, ctxtEnc_limbs);
            ctxtEnc=raised;
            EvalCoeffsToSlots(ctxtEnc, precom.m_U0hatTPreFFT, raised.slots, raised, scheme);
        }

//        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ctxtEnc->GetKeyTag());
//        auto conj       = Conjugate(ctxtEnc, evalKeyMap);
        Ciphertext conj = ctxtEnc;
        Conjugate_KeyGen(scheme.secretKey, scheme);//todo: 将密钥生成挪到预计算中去
        Conjugate_demo(conj, scheme);

        Ciphertext ctxtEncI=scheme.sub(ctxtEnc,conj); //auto ctxtEncI   = cc->EvalSub(ctxtEnc, conj);
        scheme.addAndEqual(ctxtEnc,conj);   //cc->EvalAddInPlace(ctxtEnc, conj);
        scheme.multByMonomialAndEqual(ctxtEncI,3*M/4); //algo->MultByMonomialInPlace(ctxtEncI, 3 * M / 4);

        if (rescaleTech == FIXEDMANUAL) {
//            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
//                cc->ModReduceInPlace(ctxtEnc);
//                cc->ModReduceInPlace(ctxtEncI);
            //根据单步调试知此处只调用一次 rescale
            ModReduceInternalInPlace(ctxtEnc,BASE_NUM_LEVELS_TO_DROP,scheme);
            ModReduceInternalInPlace(ctxtEncI,BASE_NUM_LEVELS_TO_DROP,scheme);
//            }
        }
        else {
//            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
//                algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
//                algo->ModReduceInternalInPlace(ctxtEncI, BASE_NUM_LEVELS_TO_DROP);
//            }
        }

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------

        // Evaluate Chebyshev series for the sine wave
//        ctxtEnc  = cc->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);
//        ctxtEncI = cc->EvalChebyshevSeries(ctxtEncI, coefficients, coeffLowerBound, coeffUpperBound);
        Ciphertext ctxtEnc_copy=ctxtEnc;
        Ciphertext ctxtEncI_copy=ctxtEncI;
        EvalChebyshevSeries(ctxtEnc,ctxtEnc_copy,coefficients,coeffLowerBound,coeffUpperBound,coefficients_len,scheme);
        EvalChebyshevSeries(ctxtEncI,ctxtEncI_copy,coefficients,coeffLowerBound,coeffUpperBound,coefficients_len,scheme);

        // Double-angle iterations are applied in the case of OPTIMIZED/uniform secrets
        if (secretKeyDist == UNIFORM_TERNARY) {//if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) {
            if (rescaleTech != FIXEDMANUAL) {//if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                ModReduceInternalInPlace(ctxtEnc,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                ModReduceInternalInPlace(ctxtEncI,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(ctxtEncI, BASE_NUM_LEVELS_TO_DROP);
            }
            ApplyDoubleAngleIterations(ctxtEnc, scheme);//ApplyDoubleAngleIterations(ctxtEnc);
            ApplyDoubleAngleIterations(ctxtEncI, scheme);//ApplyDoubleAngleIterations(ctxtEncI);
        }

        scheme.multByMonomialAndEqual(ctxtEncI,M/4); //algo->MultByMonomialInPlace(ctxtEncI, M / 4);
        scheme.addAndEqual(ctxtEnc,ctxtEncI); //cc->EvalAddInPlace(ctxtEnc, ctxtEncI);

        // scale the message back up after Chebyshev interpolation
        scheme.MultByIntegerInPlace(ctxtEnc,scalar); //algo->MultByIntegerInPlace(ctxtEnc, scalar);

        //------------------------------------------------------------------------------
        // Running SlotToCoeff
        //------------------------------------------------------------------------------

        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (rescaleTech != FIXEDMANUAL) {
            ModReduceInternalInPlace(ctxtEnc,BASE_NUM_LEVELS_TO_DROP,scheme); //algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
        }

        // Only one linear transform is needed
//        ctxtDec = (isLTBootstrap) ? EvalLinearTransform(precom->m_U0Pre, ctxtEnc) :
//                  EvalSlotsToCoeffs(precom->m_U0PreFFT, ctxtEnc);
        if (isLTBootstrap) {
            long ctxtDec_curr_limbs = ctxtEnc.curr_limbs;
            uint64_t *ax = new uint64_t[ctxtDec_curr_limbs << logN];
            uint64_t *bx = new uint64_t[ctxtDec_curr_limbs << logN];
            ctxtDec = Ciphertext(ax, bx, ctxtEnc.N, ctxtEnc.slots, ctxtDec_curr_limbs);
            EvalLinearTransform(ctxtDec, precom.m_U0Pre, LTMatrix_Row, ctxtEnc, scheme);
        } else {
            long ctxtDec_curr_limbs = ctxtEnc.curr_limbs-precom.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET]+1;
            uint64_t *ax = new uint64_t[ctxtDec_curr_limbs << logN];
            uint64_t *bx = new uint64_t[ctxtDec_curr_limbs << logN];
            ctxtDec = Ciphertext(ax, bx, ctxtEnc.N, ctxtEnc.slots, ctxtDec_curr_limbs);
            EvalSlotsToCoeffs(ctxtDec, precom.m_U0PreFFT, ctxtEnc.slots, ctxtEnc, scheme);
        }

    } else {
        //------------------------------------------------------------------------------
        // SPARSELY PACKED CASE
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // Running PartialSum
        //------------------------------------------------------------------------------
        for (int j = 1; j < N / (2 * slots); j <<= 1) {
            Ciphertext temp = raised;
            //            auto temp = cc->EvalRotate(raised, j * slots);
            FastRotate_KeyGen(scheme.secretKey, j * slots, scheme);//todo: 将密钥生成挪到预计算中去
            FastRotate_demo(temp, j * slots, scheme);

            scheme.addAndEqual(raised, temp);//cc->EvalAddInPlace(raised, temp);
        }

#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,raised,"after PartialSum");
        StringUtils::CompareAndSet(raised.curr_limbs,N,raised,
                                   reinterpret_cast<uint64_t *>(After_PartialSum[0]),
                                   "after PartialSum");
#endif


        //------------------------------------------------------------------------------
        // Running CoeffsToSlots
        //------------------------------------------------------------------------------

        ModReduceInternalInPlace(raised,BASE_NUM_LEVELS_TO_DROP,scheme);

#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,raised,"before LT");
        StringUtils::CompareAndSet(raised.curr_limbs,N,raised,
                                   reinterpret_cast<uint64_t *>(Before_LT),
                                   "Before LT");
#endif

        //        auto ctxtEnc = (isLTBootstrap) ? EvalLinearTransform(precom->m_U0hatTPre, raised) :
        //                       EvalCoeffsToSlots(precom->m_U0hatTPreFFT, raised);
        Ciphertext ctxtEnc;
        uint64_t *ctxtEnc_ax;
        uint64_t *ctxtEnc_bx;
        if (isLTBootstrap) {
            long ctxtEnc_limbs = raised.curr_limbs;
            ctxtEnc_ax = new uint64_t[ctxtEnc_limbs << logN];
            ctxtEnc_bx = new uint64_t[ctxtEnc_limbs << logN];
            ctxtEnc = Ciphertext(ctxtEnc_ax, ctxtEnc_bx, raised.N, raised.slots, ctxtEnc_limbs);
            EvalLinearTransform(ctxtEnc, precom.m_U0hatTPre, LTMatrix_Row, raised, scheme);
        } else {
            ctxtEnc=raised;
            EvalCoeffsToSlots(ctxtEnc, precom.m_U0hatTPreFFT, raised.slots, raised, scheme);
        }

#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,ctxtEnc,"after LT");
        StringUtils::CompareAndSet(ctxtEnc.curr_limbs,N,ctxtEnc,
                                   reinterpret_cast<uint64_t *>(After_LT),
                                   "after LT");
#endif

        // auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ctxtEnc->GetKeyTag());
        // auto conj = Conjugate(ctxtEnc, evalKeyMap);
        Ciphertext conj = ctxtEnc;
        Conjugate_KeyGen(scheme.secretKey, scheme);//todo: 将密钥生成挪到预计算中去
        Conjugate_demo(conj, scheme);
        scheme.addAndEqual(ctxtEnc, conj);//cc->EvalAddInPlace(ctxtEnc, conj);

        if (rescaleTech == FIXEDMANUAL) {   //if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
            //            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
            //                cc->ModReduceInPlace(ctxtEnc);
            //            }
            //根据单步调试知此处只调用一次 rescale
            ModReduceInternalInPlace(ctxtEnc,1,scheme);
        }
        else {
            //            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
            //                algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
            //            }
        }

#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,ctxtEnc,"after_add_conj");
#endif

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------

        // Evaluate Chebyshev series for the sine wave
        //ctxtEnc = cc->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);
        Ciphertext ctxtEnc_copy=ctxtEnc;
        EvalChebyshevSeries(ctxtEnc,ctxtEnc_copy,coefficients,coeffLowerBound,coeffUpperBound,coefficients_len,scheme);

#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,ctxtEnc,"after_Chebyshev");
        StringUtils::CompareAndSet(ctxtEnc.curr_limbs,N,ctxtEnc,
                                   reinterpret_cast<uint64_t *>(afterChebyshev),
                                   "after_Chebyshev");
#endif

        // Double-angle iterations are applied in the case of OPTIMIZED/uniform secrets
        if (secretKeyDist == UNIFORM_TERNARY) {//if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) {
            if (rescaleTech != FIXEDMANUAL) {//if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
                //algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
                ModReduceInternalInPlace(ctxtEnc,BASE_NUM_LEVELS_TO_DROP,scheme);
            }
            ApplyDoubleAngleIterations(ctxtEnc, scheme);//ApplyDoubleAngleIterations(ctxtEnc);
        }

#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,ctxtEnc,"after_ApplyDoubleAngleIterations");
        StringUtils::CompareAndSet(ctxtEnc.curr_limbs,N,ctxtEnc,
                                   reinterpret_cast<uint64_t *>(afterApplyDoubleAngle),
                                   "after_ApplyDoubleAngleIterations");
#endif

        // scale the message back up after Chebyshev interpolation
        scheme.MultByIntegerInPlace(ctxtEnc, scalar);//algo->MultByIntegerInPlace(ctxtEnc, scalar);


        //------------------------------------------------------------------------------
        // Running SlotsToCoeffs
        //------------------------------------------------------------------------------

        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (rescaleTech != FIXEDMANUAL) {//if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
            ModReduceInternalInPlace(ctxtEnc,BASE_NUM_LEVELS_TO_DROP,scheme);//algo->ModReduceInternalInPlace(ctxtEnc, BASE_NUM_LEVELS_TO_DROP);
        }

        // linear transform for decoding
        //        ctxtDec = (isLTBootstrap) ? EvalLinearTransform(precom->m_U0Pre, ctxtEnc) :
        //                  EvalSlotsToCoeffs(precom->m_U0PreFFT, ctxtEnc);
        if (isLTBootstrap) {
            long ctxtDec_curr_limbs = ctxtEnc.curr_limbs;
            uint64_t *ax = new uint64_t[ctxtDec_curr_limbs << logN];
            uint64_t *bx = new uint64_t[ctxtDec_curr_limbs << logN];
            ctxtDec = Ciphertext(ax, bx, ctxtEnc.N, ctxtEnc.slots, ctxtDec_curr_limbs);
            EvalLinearTransform(ctxtDec, precom.m_U0Pre, LTMatrix_Row, ctxtEnc, scheme);
        } else {
            long ctxtDec_curr_limbs = ctxtEnc.curr_limbs-precom.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET]+1;
            uint64_t *ax = new uint64_t[ctxtDec_curr_limbs << logN];
            uint64_t *bx = new uint64_t[ctxtDec_curr_limbs << logN];
            ctxtDec = Ciphertext(ax, bx, ctxtEnc.N, ctxtEnc.slots, ctxtDec_curr_limbs);
            EvalSlotsToCoeffs(ctxtDec, precom.m_U0PreFFT, ctxtEnc.slots, ctxtEnc, scheme);
        }


#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,ctxtDec,"afterLTS2C");
#endif

        //cc->EvalAddInPlace(ctxtDec, cc->EvalRotate(ctxtDec, slots));
        Ciphertext ctxtDec_rot = ctxtDec;
        FastRotate_KeyGen(scheme.secretKey, slots, scheme);//todo: 将密钥生成挪到预计算中去
        FastRotate_demo(ctxtDec_rot, slots, scheme);
        scheme.addAndEqual(ctxtDec, ctxtDec_rot);


#ifdef DEBUG_EvalBootstrap
        StringUtils::DecAndShow(scheme,ctxtDec,"afterLTS2CAdd");
#endif

    }


#if NATIVEINT != 128
    // 64-bit only: scale back the message to its original scale.
    uint64_t corFactor = (uint64_t) 1 << std::llround(correction);
    scheme.MultByIntegerInPlace(ctxtDec, corFactor);//algo->MultByIntegerInPlace(ctxtDec, corFactor);
#endif


    //yhh added
    ModReduceInternalInPlace(ctxtDec,1,scheme);
    result = ctxtDec; //return ctxtDec;

}

