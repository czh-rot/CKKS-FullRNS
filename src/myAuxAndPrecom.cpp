//
// Created by EYx on 2023/3/15.
//
#include <stdint.h>

#include "Context.h"
#include "Scheme.h"
#include "../data/testConstValue.h"
#include "../data/BSConstValue.h"
#include "Ciphertext.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/number.hpp>
using namespace boost::multiprecision;

void SetQiRelative(long L, uint64_t *qVec, uint64_t *qRoots);
void SetPiRelative(long K, uint64_t *pVec, uint64_t *pRoots);


void SetQiRelative(long L, uint64_t *qVec, uint64_t *qRoots )
{
    uint64_t *moduliQ_ptr;
    uint64_t *rootsQ_ptr;

    if(L==18){
        static uint64_t moduliQ18[ ] = {
                1152921504606844417, 576460752303458689, 576460752303406849, 576460752303447041, 576460752303408001,
                576460752303444097, 576460752303408257, 576460752303440129, 576460752303408641, 576460752303439873,
                576460752303415297, 576460752303436801, 576460752303418369, 576460752303434881, 576460752303419393,
                576460752303434497, 576460752303421441, 576460752303430529,};
        static uint64_t rootsQ18[ ] = {
                42988700452716623,407685219860279,13528086084326023,7020526219154243,1567405839084218,
                8104022074145904,16558719127045077,9410206673896024,16573956795099779,14776802342837503,
                9314221770975215,3879911350317934,9101066315480786,8834849682462566,20244755546499600,
                1249116554351241,1635442540500812,10419305971193469,};
        moduliQ_ptr=moduliQ18;
        rootsQ_ptr=rootsQ18;
    }
    else if(L==24){
        static uint64_t moduliQ24 [  ] = {
                1152921504606844417,576460752303471617,576460752303399553,576460752303464321,576460752303401729,576460752303460097,576460752303402881,576460752303458689,576460752303406849,576460752303447041,576460752303408001,576460752303444097,576460752303408257,576460752303440129,576460752303408641,576460752303439873,576460752303415297,576460752303436801,576460752303418369,576460752303434881,576460752303419393,576460752303434497,576460752303421441,576460752303430529,};
        static uint64_t rootsQ24 [  ] = {
                42988700452716623,13676418398379796,798579006277980,3253175203072499,29556602934968443,6147396647801927,5141949894743272,407685219860279,13528086084326023,7020526219154243,1567405839084218,8104022074145904,16558719127045077,9410206673896024,16573956795099779,14776802342837503,9314221770975215,3879911350317934,9101066315480786,8834849682462566,20244755546499600,1249116554351241,1635442540500812,10419305971193469,};
        moduliQ_ptr=moduliQ24;
        rootsQ_ptr=rootsQ24;
    }
    else{
        cout<<"---------------------------------------------\n"
              "error! moduliQ and rootsQ is not precomputed!\n"
              "---------------------------------------------\n";
    }

    for (int i = 0; i < L; ++i) {
        qVec[i]=moduliQ_ptr[i];
        qRoots[i]=rootsQ_ptr[i];
    }

}

void SetPiRelative(long K, uint64_t *pVec, uint64_t *pRoots )
{
    uint64_t *moduliP_ptr;
    uint64_t *rootsP_ptr;
    if(K==18){ //L=18
        static uint64_t moduliP18[] = {
                1152921504606844289, 1152921504606842753, 1152921504606837377, 1152921504606832769, 1152921504606832001,
                1152921504606831233, 1152921504606830593, 1152921504606830209, 1152921504606829697, 1152921504606827009,
                1152921504606823681, 1152921504606823297, 1152921504606815233, 1152921504606811393, 1152921504606808193,
                1152921504606807937, 1152921504606798721, 1152921504606798337,};
        static uint64_t rootsP18[] = {
                417413482957699, 2906473766327440, 5401323409338715, 32314490235789777, 7657528677974602, 27901166007825364,
                22267544142122811, 43994107232667, 7669287731435520, 32998108013266336, 11087381597658392, 13440292497346398,
                9444566094686023, 3583187404038740, 6738027634249363, 5866424793477846, 6962540504886632, 28972929474317221,};
        moduliP_ptr=moduliP18;
        rootsP_ptr=rootsP18;
    }
    else if(K==24){ //L=24
        static uint64_t moduliP24 [  ] = {
                1152921504606844289,1152921504606842753,1152921504606837377,1152921504606832769,1152921504606832001,1152921504606831233,1152921504606830593,1152921504606830209,1152921504606829697,1152921504606827009,1152921504606823681,1152921504606823297,1152921504606815233,1152921504606811393,1152921504606808193,1152921504606807937,1152921504606798721,1152921504606798337,1152921504606796289,1152921504606791809,1152921504606791681,1152921504606790913,1152921504606790657,1152921504606789761,};
        static uint64_t rootsP24 [  ] = {
                417413482957699,2906473766327440,5401323409338715,32314490235789777,7657528677974602,27901166007825364,22267544142122811,43994107232667,7669287731435520,32998108013266336,11087381597658392,13440292497346398,9444566094686023,3583187404038740,6738027634249363,5866424793477846,6962540504886632,28972929474317221,16609672961740596,6680450479684054,23492291655152123,17922869716876657,8987812706034114,13733792052411770,};
        moduliP_ptr=moduliP24;
        rootsP_ptr=rootsP24;
    }
    else{
        cout<<"---------------------------------------------\n"
              "error! moduliP and rootsP is not precomputed!\n"
              "---------------------------------------------\n";
    }

    for (int i = 0; i < K; ++i) {
        pVec[i]=moduliP_ptr[i];
        pRoots[i]=rootsP_ptr[i];
    }
}

Ciphertext::~Ciphertext(void) {
    if (ax != nullptr)
        delete[] ax;

    if (bx != nullptr)
        delete[] bx;
}

void Context::myEvalAndEqual(uint64_t *ra, uint64_t *a, long l, int j) {
    INTTAndEqual(a, l);
    for (long i = 0; i < l; ++i) {
        uint64_t *ai = a + (i << logN);
        uint64_t *rai = ra + (i << logN);
        uint64_t factor = 0;
        mulModBarrett(factor, QjHatModqi[j][i], PModq[i], qVec[i], qrVec[i], qTwok[i]);
        for (long n = 0; n < N; ++n) {
            mulModBarrett(rai[n], ai[n], factor, qVec[i], qrVec[i], qTwok[i]);
        }
    }
    NTTAndEqual(a, l);
    NTTAndEqual(ra, l);
}

void Context::myPrecomputation()
{

    // compute Q by czh
    ModulusBigint = cpp_int(1);
    for (uint64_t i = 0; i < L; ++i)
    {
        auto qi = qVec[i];
        ModulusBigint *= qi;
    }

    //compute bigintchain
    bigintChain = new cpp_int[L];
    bigintChain[0] = cpp_int(qVec[0]);
    for (uint64_t i = 1; i < L; i++)
    {
        bigintChain[i] = cpp_int(uint64_t(qVec[i]));
        bigintChain[i] *= bigintChain[i-1];
    }

    // czh: New pre-compute for LattigoMultByConst
    nttPsi = new uint64_t*[L]();
    nttPsiInv = new uint64_t*[L]();
    psiMont = new uint64_t[L]();
    psiInvMont = new uint64_t[L]();
    auto bitLenofN = uint64_t(Len64(N) - 1);

    for (int i = 0; i < L; ++i)
    {
        nttPsi[i] = new uint64_t[N]();
        nttPsiInv[i] = new uint64_t[N]();

        auto g = findPrimitiveRoot(qVec[i]);
        auto _2n = uint64_t (N << 1);
        auto power = (qVec[i] - 1) / _2n;
        auto powerInv = (qVec[i] - 1) - power;

        uint64_t bredparams[] = {0, 0};
        ComputeBRedParameters(qVec[i], bredparams[0], bredparams[1]);
        auto res = ModExp(g, power, qVec[i]);
        auto resInv = ModExp(g, powerInv, qVec[i]);
        auto PsiMont = MForm(res, qVec[i], bredparams);
        auto PsiInvMont = MForm(resInv, qVec[i], bredparams);
        psiMont[i] = PsiMont;
        psiInvMont[i] = PsiInvMont;

        nttPsi[i][0] = MForm(1, qVec[i], bredparams);
        nttPsiInv[i][0] = MForm(1, qVec[i], bredparams);

        for (long j = 1; j < N; ++j)
        {
            auto indexReversePrev = BitReverse64(j-1, bitLenofN);
            auto indexReverseNext = BitReverse64(j, bitLenofN);

            nttPsi[i][indexReverseNext] = MRed(nttPsi[i][indexReversePrev], PsiMont, qVec[i], qInvVec[i]);
            nttPsiInv[i][indexReverseNext] = MRed(nttPsiInv[i][indexReversePrev], PsiInvMont, qVec[i], qInvVec[i]);
        }
    }

    //KeySwitch params
//    alpha=(L+1)/dnum;
    alpha=L/dnum;
    //KSGen params
    QjHatModqi= new uint64_t*[dnum];
    for (long j = 0; j < dnum; ++j) {
        QjHatModqi[j]= new uint64_t [L]();
        long QjLowerBound=j*alpha;
        long QjUpperBound=(j+1)*alpha;
        for (long i = 0; i < L ; ++i) {
            QjHatModqi[j][i]=1;
            for (long k = 0; k < QjLowerBound; ++k) {
                uint64_t temp = qVec[k] % qVec[i];
                mulMod(QjHatModqi[j][i], QjHatModqi[j][i], temp, qVec[i]);
            }
            for (long k = QjUpperBound; k < L; ++k) {
                uint64_t temp = qVec[k] % qVec[i];
                mulMod(QjHatModqi[j][i], QjHatModqi[j][i], temp, qVec[i]);
            }
//            cout<<"QjHatModqi["<<j<<"]["<<i<<"]"<<QjHatModqi[j][i]<<endl;
        }
    }

    //cout
//    for (int i = 0; i < L; ++i) {
//        cout<<"qVec["<<i<<"]"<<qVec[i]<<endl;
//    }

    QjHatInvModqi= new uint64_t * [dnum];
    for (long j = 0; j<dnum ; ++j) {
        QjHatInvModqi[j]= new uint64_t [alpha]();
        for (long i = 0; i < alpha ; ++i) {
            int index=j*alpha+i;
            QjHatInvModqi[j][i]=QjHatModqi[j][index];
            QjHatInvModqi[j][i] = invMod(QjHatInvModqi[j][i], qVec[index]);
        }
    }

//    //Basis conversion params
//    QjHatDprimeModqi= new uint64_t*[L];
//    for (long j = 0; j < L; ++j) {
//        QjHatDprimeModqi[j]= new uint64_t []();
//        int tmp_beta=j/alpha;
//        long QjLowerBound=tmp_beta*alpha;
//        long QjUpperBound=(tmp_beta+1)*alpha;
//        for (long k = 0, index = QjLowerBound ; k < alpha ; ++k,++index) {
//            QjHatDprimeModqi[j][k] = PModq[index];
//            for (int i=QjLowerBound; i <QjUpperBound ; ++i) {
//                    uint64_t temp = qVec[i]%qVec[index];
//                    mulMod(QjHatDprimeModqi[j][k], QjHatDprimeModqi[j][k], temp, qVec[index]);
//            }
//        }
//        for (long k = 0; k < QjLowerBound; ++k) {
//            uint64_t temp = qVec[k] % qVec[i];
//            mulMod(QjHatDprimeModqi[j][i], QjHatDprimeModqi[j][i], temp, qVec[i]);
//        }
//        for (long k = QjUpperBound; k < L; ++k) {
//            uint64_t temp = qVec[k] % qVec[i];
//            mulMod(QjHatDprimeModqi[j][i], QjHatDprimeModqi[j][i], temp, qVec[i]);
//        }
//
//    }

    //rescale param
    // sizeQ in openFHE equals to L here.
    QlQlInvModqlDivqlModq=new uint64_t* [L-1];
    for (int k = 0; k < L - 1; ++k) {
        int l=L-(k+1);
        QlQlInvModqlDivqlModq[k]= new uint64_t [l];
        for (int i = 0; i < l; ++i) {
            uint64_t QlInvModql=1;
            for (int j = 0; j < l; ++j) {
                uint64_t temp= invMod(qVec[j],qVec[l]);
                mulMod(QlInvModql,QlInvModql,temp,qVec[l]);
            }

            cpp_int modulusQ=1;
            for (int j = 0; j < l; ++j)
                modulusQ *=qVec[j];

            cpp_int result=(QlInvModql*modulusQ)/qVec[l];
            result%=qVec[i];

            QlQlInvModqlDivqlModq[k][i]=result.convert_to<uint64_t>();
        }
    }
}

Context::Context(string STRING, long logN, long logp, long L, long K,  long dnum, long h, double sigma) :
        logN(logN), logp(logp), L(L), K(K), dnum(dnum), h(h), sigma(sigma) {

    N = 1L << logN;
    M = N << 1;
    logNh = logN - 1;
    Nh = N >> 1;
    p = 1L << logp;


    qVec = new uint64_t[L]();
    qrVec = new uint64_t[L]();
    qTwok = new long[L]();
    qkVec = new uint64_t[L]();
    qdVec = new uint64_t[L]();
    qInvVec = new uint64_t[L]();
    qRoots = new uint64_t[L]();
    qRootsInv = new uint64_t[L]();
    qRootPows = new uint64_t*[L];
    qRootScalePows = new uint64_t*[L];
    qRootScalePowsOverq = new uint64_t*[L];
    qRootScalePowsInv = new uint64_t*[L];
    qRootPowsInv = new uint64_t*[L];
    NInvModq = new uint64_t[L]();
    NScaleInvModq = new uint64_t[L]();

    // Generate Primes //
    long bnd = 1;
    long cnt = 1;

    bnd = 1;
    while(1) {
        uint64_t prime = (1ULL << Q0_BIT_SIZE) + bnd * M + 1;
        if(primeTest(prime)) {
            qVec[0] = prime;
            break;
        }
        bnd++;
    }

    bnd = 1;
    while(cnt < L) {
        uint64_t prime1 = (1ULL << logp) + bnd * M + 1;
        if(primeTest(prime1)) {
            qVec[cnt] = prime1;
            cnt++;
        }
        uint64_t prime2 = (1ULL << logp) - bnd * M + 1;
        if(primeTest(prime2)) {
            qVec[cnt] = prime2;
            cnt++;
        }
        bnd++;
    }

    if(logp - logN - 1 - ceil(log2(bnd)) < 10) {
        cerr << "ERROR: too small number of precision" << endl;
        cerr << "TRY to use larger logp or smaller depth" << endl;
    }

    //Debug qVec[i]=moduliQ[i];
    SetQiRelative(L,qVec,qRoots);

    for (long i = 0; i < L; ++i) {
        qTwok[i] = (2 * ((long)log2(qVec[i]) + 1)); //ModBarrett 算法里的 k，详见笔记
        qrVec[i] = (static_cast<unsigned __int128>(1) << qTwok[i]) / qVec[i];   //ModBarret算法里的 m，详见笔记
        qkVec[i] = static_cast<uint64_t>(((static_cast<unsigned __int128>(invMod(((uint64_t)(1) << 62), qVec[i])) << 62) - 1) / qVec[i]);
        qdVec[i] = qVec[i] << 1;
//        qRoots[i] = findMthRootOfUnity(M, qVec[i]);
        qRootsInv[i] = invMod(qRoots[i], qVec[i]);
        NInvModq[i] = invMod(N, qVec[i]);
        mulMod(NScaleInvModq[i], NInvModq[i], (static_cast<uint64_t>(1) << 32), qVec[i]);
        mulMod(NScaleInvModq[i], NScaleInvModq[i], (static_cast<uint64_t>(1) << 32), qVec[i]);
        qInvVec[i] = inv(qVec[i]);
        qRootPows[i] = new uint64_t[N]();
        qRootPowsInv[i] = new uint64_t[N]();
        qRootScalePows[i] = new uint64_t[N]();
        qRootScalePowsOverq[i] = new uint64_t[N]();
        qRootScalePowsInv[i] = new uint64_t[N]();
        uint64_t power = static_cast<uint64_t>(1);
        uint64_t powerInv = static_cast<uint64_t>(1);

        for (long j = 0; j < N; ++j) {
            uint64_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN);
            qRootPows[i][jprime] = power;
            unsigned __int128 tmp = (static_cast<unsigned __int128>(power) << 64);
            qRootScalePowsOverq[i][jprime] = static_cast<uint64_t>(tmp / qVec[i]);
            mulMod(qRootScalePows[i][jprime], qRootPows[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);
            mulMod(qRootScalePows[i][jprime], qRootScalePows[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);
            qRootPowsInv[i][jprime] = powerInv;
            mulMod(qRootScalePowsInv[i][jprime], qRootPowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);
            mulMod(qRootScalePowsInv[i][jprime], qRootScalePowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);

            if (j < N - 1) {
                mulMod(power, power, qRoots[i], qVec[i]);
                mulMod(powerInv, powerInv, qRootsInv[i], qVec[i]);
            }
        }
    }

    pVec = new uint64_t[K]();
    prVec = new uint64_t[K]();
    pTwok = new long[K]();
    pkVec = new uint64_t[K]();
    pdVec = new uint64_t[K]();
    pInvVec = new uint64_t[K]();
    pRoots = new uint64_t[K]();
    pRootsInv = new uint64_t[K]();
    pRootPows = new uint64_t*[K];
    pRootPowsInv = new uint64_t*[K];
    pRootScalePows = new uint64_t*[K];
    pRootScalePowsOverp = new uint64_t*[K];
    pRootScalePowsInv = new uint64_t*[K];
    NInvModp = new uint64_t[K]();
    NScaleInvModp = new uint64_t[K]();

    // Generate Special Primes //
    cnt = 0;
    while(cnt < K) {
        uint64_t prime1 = (1ULL << logp) + bnd * M + 1;
        if(primeTest(prime1)) {
            pVec[cnt] = prime1;
            cnt++;
        }
        if(cnt == K) break;
        uint64_t prime2 = (1ULL << logp) - bnd * M + 1;
        if(primeTest(prime2)) {
            pVec[cnt] = prime2;
            cnt++;
        }
        bnd++;
    }

    //Debug pVec[i]=moduliP[i];
    SetPiRelative(L, pVec, pRoots);

    for (long i = 0; i < K; ++i) {
        pTwok[i] = (2 * ((long)log2(pVec[i]) + 1));
        prVec[i] = (static_cast<unsigned __int128>(1) << pTwok[i]) / pVec[i];
        pkVec[i] = static_cast<uint64_t>(((static_cast<unsigned __int128>(invMod(((uint64_t)(1) << 62), pVec[i])) << 62) - 1) / pVec[i]);
        pdVec[i] = pVec[i] << 1;
//        pRoots[i] = findMthRootOfUnity(M, pVec[i]);
        pRootsInv[i] = invMod(pRoots[i], pVec[i]);
        NInvModp[i] = invMod(N, pVec[i]);
        mulMod(NScaleInvModp[i], NInvModp[i], (static_cast<uint64_t>(1) << 32), pVec[i]);
        mulMod(NScaleInvModp[i], NScaleInvModp[i], (static_cast<uint64_t>(1) << 32), pVec[i]);
        pRootPows[i] = new uint64_t[N]();
        pRootScalePows[i] = new uint64_t[N]();
        pRootScalePowsOverp[i] = new uint64_t[N]();
        pRootScalePowsInv[i] = new uint64_t[N]();
        pRootPowsInv[i] = new uint64_t[N]();
        pInvVec[i] = inv(pVec[i]);
        uint64_t power = static_cast<uint64_t>(1);
        uint64_t powerInv = static_cast<uint64_t>(1);
        for (long j = 0; j < N; ++j) {
            uint64_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN);
            pRootPows[i][jprime] = power;
            unsigned __int128 tmp = (static_cast<unsigned __int128>(power) << 64);
            mulMod(pRootScalePows[i][jprime], pRootPows[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
            mulMod(pRootScalePows[i][jprime], pRootScalePows[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
            pRootPowsInv[i][jprime] = powerInv;
            mulMod(pRootScalePowsInv[i][jprime], pRootPowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
            mulMod(pRootScalePowsInv[i][jprime], pRootScalePowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
            if (j < N - 1) {
                mulMod(power, power, pRoots[i], pVec[i]);
                mulMod(powerInv, powerInv, pRootsInv[i], pVec[i]);
            }
        }
    }

    qHatModq = new uint64_t*[L]; // [curr_limbs][i] (phat_i)_l mod p_i
    for (long l = 0; l < L; ++l) {
        qHatModq[l] = new uint64_t[l + 1]();
        for (long i = 0; i < l + 1; ++i) {
            qHatModq[l][i] = 1;
            for (long j = 0; j < i; ++j) {
                uint64_t temp = qVec[j] % qVec[i];
                mulMod(qHatModq[l][i], qHatModq[l][i], temp, qVec[i]);
            }
            for (long j = i + 1; j < l + 1; ++j) {
                uint64_t temp = qVec[j] % qVec[i];
                mulMod(qHatModq[l][i], qHatModq[l][i], temp, qVec[i]);
            }
        }
    }

    qHatInvModq = new uint64_t*[L]; // [curr_limbs][i] (phat_i)_l^-1 mod p_i
    for (long l = 0; l < L; ++l) {
        qHatInvModq[l] = new uint64_t[l + 1]();
        for (long i = 0; i < l + 1; ++i) {
            qHatInvModq[l][i] = invMod(qHatModq[l][i], qVec[i]);
        }
    }

    pHatModp = new uint64_t[K](); // [k] qhat_k mod q_k

    for (long k = 0; k < K; ++k) {
        pHatModp[k] = 1;
        for (long j = 0; j < k; ++j) {
            uint64_t temp = pVec[j] % pVec[k];
            mulMod(pHatModp[k], pHatModp[k], temp, pVec[k]);
        }
        for (long j = k + 1; j < K; ++j) {
            uint64_t temp = pVec[j] % pVec[k];
            mulMod(pHatModp[k], pHatModp[k], temp, pVec[k]);
        }
    }

    pHatInvModp = new uint64_t[K](); // [k] qhat_k^-1 mod q_k
    for (long k = 0; k < K; ++k) {
        pHatInvModp[k] = invMod(pHatModp[k], pVec[k]);
    }

    qHatModp = new uint64_t**[L];  // [curr_limbs] [i] [k]  (phat_i)_l mod q_k
    for (long l = 0; l < L; ++l) {
        qHatModp[l] = new uint64_t*[l + 1];
        for (long i = 0; i < l + 1; ++i) {
            qHatModp[l][i] = new uint64_t[K]();
            for (long k = 0; k < K; ++k) {
                qHatModp[l][i][k] = 1;
                for (long j = 0; j < i; ++j) {
                    uint64_t temp = qVec[j] % pVec[k];
                    mulMod(qHatModp[l][i][k], qHatModp[l][i][k], temp, pVec[k]);
                }
                for (long j = i + 1; j < l + 1; ++j) {
                    uint64_t temp = qVec[j] % pVec[k];
                    mulMod(qHatModp[l][i][k], qHatModp[l][i][k], temp, pVec[k]);
                }
            }
        }
    }

    pHatModq = new uint64_t*[K]; // [k][i] qhat_k mod p_i
    for (long k = 0; k < K; ++k) {
        pHatModq[k] = new uint64_t[L]();
        for (long i = 0; i < L; ++i) {
            pHatModq[k][i] = 1;
            for (long s = 0; s < k; ++s) {
                uint64_t temp = pVec[s] % qVec[i];
                mulMod(pHatModq[k][i], pHatModq[k][i], temp, qVec[i]);
            }
            for (long s = k + 1; s < K; ++s) {
                uint64_t temp = pVec[s] % qVec[i];
                mulMod(pHatModq[k][i], pHatModq[k][i], temp, qVec[i]);
            }
        }
    }



    PModq = new uint64_t[L](); // [i] qprod mod p_i
    for (long i = 0; i < L; ++i) {
        PModq[i] = 1;
        for (long k = 0; k < K; ++k) {
            uint64_t temp = pVec[k] % qVec[i];
            mulMod(PModq[i], PModq[i], temp, qVec[i]);
        }
    }

    PInvModq = new uint64_t[L](); // [i] qprod^-1 mod p_i
    for (long i = 0; i < L; ++i) {
        PInvModq[i] = invMod(PModq[i], qVec[i]);
    }

    QModp = new uint64_t*[L];
    for (long i = 0; i < L; ++i) {
        QModp[i] = new uint64_t[K]();
        for (long k = 0; k < K; ++k) {
            QModp[i][k] = 1;
            for (long j = 0; j < i + 1; ++j) {
                uint64_t temp = qVec[j] % pVec[k];
                mulMod(QModp[i][k], QModp[i][k], temp, pVec[k]);
            }
        }
    }
    QInvModp = new uint64_t*[L];
    for (long i = 0; i < L; ++i) {
        QInvModp[i] = new uint64_t[K]();
        for (long k = 0; k < K; ++k) {
            QInvModp[i][k] = invMod(QModp[i][k], pVec[k]);
        }
    }

    qInvModq = new uint64_t*[L]; // [i][j] p_i^-1 mod p_j
    for (long i = 0; i < L; ++i) {
        qInvModq[i] = new uint64_t[L]();
        for (long j = 0; j < i; ++j) {
            qInvModq[i][j] = invMod(qVec[i], qVec[j]);
        }
        for (long j = i + 1; j < L; ++j) {
            qInvModq[i][j] = invMod(qVec[i], qVec[j]);
        }
    }

    rotGroup = new long[Nh]();
    long fivePows = 1;
    for (long i = 0; i < Nh; ++i) {
        rotGroup[i] = fivePows;
        fivePows *= 5;
        fivePows %= M;
    }

    ksiPows = new complex<double>[M + 1];
    for (long j = 0; j < M; ++j) {
        double angle = 2.0 * M_PI * j / M;
        ksiPows[j].real(cos(angle));
        ksiPows[j].imag(sin(angle));
    }

    ksiPows[M] = ksiPows[0];

    p2coeff = new uint64_t[L << logN];

    for (long i = 0; i < L; ++i) {
        for (long n = 0; n < N; ++n) {
            mulModBarrett(p2coeff[n + (i << logN)], p, p, qVec[i], qrVec[i], qTwok[i]);
        }
    }
    p2hcoeff = new uint64_t[L << logN];

    for (long i = 0; i < L; ++i) {
        for (long n = 0; n < N; ++n) {
            mulModBarrett(p2hcoeff[n + (i << logN)], (p >> 1), p, qVec[i], qrVec[i], qTwok[i]);
        }
    }

    pccoeff = new uint64_t[L << logN]();
    for (long i = 0; i < L; ++i) {
        for (long n = 0; n < N; ++n) {
            mulModBarrett(pccoeff[n + (i << logN)], (94.2372881) * (p >> 20), (1L << 20), qVec[i], qrVec[i], qTwok[i]);
        }
    }
    negateAndEqual(pccoeff, L);


    taylorCoeffsMap.insert(pair<string, double*>(LOGARITHM, new double[11]{0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10}));
    taylorCoeffsMap.insert(pair<string, double*>(EXPONENT, new double[11]{1,1,0.5,1./6,1./24,1./120,1./720,1./5040, 1./40320,1./362880,1./3628800}));
    taylorCoeffsMap.insert(pair<string, double*>(SIGMOID, new double[11]{1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0}));

    myPrecomputation();

}

void Scheme::EvalMultExtInPlace(Ciphertext& cipher, Plaintext pt) {
    if (cipher.special_limbs != context.K) {
        throw invalid_argument("Ciphertexts are not extended");
    }

    long K = cipher.special_limbs;
    long curr_limbs=cipher.curr_limbs;
    context.mulAndEqual(cipher.ax, pt.mx, curr_limbs,K);
    context.mulAndEqual(cipher.bx, pt.mx, curr_limbs,K);
}

void Scheme::EvalMultExt(Ciphertext &result, const Ciphertext &cipher, Plaintext pt) {
    if (cipher.special_limbs != context.K
    || result.special_limbs != context.K) {
        throw invalid_argument("Ciphertexts are not extended");
    }

    long K = cipher.special_limbs;
    long curr_limbs = cipher.curr_limbs;
    context.mul(result.ax, cipher.ax, pt.mx, curr_limbs, K);
    context.mul(result.bx, cipher.bx, pt.mx, curr_limbs, K);
}

void Scheme::EvalAddExtInPlace(Ciphertext &cipher1, const Ciphertext &cipher2) {
    if (cipher1.curr_limbs != cipher2.curr_limbs) {
        throw invalid_argument("Ciphertexts are not in same level");
    }
    else if (cipher1.special_limbs != context.K || cipher2.special_limbs != context.K) {
        throw invalid_argument("Ciphertexts are not extended");
    }

    long K = cipher1.special_limbs;
    long curr_limbs = cipher1.curr_limbs;

    context.addAndEqual(cipher1.ax, cipher2.ax, curr_limbs, K);
    context.addAndEqual(cipher1.bx, cipher2.bx, curr_limbs, K);
}

//todo: 写测试函数验证正确性。对应：ckksrns-leveledshe.cpp：void LeveledSHECKKSRNS::MultByIntegerInPlace(Ciphertext<DCRTPoly>& ciphertext, uint64_t integer)
void Scheme::MultByIntegerInPlace(Ciphertext &cipher, uint64_t integer) {
    long curr_limbs = cipher.curr_limbs;
    long N = cipher.N;
    long logN = context.logN;
    uint64_t *qVec = context.qVec;
    uint64_t *qrVec = context.qrVec;
    long *qTwok = context.qTwok;

    //deal with cipher.ax
    for (long i = 0; i < curr_limbs; ++i) {
        uint64_t *cipherj = cipher.ax + (i << logN);
        for (long j = 0; j < N; ++j) {
            uint64_t tmp = 0;
            mulModBarrett(tmp, cipherj[j], integer, qVec[i], qrVec[i], qTwok[i]);
            cipherj[j] = tmp;
        }
    }

    //deal with cipher.bx
    for (long i = 0; i < curr_limbs; ++i) {
        uint64_t *cipherj = cipher.bx + (i << logN);
        for (long j = 0; j < N; ++j) {
            uint64_t tmp = 0;
            mulModBarrett(tmp, cipherj[j], integer, qVec[i], qrVec[i], qTwok[i]);
            cipherj[j] = tmp;
        }
    }
}