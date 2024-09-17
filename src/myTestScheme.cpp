//
// Created by EYx on 2023/3/20.
//
#include "TestScheme.h"
#include "Numb.h"
#include "Context.h"
#include "SecretKey.h"
#include "Scheme.h"
#include "EvaluatorUtils.h"
#include "StringUtils.h"
#include "TimeUtils.h"
#include "SchemeAlgo.h"

#include "Common.h"

#include <iostream>
#include <vector>
#include <chrono>

#include "../data/testConstValue.h"
#include "myChebyshevFunc.h"
#include "myRotation.h"
#include "myKeySwitch.h"
#include "myMult.h"
#include "ckksrns-leveledshe.h"

using namespace std;
using namespace chrono;

//Add tests
void TestScheme::testConstAdd(long logN, long L, long logp, long logSlots){
    cout << "!!! START ConstAdd !!!" << endl;
    //-----------------------------------------
    long K = L + 1;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = 2;

    complex<double> mvec1[2]={0.25,0.5};
    complex<double> mvec2[2]={5,4};
    double cvec[2]={2,4};
    uint32_t ciphertexts_num=2;
    complex<double> *golden_value = new complex<double> [ciphertexts_num];
//    golden_value[0]=mvec1[0]*cvec[0]+mvec2[0]*cvec[1];
//    golden_value[1]=mvec1[1]*cvec[0]+mvec2[1]*cvec[1];
    golden_value[0]=mvec1[0].real()+cvec[0];
    golden_value[1]=mvec1[1].real()+cvec[0];
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);

    scheme.addConstAndEqual(cipher1, cvec[0]);
    complex<double>* dvecCAdd = scheme.decrypt(secretKey, cipher1);
    StringUtils::showcompare(golden_value, dvecCAdd, slots, "ConstMult");
}

//mult tests
void TestScheme::testmySquare(long logN, long L, long logp, long logSlots, long dnum) {
    cout << "!!! START TEST mySquare !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvecMult = new complex<double>[slots];

    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i] * mvec1[i];
    }

    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    timeutils.start("Homomorphic Multiplication (Square)");
    Ciphertext multCipher(cipher1.N,cipher1.slots,cipher1.curr_limbs);
    mySquare(multCipher, cipher1, scheme);
    timeutils.stop("Homomorphic Multiplication (Square)");
    scheme.reScaleByAndEqual(multCipher, 1);
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);
    StringUtils::showcompare(mvecMult, dvecMult, slots, "mySquare");


    timeutils.start("Homomorphic Multiplication (SquareAndEqual)");
    mySquareAndEqual( cipher1, scheme);
    timeutils.stop("Homomorphic Multiplication (SquareAndEqual)");
    scheme.reScaleByAndEqual(cipher1, 1);
    dvecMult = scheme.decrypt(secretKey, cipher1);
    StringUtils::showcompare(mvecMult, dvecMult, slots, "mySquareAndEqual");
}

void TestScheme::testmyMultAndEqual(long logN, long L, long logp, long logSlots, long dnum) {
    cout << "!!! START TEST myMultAndEqual !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);

    complex<double>* mvecMult = new complex<double>[slots];

    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i] * mvec2[i];
    }

    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    timeutils.start("Homomorphic Multiplication");
    myMultAndEqual(cipher1, cipher2,scheme);
    timeutils.stop("Homomorphic Multiplication");
    timeutils.start("Rescaling");
    scheme.reScaleByAndEqual(cipher1, 1);
    timeutils.stop("Rescaling");

    timeutils.start("Decrypt batch");
    complex<double>* dvecMult = scheme.decrypt(secretKey, cipher1);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvecMult, dvecMult, slots, "testmyMultAndEqual");
}

void TestScheme::testComplexConstMult(long logN, long L, long logp, long logSlots) {
    cout << "!!! START multByComplexConst !!!" << endl;
    //-----------------------------------------
    long K = L ;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------

    complex<double> cnst(-1,-2.1);
    const long slots=1;
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
//    complex<double> mvec1[] = {{-0.86345,-0.9497},};//测试数据，对任意绝对值为 (Q0_BIT_SIZE/(1<<logp))/bound 范围内的 cnst值可以通过测试。
//    complex<double> mvec1[] = {{-3.86345,-3.9497},};//fixme: bound>1时，若slot数据为负数有问题，例如这个数据

    complex<double>* golden_value=new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        golden_value[i]=mvec1[i]*cnst;
    }
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cmultCipher = scheme.multByConst(cipher1,cnst); // 等价于，mvec1 每个元素乘以 cnst
    scheme.reScaleByAndEqual(cmultCipher, 1);
    complex<double> *dvecCMult = scheme.decrypt(secretKey, cmultCipher);
    StringUtils::showcompare(golden_value, dvecCMult, slots, "multByComplexConst",1);
}

void TestScheme::testComplexConstVecMult(long logN, long L, long logp, long logSlots) {
    //fixme: 这个测试函数的测试对象的问题太多了，还没分析清楚
    cout << "!!! START multByComplexConstVec!!!" << endl;
    //-----------------------------------------
    long K = L ;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------

    const long slots=1;
    complex<double> ComplexCnstVec[]={{-2,-2},{-1, -1},{-0.1,-0.1},  };

    double bound = 1.0;
//    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
//    complex<double> mvec1[] = {{-3.86345,-3.9497},};//fixme: bound>1时，若slot数据为负数有问题，例如这个数据
    complex<double> mvec1[] = {{3.86345,3.9497},};//fixme: bound>1时，若slot数据为正数也有问题，例如这个数据
//    complex<double> mvec1[] = {{0.86345,0.9497},};//fixme: {0.86345,0.9497}*{-2,-2} 有问题

    complex<double>* complex_golden_value=new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        complex_golden_value[i]= mvec1[i] * ComplexCnstVec[i];
    }

    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cvecmultCipher = scheme.multByConstVec(cipher1, ComplexCnstVec, slots);// 等价于，mvec2 点乘 ComplexCnstVec
    scheme.reScaleByAndEqual(cvecmultCipher, 1);
    complex<double>* dvecCVecMult = scheme.decrypt(secretKey, cvecmultCipher);
    StringUtils::showcompare(complex_golden_value, dvecCVecMult, slots, "multByComplexConstVec",1);
}

void TestScheme::testmyMult(long logN, long L, long logp, long logSlots, long dnum) {
    cout << "!!! START TEST HMult !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);

    complex<double>* mvecMult = new complex<double>[slots];

    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i] * mvec2[i];
    }

    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    timeutils.start("Homomorphic Multiplication");
    Ciphertext multCipher(cipher1.N, cipher1.slots, cipher1.curr_limbs);
    myMult(multCipher,cipher1, cipher2, scheme);
    timeutils.stop("Homomorphic Multiplication");
    timeutils.start("Rescaling");
    scheme.reScaleByAndEqual(multCipher, 1);
    timeutils.stop("Rescaling");

    timeutils.start("Decrypt batch");
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvecMult, dvecMult, slots, "HMult1");
}

void TestScheme::testmyMult2(long logN, long L, long logp, long logSlots, long dnum) {
    cout << "!!! START TEST HMult !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);

    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;

    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);

    //Debug Only
    complex<double> mvec1[2]={0.25,0.5};
    complex<double> mvec2[2]={5,4};
    complex<double>* mvecMult = new complex<double>[slots];
    for(long i = 0; i < slots; i++) {
        mvecMult[i] = (mvec1[i] * mvec2[i]) * (mvec1[i] * mvec2[i]);
    }

    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    timeutils.start("Homomorphic Multiplication & Rescaling");
    Ciphertext multCipher(cipher1.N, cipher1.slots, cipher1.curr_limbs);
    myMult(multCipher,cipher1, cipher2, scheme);
    scheme.reScaleByAndEqual(multCipher, 1);
    timeutils.stop("Homomorphic Multiplication & Rescaling");

    timeutils.start("Homomorphic Multiplication & Rescaling");
    Ciphertext multCipher2(multCipher.N, multCipher.slots, multCipher.curr_limbs);
    myMult(multCipher2,multCipher, multCipher, scheme);
    scheme.reScaleByAndEqual(multCipher2, 1);
    timeutils.stop("Homomorphic Multiplication & Rescaling");

    timeutils.start("Decrypt batch");
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher2);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvecMult, dvecMult, slots, "HMult2");
}

void TestScheme::testmyMult3(long logN, long L, long logp, long logSlots, long dnum) {
    cout << "!!! START TEST HMult3 !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 或者初始情况下要求dnum必须整除L，判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);

    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;

    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);

    //Debug Only
    complex<double> mvec1[2]={0.25,0.5};
    complex<double> mvec2[2]={5,4};
    complex<double>* mvecMult = new complex<double>[slots];
    for(long i = 0; i < slots; i++) {
        mvecMult[i] = (mvec1[i] * mvec2[i]);
        mvecMult[i] = mvecMult[i] * mvecMult[i];
        mvecMult[i] = mvecMult[i] * mvecMult[i];
    }

    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    timeutils.start("Homomorphic Multiplication & Rescaling");
    Ciphertext multCipher(cipher1.N, cipher1.slots, cipher1.curr_limbs);
    myMult(multCipher,cipher1, cipher2, scheme);
    scheme.reScaleByAndEqual(multCipher, 1);
    Ciphertext multCipher2(multCipher.N, multCipher.slots, multCipher.curr_limbs);
    myMult(multCipher2,multCipher, multCipher, scheme);
    scheme.reScaleByAndEqual(multCipher2, 1);
    Ciphertext multCipher3(multCipher2.N, multCipher2.slots, multCipher2.curr_limbs);
    myMult(multCipher3,multCipher2, multCipher2, scheme);
    scheme.reScaleByAndEqual(multCipher3, 1);
    timeutils.stop("Homomorphic Multiplication & Rescaling");

    timeutils.start("Decrypt batch");
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher3);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvecMult, dvecMult, slots, "HMult3");
}

void TestScheme::DebugmyMult(long logN, long L, long logp, long logSlots, long dnum) {//没什么用
    cout << "!!! START TEST myMult !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context("DEBUG_SK",logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context,"DEBUG_SK");
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    //Debug Only
    complex<double> mvec1[2]={(0.25+1i),(0.5+1i)};
    complex<double> mvec2[2]={(5.0+1i),(4.0+1i)};


    complex<double>* mvecMult = new complex<double>[slots];

    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i] * mvec2[i];
    }

    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    timeutils.start("Homomorphic Multiplication & Rescaling");
    Ciphertext multCipher (cipher1.N, cipher1.slots, cipher1.curr_limbs);
    myMult(multCipher, cipher1, cipher2, scheme);
    scheme.reScaleByAndEqual(multCipher, 1);
    timeutils.stop("Homomorphic Multiplication & Rescaling");


    timeutils.start("Decrypt batch");
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvecMult, dvecMult, slots, "cmult");
}

//mod up, mod down, rescale
void TestScheme::testModUp(long logN, long L, long logp, long logSlots) {
    cout << "!!! START TEST BASIC !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L ;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2(cipher1);
    timeutils.stop("Encrypt two batch");

    for (long i = 0; i < 128; ++i) {
        cipher2.ax[i] = A[i];
        cipher2.bx[i] = B[i];
    }

    long curr_limbs = cipher2.curr_limbs;

    //compute golden value
    uint64_t *d1 = cipher2.bx;
    uint64_t *d1_expand = new uint64_t[(curr_limbs + context.K) << context.logN];
    copy(d1, d1 + (cipher2.curr_limbs << logN), d1_expand);
    context.raiseAndEqual(d1_expand, cipher1.curr_limbs);

    //my computation
    uint64_t *d2 = cipher2.bx;
    uint64_t *d2_expand = new uint64_t[(curr_limbs + context.K) << context.logN];
    uint64_t *res_d2 = new uint64_t[curr_limbs << context.logN]();

    timeutils.start("ModUp");
    ApproxModUp(d2_expand, d2,
                0, context.alpha, 0, 0,
                0, curr_limbs, curr_limbs, context.K, scheme);
    timeutils.stop("ModUp");
    delete[] cipher2.bx;
    cipher2.bx = res_d2;

    StringUtils::showcompare(d1_expand, d2_expand, context.N * context.L, "ModUp");
}

void TestScheme::testModDown(long logN, long L, long logp, long logSlots) {
    cout << "!!! START TEST BASIC !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L ;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2(cipher1);
    timeutils.stop("Encrypt two batch");

    for (long i = 0; i < 128; ++i) {
        cipher2.ax[i] = A[i];
        cipher2.bx[i] = B[i];
    }

    long curr_limbs = cipher2.curr_limbs;
    long alpha = context.alpha;

    //compute golden value
    uint64_t *d1 = cipher2.bx;
    uint64_t *d1_expand = new uint64_t[(curr_limbs + K) << logN];
    copy(d1, d1 + (cipher2.curr_limbs << logN), d1_expand);
    context.raiseAndEqual(d1_expand, cipher1.curr_limbs);
    context.backAndEqual(d1_expand, cipher1.curr_limbs);

    //my computation
    uint64_t *d2 = cipher2.bx;
    uint64_t *d2_expand = new uint64_t[(curr_limbs + K) << logN];
    uint64_t *res_d2 = new uint64_t[curr_limbs << logN]();

    timeutils.start("ModUpAndDown ");
    ApproxModUp(d2_expand, d2,
                0, alpha, 0, 0,
                0, curr_limbs, curr_limbs, K, scheme);
    ApproxModDown(res_d2, d2_expand, curr_limbs,
                  curr_limbs, curr_limbs - curr_limbs, curr_limbs, K,
                  0, curr_limbs, 0, 0, scheme);
    timeutils.stop("ModUpAndDown");
    delete[] cipher2.bx;
    cipher2.bx = res_d2;

    StringUtils::showcompare(d1_expand, res_d2, (L << logN), "ModDown");
}

void TestScheme::testModUpAndDown(long logN, long L, long logp, long logSlots) {
    cout << "!!! START TEST BASIC !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L ;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double> *mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2(cipher1);
    timeutils.stop("Encrypt two batch");

    long curr_limbs = cipher2.curr_limbs;

    //compute golden value
    uint64_t *d1 = cipher1.bx;
    uint64_t *d1_expand = new uint64_t[(curr_limbs + K) << context.logN];
    copy(d1, d1 + (cipher1.curr_limbs << logN), d1_expand);
    context.raiseAndEqual(d1_expand, cipher1.curr_limbs);
    context.back(d1, d1_expand, cipher1.curr_limbs);

    //my computation
    uint64_t *d2 = cipher2.bx;
    uint64_t *d2_expand = new uint64_t[(curr_limbs + K) << context.logN];
    uint64_t *res_d2 = new uint64_t[curr_limbs << context.logN]();

    timeutils.start("ModUpAndDown ");
    ApproxModUp(d2_expand, d2,
                0, context.alpha, 0, 0,
                0, curr_limbs, curr_limbs, context.K, scheme);
    ApproxModDown(res_d2, d2_expand, curr_limbs,
                  curr_limbs, 0, curr_limbs, context.K,
                  0, curr_limbs, 0, 0, scheme);
    timeutils.stop("ModUpAndDown");
    delete[] cipher2.bx;
    cipher2.bx = res_d2;

    StringUtils::showcompare(d1, res_d2, slots, "ModDown");

    timeutils.start("Decrypt batch");
    complex<double> *dvecMod = scheme.decrypt(secretKey, cipher1);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvec1, dvecMod, slots, "ModUp and ModDown");

}

void TestScheme::KeySwitchingExtTrue(long logN, long L, long logp, long dnum, long logSlots)
{
    cout << "!!! START TEST KeySwitchingExtTrue(OpenFHE) !!!" << endl;
    //-----------------------------------------
    long K = L/dnum;
    Context context(logN, logp, L, K,dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = (1 << logSlots);
    complex<double>* mvec = new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        mvec[i]=i*0.1;
    }

    Ciphertext cipher = scheme.encrypt(mvec, slots, L);

    //----------------------------------------
    // Test for KeySwitchExt + KeySwitchDown
    //----------------------------------------
    long curr_limbs=cipher.curr_limbs;
    uint32_t tmpExt_len= ((curr_limbs+K) << logN);;
    uint64_t* tmpExt_ax = new uint64_t[tmpExt_len]();
    uint64_t* tmpExt_bx = new uint64_t[tmpExt_len]();
    Ciphertext cipherExt(tmpExt_ax, tmpExt_bx, cipher.N, cipher.slots, (curr_limbs+K));
    KeySwitchExt(cipherExt,cipher,2,true,scheme);

    uint32_t tmp_len= (curr_limbs << logN);;
    uint64_t* tmp_ax = new uint64_t[tmp_len]();
    uint64_t* tmp_bx = new uint64_t[tmp_len]();
    Ciphertext result(tmp_ax,tmp_bx,cipher.N,cipher.slots,cipher.curr_limbs);
    KeySwitchDown(result.ax, result.bx, cipherExt.ax, cipherExt.bx,
                  curr_limbs, scheme);

    complex<double>* dvec = scheme.decrypt(secretKey, result);
    StringUtils::showcompare(mvec, dvec, slots, "KeySwitchingExtTrue");
}

void TestScheme::KeySwitchingExtFalse(long logN, long L, long logp, long dnum, long logSlots)
{
    cout << "!!! START TEST KeySwitchingExtFalse(OpenFHE) !!!" << endl;
    //-----------------------------------------
    long K = L/dnum;
    Context context(logN, logp, L, K,dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = (1 << logSlots);
    complex<double>* mvec = new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        mvec[i]=i*0.1;
    }

    Ciphertext cipher = scheme.encrypt(mvec, slots, L);

    //----------------------------------------
    // Test for KeySwitchExt + KeySwitchDown w/o first element
    //----------------------------------------
    long curr_limbs=cipher.curr_limbs;
    uint32_t tmpExt_len= ((curr_limbs+K) << logN);;
    uint64_t* tmpExt_ax = new uint64_t[tmpExt_len]();
    uint64_t* tmpExt_bx = new uint64_t[tmpExt_len]();
    Ciphertext cipherExt(tmpExt_ax, tmpExt_bx, cipher.N, cipher.slots, (curr_limbs+K));
    KeySwitchExt(cipherExt,cipher,2,false,scheme);

    uint32_t tmp_len= (curr_limbs << logN);;
    uint64_t* tmp_ax = new uint64_t[tmp_len]();
    uint64_t *tmp_bx = new uint64_t[tmp_len]();
    Ciphertext result(tmp_ax, tmp_bx, cipher.N, cipher.slots, cipher.curr_limbs);
    KeySwitchDown(result.ax, result.bx, cipherExt.ax, cipherExt.bx, curr_limbs, scheme);
    for (int i = 0; i < (curr_limbs << logN); ++i) {
        result.bx[i] += cipher.bx[i];
    }

    complex<double> *dvec = scheme.decrypt(secretKey, result);
    StringUtils::showcompare(mvec, dvec, slots, "KeySwitchingExtFalse");
}

void TestScheme::testMultByIntegerInPlace() {
    long logN = 6;
    long L = 6;
    long logp = 55;
    long dnum = 2;
    long logSlots = 2;
    long K = L / dnum;
    Context context(logN, logp, L, K, dnum);
    cout << "N: " << context.N << endl
         << "L: " << context.L << endl
         << "dnum: " << context.dnum << endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = (1 << logSlots);
    complex<double> *mvec = new complex<double>[slots];
    complex<double> *golden_value = new complex<double>[slots];
    for (int i = 0; i < slots; ++i) {
        mvec[i] = i * 0.1;
        complex<double> tmp = 2;
        golden_value[i] = mvec[i] * tmp;
    }

    Ciphertext cipher = scheme.encrypt(mvec, slots, L);
    scheme.MultByIntegerInPlace(cipher, 2);
    complex<double> *dvec = scheme.decrypt(secretKey, cipher);
    StringUtils::showcompare(golden_value, dvec, slots, "MultByIntegerInPlace");
}

void TestScheme::testConstMult(long logN, long L, long logp, long logSlots) {
    cout << "!!! START multByConst  and multByConstVec!!!" << endl;
    //-----------------------------------------
    long K = L ;
    Context context(logN, logp, L, K);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------

    const long slots = 1 << logSlots;
    double bound = 1.0;
    complex<double>* mvec1 = new complex<double>[slots];
    for (long i = 0; i < slots; ++i)
    {
        mvec1[i] = complex<double>(1, 1);
    }
    complex<double> cnst(0.8,0.2);
    complex<double>* golden_value = new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        golden_value[i] = mvec1[i] * cnst;
    }
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cmultCipher = scheme.Lattigo_MultByConst(cipher1, cnst); // 等价于，mvec1 每个元素放大 cvec[0]倍
    scheme.reScaleByAndEqual(cmultCipher, 1);
    complex<double> *dvecCMult = scheme.decrypt(secretKey, cmultCipher);
    StringUtils::showcompare(golden_value, dvecCMult, slots, "multByConst");
}

void TestScheme::testmyRescale(long logN, long L, long logp, long logSlots, long dnum) {
    cout << "!!! START TEST HMult !!!" << endl;
    //-----------------------------------------
    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    long slots = (1 << logSlots);
    double bound = 1.0;
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);

    complex<double>* mvecMult = new complex<double>[slots];

    for(long i = 0; i < slots; i++) {
        mvecMult[i] = mvec1[i] * mvec2[i];
    }

    timeutils.start("Encrypt two batch");
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    timeutils.stop("Encrypt two batch");

    timeutils.start("Homomorphic Multiplication");
    Ciphertext multCipher(cipher1.N, cipher1.slots, cipher1.curr_limbs);
    myMult(multCipher,cipher1, cipher2, scheme);
    timeutils.stop("Homomorphic Multiplication");
    timeutils.start("Rescaling");
    ModReduceInternalInPlace(multCipher, 1,scheme);
    timeutils.stop("Rescaling");

    timeutils.start("Decrypt batch");
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);
    timeutils.stop("Decrypt batch");

    StringUtils::showcompare(mvecMult, dvecMult, slots, "HMult1");
}