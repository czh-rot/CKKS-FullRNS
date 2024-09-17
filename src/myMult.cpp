//
// Created by EYx on 2023/3/15.
//

#include <stdint.h>

#include "Scheme.h"
#include "TimeUtils.h"
#include "myKeySwitch.h"
#include "myMult.h"

////myMult 写法1：按Better拆分子模块
// 按照Better bootstrapping for approximate homomorphic encryption论文 p11-12的流程拆分模块
void myMult(Ciphertext& result, Ciphertext cipher1, Ciphertext cipher2, Scheme scheme)
{
    uint32_t alpha=scheme.context.alpha;
    uint32_t logN=scheme.context.logN;
    uint32_t K=scheme.context.K;

    long curr_limbs=cipher1.curr_limbs;
    uint64_t* axbx1 = new uint64_t[curr_limbs << logN]();
    uint64_t* axbx2 = new uint64_t[curr_limbs << logN]();
    uint64_t* axax = new uint64_t[curr_limbs << logN]();
    uint64_t* bxbx = new uint64_t[curr_limbs << logN]();

    //1. MultCore
    myMultCore(cipher1, cipher2, axbx1, axbx2, axax, bxbx, scheme);

    //precompute KeySwitchCore params
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //preallocate KeySwitchCore memory
    //total beta groups,
    // each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys
    // the size of each poly is N
    uint64_t *d2 = new uint64_t[curr_limbs << logN];
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *d2Tilde = new uint64_t[beta * expand_length];
    uint64_t *sumaxmult = new uint64_t[expand_length]();
    uint64_t *sumbxmult = new uint64_t[expand_length]();

    //Set result pointer
    uint64_t *res_ax = result.ax;
    uint64_t *res_bx = result.bx;

    //2. RNS-decompose:
    // including 2-1 zero-padding and split, and 2-2 RNS-decompose
    KS_RNS_Decompose(d2, axax, curr_limbs, beta, scheme);

    //3. Modulus Raise
    KS_Modulus_Raise(d2Tilde, d2, beta, curr_limbs, scheme);
    delete[] d2;

    //4. InnerProduct
    EvalFastKeySwitchCoreExt(sumaxmult, sumbxmult, d2Tilde, MULTIPLICATION, scheme.keyMap,
                             expand_length, beta, curr_limbs, scheme);
    delete[] d2Tilde;

    //5. Modulus Down
    KeySwitchDown(res_ax, res_bx, sumaxmult, sumbxmult,
                  curr_limbs, scheme);
    delete[] sumaxmult;
    delete[] sumbxmult;

    //6. output ct
    // post comput: 凑数
    scheme.context.addAndEqual(res_ax, axbx1, curr_limbs);
    scheme.context.subAndEqual(res_ax, bxbx, curr_limbs);
    scheme.context.subAndEqual(res_ax, axax, curr_limbs);
    scheme.context.addAndEqual(res_bx, bxbx, curr_limbs);
    delete[] axax;
    delete[] bxbx;
    delete[] axbx1;
    delete[] axbx2;
}

////myMult 写法2：按openfhe设计封装程度较高的 KeySwitchCore 函数。
//// 按照OpenFHE的方式对KS procedure进行一些封装，为后续针对 KS子函数设计的 eval fast rotation做准备
//void myMult(Ciphertext& result, Ciphertext cipher1, Ciphertext cipher2, Scheme scheme) {
//
//    uint32_t logN=scheme.context.logN;
//
//    long curr_limbs=cipher1.curr_limbs;
//    uint64_t* axbx1 = new uint64_t[curr_limbs << logN]();
//    uint64_t* axbx2 = new uint64_t[curr_limbs << logN]();
//    uint64_t* axax = new uint64_t[curr_limbs << logN]();
//    uint64_t* bxbx = new uint64_t[curr_limbs << logN]();
//
//    //1. MultCore
//    myMultCore(cipher1,cipher2, axbx1,axbx2,axax,bxbx,scheme);
//
//    //2.-5. KeySwitch Procedure
//    uint64_t *res_ax=result.ax;
//    uint64_t *res_bx=result.bx;
//    KeySwitchCore(res_ax, res_bx, axax, MULTIPLICATION, scheme.keyMap, curr_limbs, scheme);
//
//    //6. output ct
//    // post comput: 凑数
//    scheme.context.addAndEqual(res_ax, axbx1, curr_limbs);
//    scheme.context.subAndEqual(res_ax, bxbx, curr_limbs);
//    scheme.context.subAndEqual(res_ax, axax, curr_limbs);
//    scheme.context.addAndEqual(res_bx, bxbx, curr_limbs);
//    delete[] axax;
//    delete[] bxbx;
//    delete[] axbx1;
//    delete[] axbx2;
//}

void myMultCore(Ciphertext &cipher1, Ciphertext &cipher2,
                        uint64_t* axbx1, uint64_t* axbx2, uint64_t* axax, uint64_t* bxbx,Scheme scheme)
{
    scheme.context.add(axbx1, cipher1.ax, cipher1.bx, cipher1.curr_limbs);
    scheme.context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.curr_limbs);
    scheme.context.mulAndEqual(axbx1, axbx2, cipher1.curr_limbs);
    scheme.context.mul(bxbx, cipher1.bx, cipher2.bx, cipher1.curr_limbs);
    scheme.context.mul(axax, cipher1.ax, cipher2.ax, cipher1.curr_limbs);
}

void myMultAndEqual(Ciphertext& cipher1, Ciphertext cipher2, Scheme scheme) {
    uint32_t alpha=scheme.context.alpha;
    uint32_t logN=scheme.context.logN;
    uint32_t K=scheme.context.K;

    long curr_limbs=cipher1.curr_limbs;
    uint64_t* axbx1 = new uint64_t[curr_limbs << logN]();
    uint64_t* axbx2 = new uint64_t[curr_limbs << logN]();
    uint64_t *axax = new uint64_t[curr_limbs << logN]();
    uint64_t *bxbx = new uint64_t[curr_limbs << logN]();

    //1. MultCore
    myMultCore(cipher1, cipher2, axbx1, axbx2, axax, bxbx, scheme);

    //precompute KeySwitchCore params
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //preallocate KeySwitchCore memory
    //total beta groups, each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys, the size of each poly is N
    uint64_t *d2 = new uint64_t[curr_limbs << logN];
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *d2Tilde = new uint64_t[beta * expand_length];
    uint64_t *sumaxmult = new uint64_t[expand_length]();
    uint64_t *sumbxmult = new uint64_t[expand_length]();

    //Set result pointer
    uint64_t *res_ax = cipher1.ax;
    uint64_t *res_bx = cipher1.bx;

    //2. RNS-decompose:
    // including 2-1 zero-padding and split, and 2-2 RNS-decompose
    KS_RNS_Decompose(d2, axax, curr_limbs, beta, scheme);

    //3. Modulus Raise
    KS_Modulus_Raise(d2Tilde, d2, beta, curr_limbs, scheme);
    delete[] d2;

    //4. InnerProduct
    EvalFastKeySwitchCoreExt(sumaxmult, sumbxmult, d2Tilde, MULTIPLICATION, scheme.keyMap,
                             expand_length, beta, curr_limbs, scheme);
    delete[] d2Tilde;

    //5. Modulus Down
    KeySwitchDown(res_ax, res_bx, sumaxmult, sumbxmult,
                  curr_limbs, scheme);
    delete[] sumaxmult;
    delete[] sumbxmult;

    //6. output ct
    // post comput: 凑数
    scheme.context.addAndEqual(cipher1.ax, axbx1, cipher1.curr_limbs);  //凑数字
    scheme.context.subAndEqual(cipher1.ax, bxbx, cipher1.curr_limbs);
    scheme.context.subAndEqual(cipher1.ax, axax, cipher1.curr_limbs);
    scheme.context.addAndEqual(cipher1.bx, bxbx, cipher1.curr_limbs);

    delete[] axax;
    delete[] bxbx;
    delete[] axbx1;
    delete[] axbx2;
}

////mySquare 写法1：按Better拆分子模块
// 按照Better bootstrapping for approximate homomorphic encryption论文 p11-12的流程拆分模块
void mySquare(Ciphertext& result, Ciphertext cipher, Scheme scheme)
{
    uint32_t alpha=scheme.context.alpha;
    uint32_t logN=scheme.context.logN;
    uint32_t K=scheme.context.K;

    long curr_limbs=cipher.curr_limbs;
    uint64_t* axbx = new uint64_t[curr_limbs << logN]();

    uint64_t* axax = new uint64_t[curr_limbs << logN]();
    uint64_t* bxbx = new uint64_t[curr_limbs << logN]();

    //1. MultCore
    scheme.context.add(axbx, cipher.ax, cipher.bx, cipher.curr_limbs); // ax1 + bx1 mod P, 0 mod Q

    scheme.context.squareAndEqual(axbx, cipher.curr_limbs); // (ax1 + bx1) * (ax2 + bx2) mod P, 0 mod Q
    scheme.context.square(bxbx, cipher.bx, cipher.curr_limbs);
    scheme.context.square(axax, cipher.ax, cipher.curr_limbs);

    //precompute KeySwitchCore params
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //preallocate KeySwitchCore memory
    //total beta groups,
    // each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys
    // the size of each poly is N
    uint64_t *d2 = new uint64_t[curr_limbs << logN];
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *d2Tilde = new uint64_t[beta * expand_length];
    uint64_t *sumaxmult = new uint64_t[expand_length]();
    uint64_t *sumbxmult = new uint64_t[expand_length]();

    //Set result pointer
    uint64_t *res_ax = result.ax;
    uint64_t *res_bx = result.bx;

    //2. RNS-decompose:
    // including 2-1 zero-padding and split, and 2-2 RNS-decompose
    KS_RNS_Decompose(d2, axax, curr_limbs, beta, scheme);

    //3. Modulus Raise
    KS_Modulus_Raise(d2Tilde, d2, beta, curr_limbs, scheme);
    delete[] d2;

    //4. InnerProduct
    EvalFastKeySwitchCoreExt(sumaxmult, sumbxmult, d2Tilde, MULTIPLICATION, scheme.keyMap,
                             expand_length, beta, curr_limbs, scheme);
    delete[] d2Tilde;

    //5. Modulus Down
    KeySwitchDown(res_ax, res_bx, sumaxmult, sumbxmult,
                  curr_limbs, scheme);
    delete[] sumaxmult;
    delete[] sumbxmult;

    //6. output ct
    // post comput: 凑数
    scheme.context.addAndEqual(res_ax, axbx, curr_limbs);
    scheme.context.subAndEqual(res_ax, bxbx, curr_limbs);
    scheme.context.subAndEqual(res_ax, axax, curr_limbs);
    scheme.context.addAndEqual(res_bx, bxbx, curr_limbs);
    delete[] axax;
    delete[] bxbx;
    delete[] axbx;
}

////mySquare 写法1：按Better拆分子模块
// 按照Better bootstrapping for approximate homomorphic encryption论文 p11-12的流程拆分模块
void mySquareAndEqual(Ciphertext& cipher, Scheme scheme)
{
    uint32_t alpha=scheme.context.alpha;
    uint32_t logN=scheme.context.logN;
    uint32_t K=scheme.context.K;

    long curr_limbs=cipher.curr_limbs;
    uint64_t* axbx = new uint64_t[curr_limbs << logN]();

    uint64_t* axax = new uint64_t[curr_limbs << logN]();
    uint64_t* bxbx = new uint64_t[curr_limbs << logN]();

    //1. MultCore
    scheme.context.add(axbx, cipher.ax, cipher.bx, cipher.curr_limbs); // ax1 + bx1 mod P, 0 mod Q

    scheme.context.squareAndEqual(axbx, cipher.curr_limbs); // (ax1 + bx1) * (ax2 + bx2) mod P, 0 mod Q
    scheme.context.square(bxbx, cipher.bx, cipher.curr_limbs);
    scheme.context.square(axax, cipher.ax, cipher.curr_limbs);

    //precompute KeySwitchCore params
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //preallocate KeySwitchCore memory
    //total beta groups,
    // each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys
    // the size of each poly is N
    uint64_t *d2 = new uint64_t[curr_limbs << logN];
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *d2Tilde = new uint64_t[beta * expand_length];
    uint64_t *sumaxmult = new uint64_t[expand_length]();
    uint64_t *sumbxmult = new uint64_t[expand_length]();

    //Set result pointer
    uint64_t *res_ax = cipher.ax;
    uint64_t *res_bx = cipher.bx;

    //2. RNS-decompose:
    // including 2-1 zero-padding and split, and 2-2 RNS-decompose
    KS_RNS_Decompose(d2, axax, curr_limbs, beta, scheme);

    //3. Modulus Raise
    KS_Modulus_Raise(d2Tilde, d2, beta, curr_limbs, scheme);
    delete[] d2;

    //4. InnerProduct
    EvalFastKeySwitchCoreExt(sumaxmult, sumbxmult, d2Tilde, MULTIPLICATION, scheme.keyMap,
                             expand_length, beta, curr_limbs, scheme);
    delete[] d2Tilde;

    //5. Modulus Down
    KeySwitchDown(res_ax, res_bx, sumaxmult, sumbxmult,
                  curr_limbs, scheme);
    delete[] sumaxmult;
    delete[] sumbxmult;

    //6. output ct
    // post comput: 凑数
    scheme.context.addAndEqual(res_ax, axbx, curr_limbs);
    scheme.context.subAndEqual(res_ax, bxbx, curr_limbs);
    scheme.context.subAndEqual(res_ax, axax, curr_limbs);
    scheme.context.addAndEqual(res_bx, bxbx, curr_limbs);
    delete[] axax;
    delete[] bxbx;
    delete[] axbx;
}