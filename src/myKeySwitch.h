//
// Created by EYx on 2023/4/20.
//

#ifndef FRNS_BOOTSTRAP_MYKEYSWITCH_H
#define FRNS_BOOTSTRAP_MYKEYSWITCH_H
#include "Common.h"
#include "Context.h"
#include "Scheme.h"

/**
* Only supported for hybrid key switching.
* Scales down the polynomial c0 from extended basis P*Q to Q.
*
* @param ciphertext input ciphertext in the extended basis
* @return resulting polynomial
*/
//KeySwitch Down cipher.bx! cipher.bx corresponds to ct.[0] in OpenFHE
void KeySwitchDownFirstElement(uint64_t *res_bx,
                               uint64_t *sumbxmult,
                               long curr_limbs, Scheme scheme);

/**
* Only supported for hybrid key switching.
* Takes a ciphertext in the normal basis Q
* and and extends it to extended basis P*Q.
*
* @param ciphertext input ciphertext in basis Q
* @return resulting ciphertext in basis P*Q
*/
void KeySwitchExt(Ciphertext &result, const Ciphertext &cipher, long cipher_size, bool addFirst, Scheme scheme);

void KS_RNS_Decompose(uint64_t *d2, uint64_t *input,
                      long curr_limbs, long beta, Scheme scheme);

void KS_Modulus_Raise(uint64_t *d2Tilde, uint64_t *d2,
                      long beta, long curr_limbs, Scheme scheme);

//执行 KeySwitchCore 流程中 InnerProduct 的功能。
void EvalFastKeySwitchCoreExt(uint64_t *sumaxmult, uint64_t *sumbxmult, uint64_t *d2Tilde,
                              long type, map<long, Key> &keyMap_,
                              long expand_length, long beta, long curr_limbs, Scheme scheme);

//执行 KeySwitchCore 流程中 Modulus Down 的功能。
void KeySwitchDown(uint64_t *res_ax, uint64_t *res_bx, uint64_t *sumaxmult, uint64_t *sumbxmult,
                   long curr_limbs, Scheme scheme);

//执行 KeySwitchCore 流程中 RNS-decompose, Modulus Raise。
void EvalKeySwitchPrecomputeCore(uint64_t *d2Tilde, uint64_t *input,
                                 long curr_limbs, long beta, Scheme scheme);

//执行 KeySwitchCore 流程中 InnerProduct, Modulus Down。
void EvalFastKeySwitchCore(uint64_t *res_ax, uint64_t *res_bx, uint64_t *d2Tilde,
                           long type, map<long, Key> &keyMap_,
                           long expand_length, long curr_limbs, long beta, Scheme scheme);

//执行 完整的 KS 流程，包括，RNS-decompose, Modulus Raise, InnerProduct, Modulus Down
//note: res_ax 和 res_bx 的空间在调用者处分配
void KeySwitchCore(uint64_t *res_ax, uint64_t *res_bx, uint64_t *axax,
                   long type, map<long, Key> &keyMap_,
                   long curr_limbs, Scheme scheme);

//in_C_L_index: 输入在C_L中的第一个模数的下标；输出同理
//默认先存 mod Q 的，再存 mod P 的。Over 100x 里先存 mod P 的，再存 mod Q 的。
void ApproxModUp(
        uint64_t *ra, uint64_t *a,
        long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_Len,
        long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_Len, Scheme scheme);

void ApproxModDown(
        uint64_t *ra, uint64_t *a, long curr_limbs,
        long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_len,
        long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_len, Scheme scheme);

void BasisConversion(
        uint64_t *ra, uint64_t *a,
        long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_len,
        long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_len, Scheme scheme
);

//long piStart, the index of the first RNS poly mod pi
//in_C_L_index 表示输入在C_L中的第一个模数的下标。
//a 表示待操作多项式起点，i.e. 函数内部不计算偏移
void partialNTTAndEqual(
        uint64_t *a,
        long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_len,
        long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_len, Scheme scheme);

#endif //FRNS_BOOTSTRAP_MYKEYSWITCH_H
