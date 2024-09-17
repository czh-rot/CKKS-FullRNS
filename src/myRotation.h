//
// Created by EYx on 2023/4/18.
//

#ifndef FRNS_BOOTSTRAP_MYROTATION_H
#define FRNS_BOOTSTRAP_MYROTATION_H

#include <stdint.h>

#include "Context.h"
#include "Scheme.h"
#include "EvaluatorUtils.h"
#include "TimeUtils.h"

/*
 * Rotaion Key Genereation procedure for Does Fully Homomorphic Encryption Need Compute Acceleration? Algo4 HRotate
 * @param rot_list rot_list[i]>0 for left rotation,rot_list[i]<0 for right rotation
 */
void HRotate_KeyGen(SecretKey secretKey, long* rot_list, long rot_list_len, Scheme& scheme);

/*
 * Does Fully Homomorphic Encryption Need Compute Acceleration? Algo4 HRotate
 * @param rot_list rot_list[i]>0 for left rotation,rot_list[i]<0 for right rotation
 */
//todo: yx: 注：该函数尚未支持 右轮转，也没有想好右轮转是通过负值表示还是增加一个 bool 表明左右，HRotate_KeyGen也需要相应处理。
void HRotate(Ciphertext*& out_ct, Ciphertext in_ct, long* rot_list, long rot_list_len, Scheme scheme);

void FastRotate_KeyGen_List(SecretKey secretKey, int* index, long index_len, Scheme& scheme);

//compute evk for fastRotation, rotslots 对应的密钥靠 autoIndex=automorph(rotslots) 索引
void FastRotate_KeyGen(SecretKey secretKey, int index, Scheme& scheme);

//先KS 将解密密钥从(s) -> (automorph(s))，再对密文 automorph
//inplace，支持右轮转，
void FastRotate_demo(Ciphertext& cipher, int index, Scheme scheme);

//Note: digits 应当在调用者处分配空间
void EvalFastRotationPrecompute(uint64_t *digits, const Ciphertext &cipher, Scheme scheme);

/*
 * EvalFastRotation implements the automorphism and key swithcing step of
 * hoisted automorphisms.
 *
 * [!!!!CAUTION!!!!]
 * This method assumes that all required rotation keys exist. This may not be
 * true if we are using baby-step/giant-step key switching. Please refer to
 * Section 5.1 of the above reference and EvalPermuteBGStepHoisted to see how
 * to deal with this issue.
 *
 * @param ct the input ciphertext to perform the automorphism on
 * @param index the index of the rotation. Positive indices correspond to left
 * rotations and negative indices correspond to right rotations.
 * @param m is the cyclotomic order
 * @param digits the digit decomposition created by EvalFastRotationPrecompute
 * at the precomputation step.
*/
//NOTE: yx: 1) result 的存储空间需要预先在外部分配; 2) 这个函数要求，必须有对应的rotation密钥！才能正确计算！
void EvalFastRotation(Ciphertext& result, Ciphertext& cipher,
                      uint64_t* digits, int index,
                      Scheme scheme);

/**
* Only supported for hybrid key switching.
* Performs fast (hoisted) rotation and returns the results
* in the extended CRT basis P*Q
*
* @param ciphertext input ciphertext
* @param index the rotation index.
* @param digits the precomputed digits for the ciphertext
* @param addFirst if true, the the first element c0 is also computed (otherwise ignored)
* @return resulting ciphertext
*/
//NOTE: yx: 1) resultExt 的存储空间需要预先在外部分配; 2) 这个函数要求，必须有对应的rotation密钥！才能正确计算！
void EvalFastRotationExt(Ciphertext& resultExt, Ciphertext cipher, int index,
                         uint64_t* digits, bool addFirst, Scheme scheme );

/**
 * Find an automorphism index for a power-of-two cyclotomic order
 * @param i is the plaintext array index
 * @param m is the cyclotomic order
 * @return the automorphism index
 */
/**
* @see FindAutomorphismIndex2n(), version for CKKS
*/
//计算 5^{i} mod m, 若 i<0, 则考虑 (5^{-1})^{i} mod m
//usint LeveledSHECKKSRNS::FindAutomorphismIndex(usint index, usint m)
long FindAutomorphismIndex2nComplex(int32_t i, uint32_t m);

void PrecomputeAutoMap(uint32_t n, uint32_t k, uint32_t* precomp);

/* Function to reverse bits of num */
inline uint64_t ReverseBits(uint64_t num, uint64_t msb);

static int shift_trick[] = {0, 7, 6, 5, 4, 3, 2, 1};


/**
 * Method to reverse bits of num and return an unsigned int, for all bits up to
 * an including the designated most significant bit.
 *
 * @param input an unsigned int
 * @param msb the most significant bit.  All larger bits are disregarded.
 *
 * @return an unsigned integer that represents the reversed bits.
 */
// precomputed reverse of a byte
inline static unsigned char reverse_byte(unsigned char x);

/**
* @brief Performs an automorphism transform operation using precomputed bit
* reversal indices.
* @param l the number of current limbs
* @param N the length of a single RNS poly
* @param i is the element to perform the automorphism transform with.
* @param &vec a vector with precomputed indices
* @return is the result of the automorphism transform.
*/
//外部函数负责分配存储结果的空间
void AutomorphismTransform(uint64_t* ra, uint64_t* a,
                           uint32_t l, uint32_t N,
                           uint32_t i, uint32_t* precomp_vec);


void Conjugate_KeyGen(SecretKey secretKey, Scheme& scheme);

//同 openfhe 的 ckksrns-fhe.cpp 中的实现
void Conjugate_demo(Ciphertext& cipher, Scheme scheme);

#endif //FRNS_BOOTSTRAP_MYROTATION_H
