/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_SCHEME_H_
#define HEAANNTT_SCHEME_H_

#include <map>
#include <chrono>

#include "Common.h"
#include "Ciphertext.h"
#include "Context.h"
#include "Plaintext.h"
#include "SecretKey.h"
#include "Key.h"
#include "Numb.h"
#include "constants.h"
#include "ckksrns-fhe.h"

using namespace std;

static long ENCRYPTION = 0;
static long MULTIPLICATION  = 1;
static long CONJUGATION = 2;

class Scheme {
public:

	Context& context;

    SecretKey& secretKey; // 为了解决multConst 需要通过先解密再加密实现这一问题。

	map<long, Key> keyMap; ///< contain Encryption, Multiplication and Conjugation keys, if generated
	map<long, Key> leftRotKeyMap; ///< contain left rotation keys, if generated
	map<long, Key> BootstrapFastRotationKeyMap; ///< contain bootstrap rotation keys, if generated

    CKKSBootstrapPrecom precom; //fixme: bootstrap 相关参数，在HLS中应该为全局参数。
    ScalingTechnique rescaleTech;
    SecretKeyDist secretKeyDist;
    uint64_t L0;

//	Scheme(Context& context);

	Scheme(SecretKey& secretKey, Context& context);

    Scheme(SecretKey& secretKey, Context& context,
           ScalingTechnique rescaleTech,SecretKeyDist secretKeyDist, uint64_t L0=0);

	/**
	 * generates key for public encryption (key is stored in keyMap)
	 */
	void addEncKey(SecretKey& secretKey);

	/**
	 * generates key for multiplication (key is stored in keyMap)
	 */
	void addMultKey(SecretKey& secretKey);

	/**
	 * generates key for conjugation (key is stored in keyMap)
	 */
	void addConjKey(SecretKey& secretKey);

	/**
	 * generates key for left rotation (key is stored in leftRotKeyMap)
	 */
	void addLeftRotKey(SecretKey& secretKey, long rot);

	/**
	 * generates all keys for power-of-two left rotations (keys are stored in leftRotKeyMap)
	 */
	void addLeftRotKeys(SecretKey& secretKey);

	/**
	 * generates all keys for power-of-two right rotations (keys are stored in leftRotKeyMap)
	 */
	void addRightRotKeys(SecretKey& secretKey);

	Plaintext encode(double* vals, long slots, long l);

	Plaintext encode(complex<double>* vals, long slots, long l);

	Plaintext encodeSingle(complex<double> val, long l);

	complex<double>* decode(Plaintext& msg);

	complex<double> decodeSingle(Plaintext& msg);

	// Encryption (secret and public version)
	Ciphertext encryptMsg(SecretKey& secretkey, Plaintext& message);

	Ciphertext encryptMsg(Plaintext& message);

	Plaintext decryptMsg(SecretKey& secretkey, Ciphertext& cipher);

	Ciphertext encrypt(double* vals, long slots, long l);

	Ciphertext encrypt(complex<double>* vals, long slots, long l);

	Ciphertext encryptSingle(complex<double> val, long l);

	complex<double>* decrypt(SecretKey& secretKey, Ciphertext& cipher);

	complex<double> decryptSingle(SecretKey& secretKey, Ciphertext& cipher);

	// Homomorphic Negation
	Ciphertext negate(Ciphertext& cipher);
	void negateAndEqual(Ciphertext& cipher);

	// Homomorphic Addition
	Ciphertext add(Ciphertext& cipher1, Ciphertext& cipher2);
	void addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	// Homomorphic Substraction
	Ciphertext sub(Ciphertext& cipher1, Ciphertext& cipher2);

	void subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);
	void sub2AndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	// Homomorphic Multiplication
	Ciphertext mult(Ciphertext& cipher1, Ciphertext& cipher2);
	void multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);

	// Homomorphic Squaring
	Ciphertext square(Ciphertext& cipher);
	void squareAndEqual(Ciphertext& cipher);

	// Homomorphic Operations with Constant

	Ciphertext imult(Ciphertext& cipher);
	void imultAndEqual(Ciphertext& cipher);

	Ciphertext idiv(Ciphertext& cipher);
	void idivAndEqual(Ciphertext& cipher);

	Ciphertext addConst(Ciphertext& cipher, double cnst);
	Ciphertext addConst(Ciphertext& cipher, complex<double> cnst);

	void addConstAndEqual(Ciphertext& cipher, double cnst);
	void addConstAndEqual(Ciphertext& cipher, complex<double> cnst);

	void addPcAndEqual(Ciphertext& cipher);
	void addP2AndEqual(Ciphertext& cipher);
	void addP2hAndEqual(Ciphertext& cipher);

	Ciphertext multByConst(Ciphertext& cipher, double cnst);
	Ciphertext multByConst(Ciphertext& cipher, complex<double> cnst);

	Ciphertext multByConstVec(Ciphertext& cipher, double* cnstVec, long slots);
	Ciphertext multByConstVec(Ciphertext& cipher, complex<double>* cnstVec, long slots);

	void multByConstAndEqual(Ciphertext& cipher, double cnst);
	void multByConstAndEqual(Ciphertext& cipher, complex<double> cnst);

	void multByPolyAndEqual(Ciphertext& cipher, uint64_t* poly);

    void multByMonomial(Ciphertext &result, Ciphertext &cipher, long monomialDegree);
	void multByMonomialAndEqual(Ciphertext& cipher, long monomialDegree);

	Ciphertext reScaleBy(Ciphertext& cipher, long dl);
	void reScaleByAndEqual(Ciphertext& cipher, long dl);

	Ciphertext reScaleTo(Ciphertext& cipher, long l);
	void reScaleToAndEqual(Ciphertext& cipher, long l);

	Ciphertext modDownBy(Ciphertext& cipher, long dl);
	void modDownByAndEqual(Ciphertext& cipher, long dl);

	Ciphertext modDownTo(Ciphertext& cipher, long dl);
	void modDownToAndEqual(Ciphertext& cipher, long dl);

	Ciphertext leftRotateFast(Ciphertext& cipher, long rotSlots);
	void leftRotateAndEqualFast(Ciphertext& cipher, long rotSlots);

	Ciphertext leftRotateByPo2(Ciphertext& cipher, long logRotSlots);
	void leftRotateByPo2AndEqual(Ciphertext& cipher, long logRotSlots);

	Ciphertext rightRotateByPo2(Ciphertext& cipher, long logRotSlots);
	void rightRotateByPo2AndEqual(Ciphertext& cipher, long logRotSlots);

	Ciphertext leftRotate(Ciphertext& cipher, long rotSlots);

    void leftRotateAndEqual(Ciphertext &cipher, long rotSlots);

    Ciphertext rightRotate(Ciphertext &cipher, long rotSlots);

    void rightRotateAndEqual(Ciphertext &cipher, long rotSlots);

    Ciphertext conjugate(Ciphertext &cipher);

    void conjugateAndEqual(Ciphertext &cipher);

    //MY AUX procedure
    void EvalMultExtInPlace(Ciphertext &cipher, Plaintext pt);

    void EvalMultExt(Ciphertext &result, const Ciphertext &cipher, Plaintext pt);

    void EvalAddExtInPlace(Ciphertext &cipher1, const Ciphertext &cipher2);

    void MultByIntegerInPlace(Ciphertext &cipher, uint64_t integer);

    //KS key gen procedure
    //core part of switchKey generation procedure
    void KS_KeyGen(uint64_t *afterKS_sk, uint64_t *beforeKS_sk, long type, map<long, Key> &keyMap);


    void ReduceLvl(uint64_t level, uint64_t *p1, uint64_t *p2);

    complex<double> *Lattigo_decrypt(SecretKey &secretKey, Ciphertext &cipher);

    Ciphertext CZHMultByConst(Ciphertext &in, complex<double> cnst);

    Ciphertext Lattigo_MultByConst(Ciphertext &in, complex<double> cnst);

    complex<double> *fulldecrypt(SecretKey &secretKey, Ciphertext &cipher);

    Plaintext fulldecryptMsg(SecretKey &secretKey, Ciphertext &cipher);

    complex<double> *fulldecode(Plaintext &msg);

    Ciphertext fullencrypt(complex<double> *vals, long slots, long l);

    Plaintext fullencode(complex<double> *v, long slots, long l);


};

#endif
