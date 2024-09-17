/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_CONTEXT_H_
#define HEAANNTT_CONTEXT_H_

#include <complex>
#include <chrono>
#include <map>
#include <boost/multiprecision/number.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include "Common.h"
#include "Numb.h"

#define Q0_BIT_SIZE 61

using namespace std;
using namespace boost::multiprecision;

static string LOGARITHM = "Logarithm"; ///< log(x)
static string EXPONENT  = "Exponent"; ///< exp(x)
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))

class Context {
public:


	// Encryption parameters
	long logN; ///< Logarithm of Ring Dimension
	long logNh; ///< Logarithm of Ring Dimension - 1
	long L; ///< Maximum Level that we want to support
	long K; ///< The number of special modulus (usually L + 1)

	long N;
	long M;
	long Nh;

    //KeySwitch params
    long dnum; //decomposition number
    long alpha; //alpha=(L+1)/dnum
    //KSGen params
    uint64_t** QjHatModqi; //Better, p11
    uint64_t** QjHatInvModqi;
    //Basis conversion precompute params
    //todo: yx: basis conversion预计算，并确认正确性
    //Over 100x, Dprime for double prime, Tprime for triple prime
    //j对应第一维，i对应第二维
    uint64_t** QjHatDprimeModqi; //in_C_L_len>0, out_C_L_len>0
    uint64_t** QjHatDprimeModpi; //in_C_L_len>0, out_B_len>0

    uint64_t** QjHatTprimeModqi; //in_B_len>0, out_C_L_len>0
//    uint64_t** QjHatTprimeModpi; //in_B_len>0, out_B_len>0，不存在这个case

    uint64_t* QjHatDPrimeInvModqj; //in_C_L_len>0
    uint64_t* QjHatTPrimeInvModpj; //in_B_L_len>0
    //ApproxModDown
    uint64_t* QHatInvModqj;

    //Rescale openfhe
    uint64_t ** QlQlInvModqlDivqlModq;

    //Encode lattigo
    cpp_int ModulusBigint; // new add encode
    cpp_int* bigintChain; // new add encode

	long logp;
	long p;

	long h;
	double sigma;

	uint64_t* qVec;
	uint64_t* pVec;

	uint64_t* qrVec; // Barrett reduction
	uint64_t* prVec; // Barrett recution

	long* qTwok; // Barrett reduction
	long* pTwok; // Barrett reduction

	uint64_t* qkVec; // Montgomery reduction
	uint64_t* pkVec; // Montgomery reduction

	uint64_t* qdVec;
	uint64_t* pdVec;

	uint64_t* qInvVec;
	uint64_t* pInvVec;

	uint64_t* qRoots;
	uint64_t* pRoots;

	uint64_t* qRootsInv;
	uint64_t* pRootsInv;

	uint64_t** qRootPows;
	uint64_t** pRootPows;

	uint64_t** qRootPowsInv;
	uint64_t** pRootPowsInv;

	uint64_t* NInvModq;
	uint64_t* NInvModp;

	uint64_t** qRootScalePows;
	uint64_t** pRootScalePows;

	uint64_t** qRootScalePowsOverq;
	uint64_t** pRootScalePowsOverp;

	uint64_t** qRootScalePowsInv;
	uint64_t** pRootScalePowsInv;

	uint64_t* NScaleInvModq; // [i]
	uint64_t* NScaleInvModp; // [k]

	uint64_t** qHatModq; // [curr_limbs][i] (phat_i)_l mod p_i
	uint64_t* pHatModp; // [k] qhat_k mod q_k

	uint64_t** qHatInvModq; // [curr_limbs][i] (qhat_i)_l^-1 mod q_i
	uint64_t* pHatInvModp; // [k] phat_k^-1 mod p_k

	uint64_t*** qHatModp; // [curr_limbs] [i] [k]  (phat_i)_l mod q_k

	uint64_t** pHatModq; // [k][i] qhat_k mod p_i

	uint64_t* PModq; // [i] qprod mod p_i
	uint64_t* PInvModq; // [i] qprod mod p_i

	uint64_t** QModp; // [i] qprod mod p_i
	uint64_t** QInvModp; // [i] qprod mod p_i

	uint64_t** qInvModq; // [i][j] p_i^-1 mod p_j

	long* rotGroup; ///< precomputed rotation group indexes

	complex<double>* ksiPows; ///< precomputed ksi powers

	map<string, double*> taylorCoeffsMap; ///< precomputed taylor coefficients

    uint64_t **nttPsi; //powers of the inverse of the 2nth primitive root in Montgomery form (in bitreversed order)
    uint64_t **nttPsiInv; //powers of the inverse of the 2nth primitive root in Montgomery form (in bitreversed order)
    uint64_t *psiMont;
    uint64_t *psiInvMont;

	uint64_t* p2coeff;
	uint64_t* pccoeff;
	uint64_t* p2hcoeff;

    //by default, the full rns variant of ckks proj dnum is 1
	Context(long logN, long logp, long L, long K, long dnum=1, long h = 64, double sigma = 3.2);
	Context(string STRING, long logN, long logp, long L, long K, long dnum, long h = 64, double sigma = 3.2);

	void arrayBitReverse(complex<double>* vals, const long size);
	void arrayBitReverse(uint64_t* vals, const long size);

	void fft(complex<double>* vals, const long size);
	void fftInvLazy(complex<double>* vals, const long size);
	void fftInv(complex<double>* vals, const long size);

	void fftSpecial(complex<double>* vals, const long size);
	void fftSpecialInvLazy(complex<double>* vals, const long size);
	void fftSpecialInv(complex<double>* vals, const long size);

	void encode(uint64_t* ax, complex<double>* vals, long slots, long l);
	void encode(uint64_t* ax, double* vals, long slots, long l);

	void encodeSingle(uint64_t* ax, complex<double>& val, long l);
	void encodeSingle(uint64_t* ax, double val, long l);

	void decode(uint64_t* ax, complex<double>* vals, long slots, long l);
	void decode(uint64_t* ax, double* vals, long slots, long l);

	void decodeSingle(uint64_t* ax, complex<double>& val, long l);
	void decodeSingle(uint64_t* ax, double val, long l);

	void qiNTT(uint64_t* res, uint64_t* a, long index);
	void piNTT(uint64_t* res, uint64_t* a, long index);

	void NTT(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiNTTAndEqual(uint64_t* a, long index);
	void piNTTAndEqual(uint64_t* a, long index);

	void NTTAndEqual(uint64_t* a, long l, long k = 0);

	void qiINTT(uint64_t* res, uint64_t* a, long index);
	void piINTT(uint64_t* res, uint64_t* a, long index);

	void INTT(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiINTTAndEqual(uint64_t* a, long index);
	void piINTTAndEqual(uint64_t* a, long index);

	void INTTAndEqual(uint64_t* a, long l, long k = 0);

	void qiNegate(uint64_t* res, uint64_t* a, long index);
	void piNegate(uint64_t* res, uint64_t* a, long index);

	void negate(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiNegateAndEqual(uint64_t* a, long index);
	void piNegateAndEqual(uint64_t* a, long index);

	void negateAndEqual(uint64_t* a, long l, long k = 0);

	void qiAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void piAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index);

	void addConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k = 0);

	void qiAddConstAndEqual(uint64_t* a, uint64_t c, long index);
	void piAddConstAndEqual(uint64_t* a, uint64_t c, long index);

	void addConstAndEqual(uint64_t* a, uint64_t c, long l, long k = 0);

	void qiSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void piSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index);

	void subConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k = 0);

	void qiSubConstAndEqual(uint64_t* a, uint64_t c, long index);
	void piSubConstAndEqual(uint64_t* a, uint64_t c, long index);

	void subConstAndEqual(uint64_t* a, uint64_t c, long l, long k = 0);

	void qiAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index);

	void add(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiAddAndEqual(uint64_t* a, uint64_t* b, long index);
	void piAddAndEqual(uint64_t* a, uint64_t* b, long index);

	void addAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSub(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piSub(uint64_t* res, uint64_t* a, uint64_t* b, long index);

	void sub(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSubAndEqual(uint64_t* a, uint64_t* b, long index);
	void piSubAndEqual(uint64_t* a, uint64_t* b, long index);

	void subAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSub2AndEqual(uint64_t* a, uint64_t* b, long index);
	void piSub2AndEqual(uint64_t* a, uint64_t* b, long index);

	void sub2AndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index);
	void piMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index);

	void mulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long l, long k = 0);

	void qiMulConstAndEqual(uint64_t* res, uint64_t cnst, long index);
	void piMulConstAndEqual(uint64_t* res, uint64_t cnst, long index);

	void mulConstAndEqual(uint64_t* res, uint64_t cnst, long l, long k = 0);

	void qiMul(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piMul(uint64_t* res, uint64_t* a, uint64_t* b, long index);

	void mul(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);

	void mulKey(uint64_t* res, uint64_t* a, uint64_t* b, long l);

	void qiMulAndEqual(uint64_t* a, uint64_t* b, long index);
	void piMulAndEqual(uint64_t* a, uint64_t* b, long index);

	void mulAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSquare(uint64_t* res, uint64_t* a, long index);
	void piSquare(uint64_t* res, uint64_t* a, long index);

	void square(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiSquareAndEqual(uint64_t* a, long index);
	void piSquareAndEqual(uint64_t* a, long index);

	void squareAndEqual(uint64_t* a, long l, long k = 0);

	void evalAndEqual(uint64_t* a, long l);

	void raise(uint64_t* res, uint64_t* a, long l);
    void raiseAndEqual(uint64_t*& a, long l);

    void back(uint64_t* res, uint64_t* a, long l);
    void backAndEqual(uint64_t*& a, long input_l);

	void reScale(uint64_t* res, uint64_t* a, long l);
    void reScaleAndEqual(uint64_t*& a, long l);

	uint64_t* modDown(uint64_t* a, long l, long dl);
    void modDownAndEqual(uint64_t*& a, long l, long dl);

	void leftRot(uint64_t* res, uint64_t* a, long l, long rotSlots);
	void leftRotAndEqual(uint64_t* a, long l, long rotSlots);

	void conjugate(uint64_t* res, uint64_t* a, long l);
	void conjugateAndEqual(uint64_t* a, long l);

	void mulByMonomial(uint64_t* res, uint64_t* a, long l, long monomialDeg);
	void mulByMonomialAndEqual(uint64_t* a, long l, long monomialDeg);

	void sampleGauss(uint64_t* res, long l, long k = 0);
	void sampleZO(uint64_t* res, long s, long l, long k = 0);
	void sampleUniform(uint64_t* res, long l, long k = 0);
	void sampleHWT(uint64_t* res, long l, long k = 0);


    void myEvalAndEqual(uint64_t* ra, uint64_t* a, long l, int j=0);
    void myPrecomputation();

    static uint64_t ScaleUpExact(double value, long n, uint64_t q);

    static uint64_t BRedAdd(uint64_t x, uint64_t q, uint64_t *u);

    uint64_t ModExp(uint64_t x, uint64_t e, uint64_t p);

    static uint64_t MRed(uint64_t x, uint64_t y, uint64_t q, uint64_t qInv);

    void MRed2(uint64_t *x, uint64_t y, uint64_t q, uint64_t qInv, uint64_t &r);

    static uint64_t CRed(uint64_t a, uint64_t q);

    static void ComputeBRedParameters(uint64_t q, uint64_t &mhi, uint64_t &mlo);

    int Len64(uint64_t x);

    static uint64_t MForm(uint64_t a, uint64_t q, const uint64_t *u);

    uint64_t BitReverse64(uint64_t index, uint64_t bitLen);

    uint64_t Reverse64(uint64_t x);


    void fulldecode(uint64_t *a, complex<double> *v, long slots, long l);

    void PolyToBigint(uint64_t *p1, cpp_int *coeffsBigint, long level);

    static cpp_int ModInverse(cpp_int g, cpp_int n);

    static cpp_int GCD(cpp_int& x, cpp_int& y, cpp_int a, cpp_int b);

    int cmp(cpp_int a, cpp_int b);

    double scaleDown(cpp_int &coeff, double n);

    void scaleUpVecExact(double *values, double n, uint64_t *moduli, uint64_t *coeffs, long l);

    void fullencode(uint64_t *a, complex<double> *v, long slots, long l);

    };

#endif
