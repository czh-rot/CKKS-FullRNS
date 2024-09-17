//
// Created by EYx on 2023/5/10.
//
#include "../data/testConstValue.h"

#include "Context.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "StringUtils.h"

#include "myTestBootstrap.h"
#include "myBootstrapExample.h"
void testModulusSwitch()
{

    long dnum=3;
    long L=6;
    long K = L/dnum;
    long logN=10;
    long N=1<<logN;
    long logp=55;
    long slots=8;
    Context context(logN,logp,L,K,dnum);
    SecretKey secretKey(context);
    Scheme scheme(secretKey,context);

    uint64_t *ax=new uint64_t [1<<logN];
    uint64_t *bx=new uint64_t [1<<logN];
    uint64_t *ax_q1=new uint64_t [1<<logN];
    uint64_t *bx_q1=new uint64_t [1<<logN];

    for (int i = 0; i < (1 << logN); ++i) {
        ax[i]=ctxtDCRT0[i];
        bx[i]=ctxtDCRT1[i];
    }
    SwitchModulus(ax_q1,36028797018918401,ax,1152921504606844417,N);
    SwitchModulus(bx_q1,36028797018918401,bx,1152921504606844417,N);

    StringUtils::showcompare(ctxtDCRT0_key, ax_q1, slots, "Modulus_Switch");
    StringUtils::showcompare(ctxtDCRT1_key, bx_q1, slots, "Modulus_Switch");
}

void testBootstrap()
{
    long dnum=3;
    long L=6;
    long K = L/dnum;
    long logN=10;
    long N=1<<logN;
    long logp=55;
    long slots=8;
    Context context(logN,logp,L,K,dnum);
    SecretKey secretKey(context);
    Scheme scheme(secretKey,context);

    uint64_t *moduli=scheme.context.qVec;
    uint64_t *roots=scheme.context.qRoots;
    uint32_t sizeQ=scheme.L0;

    std::complex<double> mvec[]{0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888};

    Ciphertext cipher = scheme.encrypt(mvec, slots, L);
    scheme.modDownTo(cipher,1);// ModRaise之前还有一步AdjustCiphertext，其中有一个明密文乘，需要一次ModReduce

//    long L0=L;
    uint64_t *ax=new uint64_t [1<<logN];
    uint64_t *bx=new uint64_t [1<<logN];
    uint64_t *ax_q1=new uint64_t [1<<logN];
    uint64_t *bx_q1=new uint64_t [1<<logN];

    for (int i = 0; i < (1 << logN); ++i) {
        ax[i]=ctxtDCRT0[i];
        bx[i]=ctxtDCRT1[i];
    }
    SwitchModulus(ax_q1,36028797018918401,ax,1152921504606844417,N);
    SwitchModulus(bx_q1,36028797018918401,bx,1152921504606844417,N);

    StringUtils::showcompare(ctxtDCRT0_key, ax_q1, slots, "Modulus_Switch");
    StringUtils::showcompare(ctxtDCRT1_key, bx_q1, slots, "Modulus_Switch");

//    //deal with cipher.ax
//    scheme.context.INTTAndEqual(raised.ax, 1);
//    for (int i = 1; i < 2; ++i) {
//        uint64_t* raised_axj=raised.ax+(i<<logN);
//    }
//    scheme.context.NTTAndEqual(raised.ax, L0);
//
//    //deal with cipher.bx
//    scheme.context.INTTAndEqual(raised.bx, 1);
//    for (int i = 1; i < L0; ++i) {
//        uint64_t* raised_bxj=raised.bx+(i<<logN);
//    }
//    scheme.context.NTTAndEqual(raised.bx, L0);

//    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
//    StringUtils::showcompare(mvec, dvec, slots, "Modulus_Switch");
}