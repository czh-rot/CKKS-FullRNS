//
// Created by EYx on 2023/4/1.
//


#include "Common.h"
#include "Context.h"
#include "Scheme.h"
#include "myKeySwitch.h"
#include "myUtils.h"

//todo: yx: to be tested
void KeySwitchDownFirstElement(uint64_t *res_bx,
                               uint64_t *sumbxmult,
                               long curr_limbs, Scheme scheme) {
    uint32_t K = scheme.context.K;
    uint32_t logN = scheme.context.logN;

    //SetZero result.ax and result.bx for KeySwitchDown
    long res_len = curr_limbs << logN;
    SetZero(res_bx, 0, res_len);//fixme: res应该是只写入的，所以理论上不需要这个初始化，但是目前不初始化 myMultAndEqual之类的乘法函数就会算错。轮转、共轭、ks辅助函数不需要初始化。

    // 将 curr_limbs+K 降成 curr_limbs 个;
    // 对照 over100x Algo4 写 ModDown
    long in_C_L_len = 0;
    ApproxModDown(res_bx, sumbxmult, curr_limbs,
                  curr_limbs, in_C_L_len, curr_limbs, K,
                  0, curr_limbs, 0, 0, scheme);
}

void KeySwitchExt(Ciphertext &result, const Ciphertext &cipher, long cipher_size, bool addFirst, Scheme scheme) {
    //NOTE: 如果 某次密文乘法后不做重线性导致这里出现三个分量的输入就会出bug，openfhe是可以解决任意个分量的密文的，但是这里没有实现这种通用性
    if (cipher_size > 2) {
        cout << "error, KeySwitchExt 不支持 " << cipher_size << " 个数量的密文" << endl;
    }

    uint32_t curr_limbs = cipher.curr_limbs;
    uint32_t N = cipher.N;
    uint32_t logN = scheme.context.logN;
    uint32_t K = scheme.context.K;

    // sizeCv 即为 cipher_size 一般为2，即，密文有两个分量
//    for (usint k = 0; k < sizeCv; k++) {
//        resultElements[k] = DCRTPoly(paramsQlP, Format::EVALUATION, true);
//        if ((addFirst) || (k > 0)) {
//            auto cMult = cv[k].TimesNoCheck(cryptoParams->GetPModq());
//            for (usint i = 0; i < sizeQl; i++) {
//                resultElements[k].SetElementAtIndex(i, cMult.GetElementAtIndex(i));
//            }
//        }
//    }
    //上述openfhe代码中，初始化resultElements的参数列表中的true表示将element都初始化为0
    // i.e.模 pi 的RNS多项式每项值为 0
    SetZero(result.bx, curr_limbs << logN, (curr_limbs + K) << logN);
    SetZero(result.ax, curr_limbs << logN, (curr_limbs + K) << logN);

    if ((addFirst)) {
        for (int i = 0; i < curr_limbs; i++) {
            uint64_t *result_bxj = result.bx + (i << logN);
            uint64_t *cipher_bxj = cipher.bx + (i << logN);
            uint64_t PModqi = scheme.context.PModq[i];
            for (int j = 0; j < N; ++j) {
                //compute result_bxj[i]=cipher_bxj[i]*PModqi;
                mulModBarrett(result_bxj[j], cipher_bxj[j], PModqi,
                              scheme.context.qVec[i], scheme.context.qrVec[i], scheme.context.qTwok[i]);
            }
        }
    } else {
        SetZero(result.bx, 0, curr_limbs << logN);
    }

    for (int i = 0; i < curr_limbs; i++) {
        uint64_t *result_axj = result.ax + (i << logN);
        uint64_t *cipher_axj = cipher.ax + (i << logN);
        uint64_t PModqi = scheme.context.PModq[i];
        for (int j = 0; j < N; ++j) {
            //compute cMultj[i]=cipher_bxj[i]*PModqi;
            mulModBarrett(result_axj[j], cipher_axj[j], PModqi,
                          scheme.context.qVec[i], scheme.context.qrVec[i], scheme.context.qTwok[i]);
        }
    }
}

void KS_RNS_Decompose(uint64_t *d2, uint64_t *input,
                      long curr_limbs, long beta,
                      Scheme scheme) {

    uint32_t alpha = scheme.context.alpha;
    uint32_t logN = scheme.context.logN;
    uint32_t N = scheme.context.N;
    uint64_t *qVec = scheme.context.qVec;
    uint64_t *qrVec = scheme.context.qrVec;
    long *qTwok = scheme.context.qTwok;
    uint64_t **QjHatInvModqi = scheme.context.QjHatInvModqi;
    //2. RNS-decompose
    //2-1 zero-padding and split, and 2-2 RNS-decompose
    //RNS-decompose: times QjHatInvModqi[j][i]
    for (int j = 0; j < beta; ++j) {
        uint64_t *d2j = d2 + (j * alpha << logN);
        for (int i = 0; i < alpha; ++i) {
            long j_alpha_Plus_i = j * alpha + i;
            int currPolyStart = (i << logN);
            //note: 这里是和 模数总个数 进行比较，区别与Better论文中的和当前level大小比较，（模数总个数=level+1）。所以下方判断条件应该写＞=比较合理。
            //note: 这里参考openfhe，最后一组剩余多少个就处理多少个，不填充至alpha个。
            // 因此，下方的factor中也不乘Better和over100x中提到的 (Q' mod q )
            if (j_alpha_Plus_i >= curr_limbs) { //alpha 不整除(curr_limbs), 不填充
                return;
            } else {   //ja+i < curr_limbs
                uint64_t *axaxj = input + (j * alpha << logN);
                uint64_t factor = QjHatInvModqi[j][i];
                //traverse curr poly
                for (int k = 0; k < N; ++k) {
                    //*(d2j+currPolyStart+k)=*(input+currPolyStart+k)*QjHatInvModqi;
                    uint64_t result = 1;
                    //todo: yx: skip the first mulModBarett if (alpha|(curr_limbs+1))
                    mulModBarrett(result, factor, *(axaxj + currPolyStart + k),
                                  qVec[j_alpha_Plus_i], qrVec[j_alpha_Plus_i], qTwok[j_alpha_Plus_i]);
                    *(d2j + currPolyStart + k) = result;
                }
            }
        }
    }
}

void KS_Modulus_Raise(uint64_t *d2Tilde, uint64_t *d2,
                      long beta, long curr_limbs, Scheme scheme) {
    uint32_t alpha = scheme.context.alpha;
    uint32_t logN = scheme.context.logN;
    uint32_t K = scheme.context.K;

    //每 (alpha 个乘法modulus) 转为 ( curr_limbs个乘法modulus + K个special modulus) 上，
    // i.e. 由 currlimbs->beta*(curr_limbs+K)
    //[!]先存储 (currlimbs个乘法modulus), 再存储 (K个special modulus)
    //本工程默认的raiseAndEqual函数是 a 个升成 a+k 个；我们需要将指定的 连续 alpha 个RNS多项式，升成 curr_limbs+K 个；
    // 对照 Over100x Algo2 写 ModUp
    // out_C_L_len 指升模操作完成后所有模qi的RNS多项式个数
    long expand_length = ((curr_limbs + K) << logN);
    for (int j = 0; j < beta; j++) {
        if (j == (beta - 1)) {
            long in_C_L_len = curr_limbs - j * alpha;
            ApproxModUp(d2Tilde + (j * expand_length),
                        d2 + (j * alpha << logN),
                        j * alpha, in_C_L_len, 0, 0,
                        0, curr_limbs, curr_limbs, K, scheme);
        } else {
            ApproxModUp(d2Tilde + (j * expand_length),
                        d2 + (j * alpha << logN),
                        j * alpha, alpha, 0, 0,
                        0, curr_limbs, curr_limbs, K, scheme);
        }
    }
}

void EvalFastKeySwitchCoreExt(uint64_t *sumaxmult, uint64_t *sumbxmult, uint64_t *d2Tilde,
                              long type, map<long, Key> &keyMap_,
                              long expand_length, long beta, long curr_limbs, Scheme scheme) {
    uint64_t dnum = scheme.context.dnum;
    uint64_t K = scheme.context.K;

    uint64_t *axmult = new uint64_t[expand_length]();
    uint64_t *bxmult = new uint64_t[expand_length]();

    SetZero(sumaxmult, 0, expand_length);
    SetZero(sumbxmult, 0, expand_length);

    for (int j = 0; j < beta; ++j) {
        Key key = keyMap_.at(type * dnum + j); //compute key index
        uint64_t *d2Tildej = d2Tilde + j * expand_length; //compute d2Tildej starting point
        //evk是一个二维向量，在论文 Better xxx 的 Key-Switch 第三步中，evk的两个分量分别和升模后的d2(j)相乘
        //[!]context.mulKey函数默认处理 curr_limbs+k 个元素，i.e.此处的 (curr_limbs+context.K) 个 RNS多项式
        scheme.context.mulKey(axmult, d2Tildej, key.ax, curr_limbs);
        scheme.context.mulKey(bxmult, d2Tildej, key.bx, curr_limbs);
        scheme.context.addAndEqual(sumaxmult, axmult, curr_limbs, K);
        scheme.context.addAndEqual(sumbxmult, bxmult, curr_limbs, K);
    }
    delete[] axmult;
    delete[] bxmult;
}

void KeySwitchDown(uint64_t *res_ax, uint64_t *res_bx,
                   uint64_t *sumaxmult, uint64_t *sumbxmult,
                   long curr_limbs, Scheme scheme) {
    uint32_t K = scheme.context.K;
    uint32_t logN = scheme.context.logN;

    //SetZero result.ax and result.bx for KeySwitchDown
    long res_len = curr_limbs << logN;
    SetZero(res_ax, 0, res_len);//fixme: res应该是只写入的，所以理论上不需要这个初始化，但是目前不初始化 myMultAndEqual之类的乘法函数就会算错。轮转、共轭、ks辅助函数不需要初始化。
    SetZero(res_bx, 0, res_len);//fixme: res应该是只写入的，所以理论上不需要这个初始化，但是目前不初始化 myMultAndEqual之类的乘法函数就会算错。轮转、共轭、ks辅助函数不需要初始化。

    // 将 curr_limbs+K 降成 curr_limbs 个;
    // 对照 over100x Algo4 写 ModDown
    long in_C_L_len = 0;
    ApproxModDown(res_ax, sumaxmult, curr_limbs,
                  curr_limbs, in_C_L_len, curr_limbs, K,
                  0, curr_limbs, 0, 0, scheme);
    ApproxModDown(res_bx, sumbxmult, curr_limbs,
                  curr_limbs, in_C_L_len, curr_limbs, K,
                  0, curr_limbs, 0, 0, scheme);
}


void EvalKeySwitchPrecomputeCore(uint64_t *d2Tilde, uint64_t *input,
                                 long curr_limbs, long beta, Scheme scheme) {
    uint32_t logN = scheme.context.logN;

    uint64_t *d2 = new uint64_t[curr_limbs << logN];
    //2. RNS-decompose:
    // including 2-1 zero-padding and split, and 2-2 RNS-decompose
    KS_RNS_Decompose(d2, input, curr_limbs, beta, scheme);

    //3. Modulus Raise
    KS_Modulus_Raise(d2Tilde, d2, beta, curr_limbs, scheme);

    delete[] d2;
}


void EvalFastKeySwitchCore(uint64_t *res_ax, uint64_t *res_bx, uint64_t *d2Tilde,
                           long type, map<long, Key> &keyMap_,
                           long expand_length, long curr_limbs, long beta, Scheme scheme) {

    uint64_t *sumaxmult = new uint64_t[expand_length]();
    uint64_t *sumbxmult = new uint64_t[expand_length]();

    //4. InnerProduct
    EvalFastKeySwitchCoreExt(sumaxmult, sumbxmult, d2Tilde, type, keyMap_, expand_length, beta, curr_limbs, scheme);

    //5. Modulus Down
    //fixme: 比较与
    // scheme.context.back(res_ax, sumaxmult, curr_limbs);
    // scheme.context.back(res_bx, sumbxmult, curr_limbs);实现效率的差距
    KeySwitchDown(res_ax, res_bx, sumaxmult, sumbxmult,
                  curr_limbs, scheme);
    delete[] sumaxmult;
    delete[] sumbxmult;
}


void KeySwitchCore(uint64_t *res_ax, uint64_t *res_bx, uint64_t *axax,
                   long type, map<long, Key> &keyMap_,
                   long curr_limbs, Scheme scheme) {
    uint32_t alpha = scheme.context.alpha;
    uint32_t logN = scheme.context.logN;
    uint32_t K = scheme.context.K;

    //precompute KeySwitchCore params
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //preallocate KeySwitchCore memory
    //total beta groups,
    // each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys
    // the size of each poly is 1<<context.logN
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *d2Tilde = new uint64_t[beta * expand_length];

    EvalKeySwitchPrecomputeCore(d2Tilde, axax, curr_limbs, beta, scheme);
    EvalFastKeySwitchCore(res_ax, res_bx, d2Tilde, type, keyMap_,
                          expand_length, curr_limbs, beta, scheme);
    delete[] d2Tilde;//fixme: 事实上，在完成 EvalFastKeySwitchCore 中的 InnerP 之后即可删除该变量，而无需等到完整的 KS_core 结束。此处目的是在new变量处，进行delete。

}

void partialNTTAndEqual(uint64_t* a,
                                 long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_len,
                                 long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_len, Scheme scheme)
{
    uint32_t logN=scheme.context.logN;

    for (long qi_modulus_index = 0,poly_index = 0; qi_modulus_index < out_C_L_len; ++qi_modulus_index,++poly_index) {
        if ((qi_modulus_index < in_C_L_index) || (qi_modulus_index >= in_C_L_index + in_C_L_len)) {
            uint64_t* ai = a + (poly_index << logN);
            scheme.context.qiNTTAndEqual(ai, qi_modulus_index);
        }
    }

    uint64_t modulus_index = 0;
    for (long poly_index = out_B_Start; poly_index < out_B_Start + out_B_len; ++modulus_index, ++poly_index) {
        uint64_t* ai = a + (poly_index << logN);
        scheme.context.piNTTAndEqual(ai, modulus_index);
    }
}

void ApproxModUp(uint64_t *ra, uint64_t *a,
                 long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_Len,
                 long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_Len, Scheme scheme) {
    uint32_t logN = scheme.context.logN;

    copy(a, a + (in_C_L_len << logN), ra + (in_C_L_index << logN));

    scheme.context.INTTAndEqual(a, in_C_L_len);

    BasisConversion(ra, a,
                    in_C_L_index, in_C_L_len, in_B_Start, in_B_Len,
                    out_C_L_index, out_C_L_len, out_B_Start, out_B_Len, scheme);

    partialNTTAndEqual(ra,
                       in_C_L_index, in_C_L_len, in_B_Start, in_B_Len,
                       out_C_L_index, out_C_L_len, out_B_Start, out_B_Len,scheme);

}

//todo: yx: FAB Algo1 好像有新算法
void ApproxModDown(uint64_t *ra, uint64_t *a, long curr_limbs,
                   long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_len, long out_C_L_index,
                   long out_C_L_len, long out_B_Start, long out_B_len, Scheme scheme) {
    uint32_t logN = scheme.context.logN;
    uint32_t N = scheme.context.N;
    uint64_t *qVec = scheme.context.qVec;
    uint64_t *QHatInvModqj = scheme.context.QHatInvModqj;
    uint64_t *PInvModq = scheme.context.PInvModq;

    scheme.context.INTTAndEqual(a, curr_limbs, in_B_len);

    uint64_t *tmp3 = new uint64_t[out_C_L_len << logN];
    copy(a, a + (out_C_L_len << logN), tmp3);

    BasisConversion(ra, a,
                    in_C_L_index, in_C_L_len, in_B_Start, in_B_len,
                    out_C_L_index, out_C_L_len, out_B_Start, out_B_len,scheme);

    //fixme: yx: move to precompute, 删除delete
    QHatInvModqj=new uint64_t [out_C_L_len]();
    for (int modulus_index = 0; modulus_index < out_C_L_len ; ++modulus_index) {
        uint64_t tmp = invMod(1, qVec[modulus_index]);
        mulMod(QHatInvModqj[modulus_index],tmp,PInvModq[modulus_index],qVec[modulus_index]);
    }

    //operate on RNS polys
    for (int j = 0; j < out_C_L_len; ++j) {
        uint64_t *rai = ra + (j << logN);
        uint64_t *tmp3i = tmp3 + (j << logN);
        for (int i = 0; i < N; ++i) {
            uint64_t rs = 0;
            subMod(rs, tmp3i[i], rai[i], qVec[j]); //Over 100x Algo4 ModDown Line7 subtraction
            mulMod(rai[i], rs, QHatInvModqj[j], qVec[j]); ////Over 100x Algo4 ModDown Line7 mult
        }
    }

    scheme.context.NTTAndEqual(ra, out_C_L_len);
    scheme.context.NTTAndEqual(a, curr_limbs, in_B_len);  // 还原输入a的样子

    delete[] tmp3;
    delete[] QHatInvModqj;
}

void BasisConversion(
        uint64_t* ra, uint64_t* a,
        long in_C_L_index, long in_C_L_len, long in_B_Start, long in_B_len,
        long out_C_L_index, long out_C_L_len, long out_B_Start, long out_B_len, Scheme scheme)
{

    uint32_t logN=scheme.context.logN;
    uint32_t N=scheme.context.N;
    uint64_t* qVec=scheme.context.qVec;
    uint64_t* qrVec=scheme.context.qrVec;
    uint64_t* pVec=scheme.context.pVec;
    uint64_t* prVec=scheme.context.prVec;
    long* pTwok=scheme.context.pTwok;
    uint64_t* PInvModq=scheme.context.PInvModq;
    uint64_t* PModq=scheme.context.PModq;
    long* qTwok=scheme.context.qTwok;
    uint64_t *QjHatDPrimeInvModqj = scheme.context.QjHatDPrimeInvModqj;
    uint64_t **QjHatDprimeModqi = scheme.context.QjHatDprimeModqi;
    uint64_t **QjHatDprimeModpi = scheme.context.QjHatDprimeModpi;
    uint64_t *QjHatTPrimeInvModpj = scheme.context.QjHatTPrimeInvModpj;
    uint64_t **QjHatTprimeModqi = scheme.context.QjHatTprimeModqi;

    //计算 Over100x Algo 1

    long in_C_L_end = in_C_L_index + in_C_L_len;
    long out_C_L_end = out_C_L_index + out_C_L_len;

    if (in_C_L_len > 0) {
        //fixme：yx: 相关辅助数组move to precompute，
        // 并删除最后的delete
        //fixme: 挪至预计算，需要设计通过index访问的方式
        QjHatDPrimeInvModqj = new uint64_t[in_C_L_end]();
        long modulus_index = in_C_L_index;
        for (; modulus_index < in_C_L_end; ++modulus_index) {
            uint64_t tmp = 1;
            for (int op_index = in_C_L_index; op_index < in_C_L_end; ++op_index) {
                if (op_index != modulus_index) {
                    mulModBarrett(tmp, tmp, qVec[op_index], qVec[modulus_index],
                                  qrVec[modulus_index], qTwok[modulus_index]);
                }
            }
            tmp = invMod(tmp, qVec[modulus_index]);
            if (in_B_len > 0) {
                mulMod(QjHatDPrimeInvModqj[modulus_index], tmp, PInvModq[modulus_index], qVec[modulus_index]);
            } else {
                QjHatDPrimeInvModqj[modulus_index] = tmp;
            }
        }

        //fixme: 挪至预计算，需要设计通过index访问的方式
        QjHatDprimeModqi = new uint64_t *[in_C_L_end]; //traverse Qj
        QjHatDprimeModpi = new uint64_t *[in_C_L_end]; //traverse Qj
        for (int i = 0; i < in_C_L_end; ++i) {
            QjHatDprimeModqi[i] = new uint64_t[out_C_L_len]();    //traverse qi
            QjHatDprimeModpi[i] = new uint64_t[out_B_len]();    //traverse pi
        }
        //out_C_L_index default is 0
        for (long qi_modulus_index = 0; qi_modulus_index <
                                        out_C_L_end; ++qi_modulus_index) { //qi_modulus_index=>每一个 q_i in L1 of Algo1 in Over100x
            if ((qi_modulus_index < in_C_L_index) ||
                (qi_modulus_index >= in_C_L_index + in_C_L_len)) { //QjHatDprimeModqi在 qi_modulus_index 的 这些取值下有意义
                for (long Qj_index = in_C_L_index; Qj_index < in_C_L_end; ++Qj_index) {
                    uint64_t multModqi = 1;
                    for (int op_index = in_C_L_index; op_index < in_C_L_index + in_C_L_len; ++op_index) {
                        if (op_index != Qj_index) {
                            mulModBarrett(multModqi, multModqi, qVec[op_index], qVec[qi_modulus_index],
                                          qrVec[qi_modulus_index], qTwok[qi_modulus_index]);
                        }
                    }
                    if (in_B_len > 0) {
                        mulModBarrett(QjHatDprimeModqi[Qj_index][qi_modulus_index], multModqi, PModq[qi_modulus_index],
                                      qVec[qi_modulus_index], qrVec[qi_modulus_index], qTwok[qi_modulus_index]);
                    } else {
                        QjHatDprimeModqi[Qj_index][qi_modulus_index] = multModqi;
                    }
                }
            }
        }

        //out_B_len>0, 是 ModUp 的 case，因此一定有 in_B_len=0, 所以计算QjHatDprimeModpi 最下面不需要像前面一样判断 if (in_B_len > 0)
        for (long pi_modulus_index = 0;
             pi_modulus_index < out_B_len; ++pi_modulus_index) { //pi_modulus_index=>每一个 q_i in L1 of Algo1 in Over100x
            for (long Qj_index = in_C_L_index; Qj_index < in_C_L_end; ++Qj_index) {
                uint64_t multModpi = 1;
                for (int op_index = in_C_L_index; op_index < in_C_L_index + in_C_L_len; ++op_index) {
                    if (op_index != Qj_index) {
                        uint64_t temp = qVec[op_index] % pVec[pi_modulus_index];
                        mulModBarrett(multModpi, multModpi, temp, pVec[pi_modulus_index],
                                      prVec[pi_modulus_index], pTwok[pi_modulus_index]);
                    }
                }
                QjHatDprimeModpi[Qj_index][pi_modulus_index] = multModpi;
            }
        }

        //todo: yx: a(x)mod q_j 或者 a(x)mod p_j first before mulModBarrett
        // (Over100x Algo 1, step2\3)

        //line2 待求和项的左边第一项乘积
        uint64_t *tmp3 = new uint64_t[in_C_L_len << logN];
        modulus_index = in_C_L_index;
        for (long poly_index = 0; poly_index < in_C_L_len; ++poly_index, ++modulus_index) {
            uint64_t *tmp3i = tmp3 + (poly_index << logN);
            uint64_t *ai = a + (poly_index << logN);
            for (long n = 0; n < N; ++n) {
                mulModBarrett(tmp3i[n], ai[n], QjHatDPrimeInvModqj[modulus_index],
                              qVec[modulus_index], qrVec[modulus_index], qTwok[modulus_index]);
            }
        }

        ////////
        //计算 RNS poly 模 q_0 -> q_{in_C_L_index-1} 的部分
        //计算 RNS poly 模 q_{in_C_L_end} -> q_alppha*beta 的部分
        modulus_index = 0;
        //out_C_L_index default is 0
        for (long poly_index = 0;
             poly_index < out_C_L_end; ++poly_index, ++modulus_index) { //poly_index=>每一个 q_i in L1 of Algo1 in Over100x
            if ((poly_index < in_C_L_index) || (poly_index >= in_C_L_index + in_C_L_len)) {
                uint64_t *rak = ra + (poly_index << logN);
                for (long n = 0; n < N; ++n) {
                    uint64_t tt;
                    unsigned __int128 sum = static_cast<unsigned __int128>(0);
                    for (long i = 0; i < in_C_L_len; ++i) { //i=>j in Algo1 in Over100x
                        tt = tmp3[n + (i << logN)];
                        sum += static_cast<unsigned __int128>(tt) * QjHatDprimeModqi[i][modulus_index];
                    }
                    modBarrett(rak[n], sum, qVec[modulus_index], qrVec[modulus_index], qTwok[modulus_index]);
                }
            }

        }

        //计算 RNS poly 模 P=Πp_i 的部分
        modulus_index = 0;
        for (long poly_index = out_B_Start; poly_index < out_B_Start + out_B_len; ++poly_index, ++modulus_index) {
            uint64_t *rak = ra + (poly_index << logN);
            for (long n = 0; n < N; ++n) {
                uint64_t tt;
                unsigned __int128 sum = static_cast<unsigned __int128>(0);
                for (long i = 0; i < in_C_L_len; ++i) {
                    tt = tmp3[n + (i << logN)];
                    sum += static_cast<unsigned __int128>(tt) * QjHatDprimeModpi[i][modulus_index];
                }
                modBarrett(rak[n], sum, pVec[modulus_index], prVec[modulus_index], pTwok[modulus_index]);
            }
        }
        delete[] tmp3;
        delete[] QjHatDPrimeInvModqj;
        for (int i = 0; i < in_C_L_len; ++i) {
            delete[] QjHatDprimeModqi[i];
            delete[] QjHatDprimeModpi[i];
        }
        delete[] QjHatDprimeModqi;
        delete[] QjHatDprimeModpi;
    }


    if (in_B_len>0) {
        //fixme: yx: 挪到预计算里去，并删除最后的 delete[] xxx;
        QjHatTPrimeInvModpj = new uint64_t [in_B_len]();
        long modulus_index = 0;
        for (; modulus_index < in_B_len; ++modulus_index) {
            uint64_t mult_pi_Inv=1;
            for (int op_index = 0; op_index < in_B_len; ++op_index) {
                if (op_index!=modulus_index) {
                    mulModBarrett(mult_pi_Inv, mult_pi_Inv, pVec[op_index], pVec[modulus_index],
                                  prVec[modulus_index], pTwok[modulus_index]);
                }
            }
            mult_pi_Inv=invMod(mult_pi_Inv, pVec[modulus_index]);

            uint64_t mult_qi_Inv=1;
            for (int op_index = in_C_L_index; op_index < in_C_L_end; ++op_index) {
                mulModBarrett(mult_qi_Inv, mult_qi_Inv, qVec[op_index], pVec[modulus_index],
                              prVec[modulus_index], pTwok[modulus_index]);
            }
            mult_qi_Inv=invMod(mult_qi_Inv, pVec[modulus_index]);

            mulMod(QjHatTPrimeInvModpj[modulus_index], mult_pi_Inv, mult_qi_Inv, pVec[modulus_index]);
        }

        QjHatTprimeModqi = new uint64_t* [in_B_len]; //traverse Qj
        for (int i = 0; i < in_B_len; ++i) {
            QjHatTprimeModqi[i]= new uint64_t [out_C_L_len]();    //traverse pi
        }
        for (long qi_modulus_index = 0; qi_modulus_index < out_C_L_len; ++qi_modulus_index) { //qi_modulus_index=>每一个 q_i in L1 of Algo1 in Over100x
            if ((qi_modulus_index < in_C_L_index)||(qi_modulus_index >= in_C_L_index + in_C_L_len)) { //QjHatTprimeModqi  qi_modulus_index 的 这些取值下有意义
                modulus_index= 0 ;
                for (int Qj_index=0; Qj_index < in_B_len; ++Qj_index) {
                    uint64_t mult_pi=1;
                    for (int op_index = 0; op_index < in_B_len; ++op_index) {
                        if (op_index!=Qj_index) {
                            mulModBarrett(mult_pi, mult_pi, pVec[op_index], qVec[qi_modulus_index],
                                          qrVec[qi_modulus_index], qTwok[qi_modulus_index]);
                        }
                    }

                    uint64_t mult_qi=1;
                    for (int op_index = in_C_L_index; op_index < in_C_L_end; ++op_index) {

                        mulModBarrett(mult_qi, mult_qi, qVec[op_index], qVec[qi_modulus_index],
                                      qrVec[qi_modulus_index], qTwok[qi_modulus_index]);
                    }
                    mulMod(QjHatTprimeModqi[Qj_index][qi_modulus_index], mult_pi, mult_qi, qVec[qi_modulus_index]);
                }
            }
        }

        //line3 待求和项的左边第一项乘积
        uint64_t *tmp4 = new uint64_t[in_B_len << logN];
        modulus_index = 0;
        for (long poly_index = 0; poly_index < in_B_len; ++poly_index, ++modulus_index) {
            uint64_t *tmp4i = tmp4 + (poly_index << logN);
            uint64_t *ai = a + ((in_B_Start + poly_index) << logN);
            for (long n = 0; n < N; ++n) {
                mulModBarrett(tmp4i[n], ai[n], QjHatTPrimeInvModpj[modulus_index],
                              pVec[modulus_index], prVec[modulus_index], pTwok[modulus_index]);
            }
        }

        ////////
        //计算 RNS poly 模 q_0 -> q_{in_C_L_index-1} 的部分
        //计算 RNS poly 模 q_{in_C_L_end} -> q_beta 的部分
        //和上面的区别，最内层 for-loop 的循环次数，QjHatTprimeModqi
        modulus_index=0;
        for (long poly_index = 0; poly_index < out_C_L_end; ++poly_index,++modulus_index) { //poly_index=>每一个 q_i in L1 of Algo1 in Over100x
            if ((poly_index < in_C_L_index)||(poly_index >= in_C_L_index + in_C_L_len)) {
                uint64_t* rak = ra + (poly_index << logN);
                for (long n = 0; n < N; ++n) {
                    uint64_t tt;
                    unsigned __int128 sum = static_cast<unsigned __int128>(0);
                    for (long i = 0; i < in_B_len; ++i) { //i=>j in Algo1 in Over100x
                        tt  = tmp4[n + (i << logN)];
                        sum+= static_cast<unsigned __int128>(tt) * QjHatTprimeModqi[i][modulus_index];
                    }
                    uint64_t a_i_Mod_q_i;
                    modBarrett(a_i_Mod_q_i, sum, qVec[modulus_index], qrVec[modulus_index], qTwok[modulus_index]);
                    addMod(rak[n], rak[n], a_i_Mod_q_i, qVec[modulus_index]);// Algo1 Line4
                }
            }
        }

        delete [] tmp4;
        delete [] QjHatTPrimeInvModpj;
        for (int i = 0; i < in_B_len; ++i) {
            delete [] QjHatTprimeModqi[i];    //traverse pi
        }
        delete [] QjHatTprimeModqi; //traverse Qj
    }

}
