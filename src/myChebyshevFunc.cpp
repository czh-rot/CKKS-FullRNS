//
// Created by EYx on 2023/4/2.
//

#include "myChebyshevFunc.h"
//
// Created by EYx on 2023/4/2.
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

//my implementation
#include "myMult.h"

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>

#include "../data/testConstValue.h"
#include "ckksrns-leveledshe.h"

using namespace std;
using namespace chrono;

void EvalChebyshevSeries(Ciphertext& result, Ciphertext x, double* coefficients,
                               double a, double b, uint32_t poly_degree,
                               Scheme scheme)
{

    uint32_t n = Degree(coefficients,poly_degree);

//    if (n < 5) {
//        return EvalChebyshevSeriesLinear(x, coefficients, a, b);
//    }

//    return EvalChebyshevSeriesPS(x, coefficients, a, b, poly_degree, scheme);
    EvalChebyshevSeriesPS(result, x, coefficients, a, b, poly_degree, scheme);
}

void EvalChebyshevSeriesPS(Ciphertext& result, Ciphertext x, double* coefficients,
                                 double a, double b, uint32_t coefficients_len,
                                 Scheme scheme)
{
    Ciphertext tmpct1,tmpct2;//buffer，用于存放调整level后的密文
    uint32_t n = Degree(coefficients, coefficients_len);

    //std::vector<double> f2 = coefficients;
    uint32_t f2_len=coefficients_len;
    double *f2=new double[f2_len];
    for (int i = 0; i < f2_len; ++i) {
        f2[i]=coefficients[i];
    }

    // Make sure the coefficients do not have the zero dominant terms
    if(coefficients[coefficients_len - 1] == 0) //    if (coefficients[coefficients.size() - 1] == 0)
        f2_len=n+1; //f2.resize(n + 1);

    std::vector<uint32_t> degs = ComputeDegreesPS(n);
    uint32_t k                 = degs[0];
    uint32_t m                 = degs[1];

    // computes linear transformation y = -1 + 2 (x-a)/(b-a)
    // consumes one level when a <> -1 && b <> 1
//    auto cc = x->GetCryptoContext();
    Ciphertext *T=new Ciphertext[k]; ////    std::vector<Ciphertext<DCRTPoly>> T(k);
    if ((a - std::round(a) < 1e-10) && (b - std::round(b) < 1e-10) && (std::round(a) == -1) && (std::round(b) == 1)) {
        // no linear transformation is needed if a = -1, b = 1
        // T_1(y) = y
        T[0]=x;//T[0] = x->Clone();
    }
    else {
        // linear transformation is needed
        double alpha = 2 / (b - a);
        double beta  = 2 * a / (b - a);

//        T[0]=scheme.multByConst(x,alpha);//T[0] = cc->EvalMult(x, alpha);
        T[0]=EvalMult(x,alpha,scheme);//T[0] = cc->EvalMult(x, alpha);
//        T[0]=scheme.Lattigo_MultByConst(x,alpha);//T[0] = cc->EvalMult(x, alpha);
//        scheme.reScaleByAndEqual(T[0], 1); //cc->ModReduceInPlace(T[0]);
        ModReduceInternalInPlace(T[0],1,scheme);

        scheme.addConstAndEqual(T[0], -1.0 - beta); // cc->EvalAddInPlace(T[0], -1.0 - beta);

    }

    Ciphertext y = T[0];

    // Computes Chebyshev polynomials up to degree k
    // for y: T_1(y) = y, T_2(y), ... , T_k(y)
    // uses binary tree multiplication
    for (uint32_t i = 2; i <= k; i++) {
        // if i is a power of two
        if (!(i & (i - 1))) {
            // compute T_{2i}(y) = 2*T_i(y)^2 - 1
            Ciphertext square(T[i / 2 - 1].N,T[i / 2 - 1].slots,T[i / 2 - 1].curr_limbs);
            mySquare(square,T[i / 2 - 1],scheme);
            T[i - 1]    = scheme.add(square, square); //T[i - 1]    = cc->EvalAdd(square, square);
            //fixme: yx: high-precision rescale technique可能需要改变下两行的执行顺序
//            scheme.reScaleByAndEqual(T[i - 1], 1);//cc->ModReduceInPlace(T[i - 1]);
            ModReduceInternalInPlace(T[i - 1],1,scheme);
            scheme.addConstAndEqual(T[i-1],-1.0);//cc->EvalAddInPlace(T[i - 1], -1.0);
        }
        else {
            // non-power of 2
            if (i % 2 == 1) {
                // if i is odd
                // compute T_{2i+1}(y) = 2*T_i(y)*T_{i+1}(y) - y
                CheckAndAdjusLevel(tmpct1,tmpct2,T[i / 2 - 1], T[i / 2],scheme);
                Ciphertext prod(tmpct1.N,tmpct1.slots,tmpct1.curr_limbs);
                myMult(prod,tmpct1,tmpct2,scheme);//auto prod = cc->EvalMult(T[i / 2 - 1], T[i / 2]);//auto prod = cc->EvalMult(T[i / 2 - 1], T[i / 2]);
                T[i-1]=scheme.add(prod,prod);//T[i - 1]  = cc->EvalAdd(prod, prod);
//                scheme.reScaleByAndEqual(T[i - 1], 1);//cc->ModReduceInPlace(T[i - 1]);
                ModReduceInternalInPlace(T[i-1],1,scheme);
                CheckAndAdjusLevel(T[i - 1],tmpct2,T[i - 1], y,scheme);
                scheme.subAndEqual(T[i - 1], tmpct2);//cc->EvalSubInPlace(T[i - 1], y);
            }
            else {
                // i is even but not power of 2
                // compute T_{2i}(y) = 2*T_i(y)^2 - 1
                Ciphertext square(T[i / 2 - 1].N,T[i / 2 - 1].slots,T[i / 2 - 1].curr_limbs);
                mySquare(square,T[i / 2 - 1],scheme);//auto square = cc->EvalSquare(T[i / 2 - 1]);
                T[i - 1] = scheme.add(square, square);//T[i - 1] = cc->EvalAdd(square, square);
//                scheme.reScaleByAndEqual(T[i - 1], 1);//cc->ModReduceInPlace(T[i - 1]);
                ModReduceInternalInPlace(T[i-1],1,scheme);
                scheme.addConstAndEqual(T[i-1],-1.0);//cc->EvalAddInPlace(T[i - 1], -1.0);
            }
        }
    }
    //fixme: yx: 现在默认 rescaleTech 是 FIXEDMANUAL，之后支持 FIXEDAUTO 模式后这里也要改
//    if (cryptoParams->GetScalingTechnique() == FIXEDMANUAL) {
        // brings all powers of x to the same level
        for (size_t i = 1; i < k; i++) {
            uint64_t levelDiff = T[i - 1].curr_limbs - T[k - 1].curr_limbs ;//usint levelDiff = T[k - 1]->GetLevel() - T[i - 1]->GetLevel();
            scheme.modDownByAndEqual(T[i-1],levelDiff);//cc->LevelReduceInPlace(T[i - 1], nullptr, levelDiff);
        }
//    }
//    else{
//        for (size_t i = 1; i < k; i++) {
//          // call AdjustLevelsAndDepthInPlace，这个函数我怎么实现
//          algo->AdjustLevelsAndDepthInPlace(T[i - 1], T[k - 1]);
//        }
//    }

    Ciphertext *T2=new Ciphertext[m]; //std::vector<Ciphertext<DCRTPoly>> T2(m);
    // Compute the Chebyshev polynomials T_{2k}(y), T_{4k}(y), ... , T_{2^{m-1}k}(y)
    T2[0]=T[k-1];//T2.front() = T.back();
    for (uint32_t i = 1; i < m; i++) {
        Ciphertext square(T2[i - 1].N,T2[i - 1].slots,T2[i - 1].curr_limbs);
        mySquare(square,T2[i - 1],scheme); //auto square = cc->EvalSquare(T2[i - 1]);
        T2[i]=scheme.add(square, square);           //T2[i] = cc->EvalAdd(square, square);
//        scheme.reScaleByAndEqual(T2[i], 1);         //cc->ModReduceInPlace(T2[i]);
        ModReduceInternalInPlace(T2[i],1,scheme);
        scheme.addConstAndEqual(T2[i],-1.0);     //cc->EvalAddInPlace(T2[i], -1.0);
    }

    // computes T_{k(2*m - 1)}(y)
    Ciphertext T2km1 = T2[0];   //auto T2km1 = T2.front();
    for (uint32_t i = 1; i < m; i++) {
        // compute T_{k(2*m - 1)} = 2*T_{k(2^{m-1}-1)}(y)*T_{k*2^{m-1}}(y) - T_k(y)
        CheckAndAdjusLevel(tmpct1,tmpct2,T2km1, T2[i],scheme);
        Ciphertext prod(tmpct1.N,tmpct1.slots,tmpct1.curr_limbs);
        myMult(prod,tmpct1,tmpct2,scheme);//auto prod = cc->EvalMult(T2km1, T2[i]);
        T2km1=scheme.add(prod,prod);               //T2km1 = cc->EvalAdd(prod, prod);
//        scheme.reScaleByAndEqual(T2km1, 1);        //cc->ModReduceInPlace(T2km1);
        ModReduceInternalInPlace(T2km1,1,scheme);
        CheckAndAdjusLevel(T2km1,tmpct2,T2km1, T2[0],scheme);
        scheme.subAndEqual(T2km1,tmpct2);//cc->EvalSubInPlace(T2km1, T2.front());
    }

    // We also need to reduce the number of levels of T[k-1] and of T2[0] by another level.
    //  cc->LevelReduceInPlace(T[k-1], nullptr);
    //  cc->LevelReduceInPlace(T2.front(), nullptr);

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Add T^{k(2^m - 1)}(y) to the polynomial that has to be evaluated
    //f2.resize(2 * k2m2k + k + 1, 0.0);
    uint32_t new_f2_len=  2 * k2m2k + k + 1;
    if (f2_len<new_f2_len){
//        cout<<"error!! 应该"
//              "1）在最初计算出 k2m2k；"
//              "2）初始时就分配更大的f2；"
//              "3）初始时用不到的空间应该置零！"<<endl;
        double *new_f2=new double[new_f2_len];
        for (int i = 0; i < f2_len; ++i) {
            new_f2[i]=f2[i];
        }
        delete[] f2; //！！！注意：这个delete不能随意移动！
        for (uint32_t i = f2_len; i < new_f2_len-1; ++i) {
            new_f2[i]=0.0;
        }
        new_f2[new_f2_len - 1]=1;//f2.back() = 1;
        f2_len=new_f2_len;
        f2=new_f2;
    } else{//f2_len>=new_f2_len
        f2_len=new_f2_len;
        f2[f2_len-1]=1;
    }

    // Divide f2 by T^{k*2^{m-1}}
    long Tkm_len=int32_t(k2m2k + k) + 1;
    double* Tkm=new double[Tkm_len]; //std::vector<double> Tkm(int32_t(k2m2k + k) + 1, 0.0);
    for (int i = 0; i < Tkm_len; ++i) {
        Tkm[i]=0.0;
    }
    Tkm[Tkm_len-1] = 1;//Tkm.back() = 1;
//    auto divqr = LongDivisionChebyshev(f2, Tkm);
    double *divqr_q, *divqr_r;
    uint32_t divqr_q_len, divqr_r_len;
    LongDivisionChebyshev(divqr_q, divqr_q_len, divqr_r, divqr_r_len, f2, f2_len, Tkm, Tkm_len);
    delete []Tkm;
    delete []f2;

    // Subtract x^{k(2^{m-1} - 1)} from r
//    std::vector<double> r2 = divqr->r;
//    if (int32_t(k2m2k - Degree(divqr->r)) <= 0) {
//        r2[int32_t(k2m2k)] -= 1;
//        r2.resize(Degree(r2) + 1);
//    }
//    else {
//        r2.resize(int32_t(k2m2k + 1), 0.0);
//        r2.back() = -1;
//    }

    uint32_t r2_len;
    double *r2;
    if (int32_t(k2m2k - Degree(divqr_r,divqr_r_len)) <= 0){
        divqr_r[int32_t(k2m2k)]-=1;
        r2_len = Degree(divqr_r,divqr_r_len)+1;
        r2 = new double [r2_len];
        for (int i = 0;(i < divqr_r_len) && (i<r2_len); ++i) {
            r2[i]=divqr_r[i];
        }
        divqr_r[int32_t(k2m2k)]+=1; //处理 r2[int32_t(k2m2k)] -= 1;
    }
    else{
        r2_len=int32_t(k2m2k + 1);
        r2 =new double [r2_len];
        for (int i = 0; (i < divqr_r_len) && (i<r2_len); ++i) {
            r2[i]=divqr_r[i];
        }
        for (int i = divqr_r_len; i < r2_len-1; ++i) {
            r2[i]=0.0;
        }
        r2[r2_len-1]=-1;
    }

    // Divide r2 by q
    //auto divcs = LongDivisionChebyshev(r2, divqr->q);
    double *divcs_q, *divcs_r;
    uint32_t divcs_q_len, divcs_r_len;
    LongDivisionChebyshev(divcs_q, divcs_q_len, divcs_r, divcs_r_len, r2, r2_len, divqr_q, divqr_q_len);

    delete[] r2;

    // Add x^{k(2^{m-1} - 1)} to s
    double *s2;//std::vector<double> s2 = divcs->r;
    uint32_t s2_len=int32_t(k2m2k + 1);//s2.resize(int32_t(k2m2k + 1), 0.0);
    s2=new double[s2_len];
    for (int i = 0; (i <divcs_r_len)&&(i<s2_len); ++i) {
        s2[i]=divcs_r[i];
    }
    for (uint32_t i = divcs_r_len; i <s2_len-1; ++i) {
        s2[i]=0.0;
    }
    s2[s2_len-1]=1;//s2.back() = 1;

    // Evaluate c at u
    Ciphertext cu;//Ciphertext<DCRTPoly> cu;
    uint32_t dc = Degree(divcs_q,divcs_q_len);//uint32_t dc = Degree(divcs->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            if (divcs_q[1] != 1) {//if (divcs->q[1] != 1) {
//                cu = scheme.multByConst(T[0], divcs_q[1]);//cu = cc->EvalMult(T.front(), divcs->q[1]);
                cu = EvalMult(T[0], divcs_q[1],scheme);//cu = cc->EvalMult(T.front(), divcs->q[1]);
//                cu = scheme.Lattigo_MultByConst(T[0], divcs_q[1]);//cu = cc->EvalMult(T.front(), divcs->q[1]);
//                scheme.reScaleByAndEqual(cu, 1);//cc->ModReduceInPlace(cu);
                ModReduceInternalInPlace(cu,1,scheme);
            } else {
                cu = T[0];//cu = T.front();
            }
        }
        else {
            Ciphertext *ctxs = new Ciphertext[dc](); //std::vector<Ciphertext<DCRTPoly>> ctxs(dc);
            double *weights = new double[dc](); //std::vector<double> weights(dc);

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i] = T[i];
                weights[i] = divcs_q[i + 1];//weights[i] = divcs->q[i + 1];
            }
            //cu = cc->EvalLinearWSumMutable(ctxs, weights);
            EvalLinearWSumMutable(cu, ctxs, dc, weights, scheme);
            delete []ctxs;
            delete []weights;

        }

        // adds the free term (at x^0)
        scheme.addConstAndEqual(cu, divcs_q[0] / 2);//cc->EvalAddInPlace(cu, divcs->q.front() / 2);
        // TODO : Andrey why not T2[m-1]->GetLevel() instead?
        // Need to reduce levels to the level of T2[m-1].
        //    usint levelDiff = y->GetLevel() - cu->GetLevel() + ceil(log2(k)) + m - 1;
        //    cc->LevelReduceInPlace(cu, nullptr, levelDiff);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    Ciphertext qu;//Ciphertext<DCRTPoly> qu;

    if (Degree(divqr_q,divqr_q_len) > k) {//if (Degree(divqr->q) > k) {
        // InnerEvalChebyshevPS待确认作用和意义
//        qu = InnerEvalChebyshevPS(x, divqr->q, k, m - 1, T, T2);
        InnerEvalChebyshevPS(qu, x, divqr_q, divqr_q_len,k, m - 1, T, T2,scheme);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        double *qcopy; //auto qcopy = divqr->q;
        uint32_t qcopy_len=k; //qcopy.resize(k);
        qcopy=new double [qcopy_len]();
        for (int i = 0; (i < divqr_q_len)&&(i < qcopy_len); ++i) {
            qcopy[i]=divqr_q[i];
        }
        for (uint32_t i = divqr_q_len; i < qcopy_len; ++i) {
            qcopy[i]=0.0;
        }

        uint32_t deg_qcopy=Degree(qcopy,qcopy_len);
        if (deg_qcopy > 0) {//if (Degree(qcopy) > 0) {
            Ciphertext* ctxs=new Ciphertext[deg_qcopy]();//std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(qcopy));
            double* weights = new double[deg_qcopy](); //std::vector<double> weights(Degree(qcopy));

            for (uint32_t i = 0; i < deg_qcopy; i++) {//for (uint32_t i = 0; i < Degree(qcopy); i++) {
                ctxs[i]    = T[i];
                weights[i] = divqr_q[i + 1];//weights[i] = divqr->q[i + 1];
            }

            //qu = cc->EvalLinearWSumMutable(ctxs, weights);
            EvalLinearWSumMutable(qu,ctxs,deg_qcopy,weights,scheme);
            // the highest order coefficient will always be 2 after one division because of the Chebyshev division rule
            Ciphertext sum= scheme.add(T[k - 1], T[k - 1]);//Ciphertext<DCRTPoly> sum = cc->EvalAdd(T[k - 1], T[k - 1]);
            CheckAndAdjusLevel(qu,sum,qu,sum,scheme);//对齐level
            scheme.addAndEqual(qu,sum);//cc->EvalAddInPlace(qu, sum);
        }
        else {
            qu = T[k - 1];

            for (uint32_t i = 1; i < divqr_q[divqr_q_len-1]; i++) {//for (uint32_t i = 1; i < divqr->q.back(); i++) {
                scheme.addAndEqual(qu, T[k - 1]);//cc->EvalAddInPlace(qu, T[k - 1]);
            }
        }

        // adds the free term (at x^0)
        scheme.addConstAndEqual(qu, divqr_q[0] / 2);//cc->EvalAddInPlace(qu, divqr->q.front() / 2);
        // The number of levels of qu is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so the number of levels of qu and T2[m-1] will be the same.
        delete []qcopy;
    }

    Ciphertext su;//Ciphertext<DCRTPoly> su;
    uint32_t deg_s2= Degree(s2,s2_len);
    if (deg_s2 > k) {//if (Degree(s2) > k) {
//        su = InnerEvalChebyshevPS(x, s2, k, m - 1, T, T2);
        InnerEvalChebyshevPS(su, x, s2,s2_len, k, m - 1, T, T2,scheme);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        double *scopy; //auto scopy = s2;
        uint32_t scopy_len=k; //scopy.resize(k);
        scopy=new double [scopy_len]();
        for (int i = 0; (i < s2_len)&&(i < scopy_len); ++i) {
            scopy[i]=s2[i];
        }
        for (uint32_t i = s2_len; i < scopy_len; ++i) {
            scopy[i]=0.0;
        }

        uint32_t deg_scopy=Degree(scopy,scopy_len);
        if (deg_scopy > 0) {
            Ciphertext * ctxs=new Ciphertext[deg_scopy](); //std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(scopy));
            double * weights=new double [deg_scopy](); //std::vector<double> weights(Degree(scopy));

            for (uint32_t i = 0; i < deg_scopy; i++) {//for (uint32_t i = 0; i < Degree(scopy); i++) {
                ctxs[i]    = T[i];
                weights[i] = s2[i + 1];
            }

            // su = cc->EvalLinearWSumMutable(ctxs, weights);
            EvalLinearWSumMutable(su,ctxs,deg_scopy,weights,scheme);

            // the highest order coefficient will always be 1 because s2 is monic.
            //cc->EvalAddInPlace(su, T[k - 1]);
            Ciphertext tmp_T;
            if (T[k-1].curr_limbs > su.curr_limbs){
                CheckAndAdjusLevel(su,tmp_T,su,T[k-1],scheme);
                scheme.addAndEqual(su, tmp_T);
            }
            else{
                CheckAndAdjusLevel(su,T[k-1],su,T[k-1],scheme);
                scheme.addAndEqual(su, T[k - 1]);
            }
        }
        else {
            su = T[k - 1];
        }

        // adds the free term (at x^0)
        scheme.addConstAndEqual(su,s2[0]/2);//cc->EvalAddInPlace(su, s2.front() / 2);
        // The number of levels of su is the same as the number of levels of T[k-1] + 1.
        // Will only get here when m = 2, so need to reduce the number of levels by 1.

        delete []scopy;
    }

    // TODO : Andrey : here is different from 895 line
    // Reduce number of levels of su to number of levels of T2km1.
    //  cc->LevelReduceInPlace(su, nullptr);

//    Ciphertext result;//Ciphertext<DCRTPoly> result;

    if (flag_c) {
        CheckAndAdjusLevel(T2[m - 1],cu,T2[m - 1],cu,scheme);
        result=scheme.add(T2[m - 1],cu);//result = cc->EvalAdd(T2[m - 1], cu);
    }
    else {
        result =scheme.addConst(T2[m - 1], divcs_q[0] / 2);//result = cc->EvalAdd(T2[m - 1], divcs->q.front() / 2);
    }

    CheckAndAdjusLevel(result,qu,result,qu,scheme);
    myMultAndEqual(result,qu,scheme);//result = cc->EvalMult(result, qu);
//    scheme.reScaleByAndEqual(result,1);//cc->ModReduceInPlace(result);
    ModReduceInternalInPlace(result,1,scheme);

    CheckAndAdjusLevel(result,su,result,su,scheme);
    scheme.addAndEqual(result,su);//cc->EvalAddInPlace(result, su);
    CheckAndAdjusLevel(result,T2km1,result,T2km1,scheme);
    scheme.subAndEqual(result,T2km1);//cc->EvalSubInPlace(result, T2km1);

    delete [] T;
    delete [] T2;

//    return result;
}

/*Return the degree of the polynomial described by coefficients,
which is the index of the last non-zero element in the coefficients - 1.
Don't throw an error if all the coefficients are zero, but return 0. */
uint32_t Degree(double *coefficients, uint32_t poly_degree) {
    uint32_t deg = 1;
    for (uint32_t i = poly_degree - 1; i > 0; i--) {
        if (coefficients[i] == 0) {
            deg += 1;
        }
        else
            break;
    }
    return poly_degree - deg;
}

/* Compute positive integers k,m such that n < k(2^m-1) and k close to sqrt(n/2) */
//todo: yx: 将下述 vector 返回值改为 array？或者在特定密码学参数和应用参数条件下，提前计算出 (k,m) 并写死？
std::vector<uint32_t> ComputeDegreesPS(const uint32_t n) {

    std::vector<uint32_t> klist;
    std::vector<uint32_t> mlist;

    double sqn2 = sqrt(n / 2);

    for (uint32_t k = 1; k <= n; k++) {
        for (uint32_t m = 1; m <= ceil(log2(n / k) + 1) + 1; m++) {
            if (int32_t(n - k * ((1 << m) - 1)) < 0) {
                if ((static_cast<double>(k - sqn2) >= -2) && ((static_cast<double>(k - sqn2) <= 2))) {
                    klist.push_back(k);
                    mlist.push_back(m);
                }
            }
        }
    }

    uint32_t minIndex = std::min_element(mlist.begin(), mlist.end()) - mlist.begin();

    return std::vector<uint32_t>{{klist[minIndex], mlist[minIndex]}};
}

void CheckAndAdjusLevel(Ciphertext& rct1, Ciphertext& rct2, Ciphertext& ct1, Ciphertext& ct2,
                        Scheme scheme)
{
    rct1=ct1;
    rct2=ct2;
    //curr_limbs means current total limbs, should choose the lower one
    if (rct1.curr_limbs > rct2.curr_limbs){
        scheme.modDownToAndEqual(rct1,rct2.curr_limbs);
    }
    else if (rct1.curr_limbs < rct2.curr_limbs){
        scheme.modDownToAndEqual(rct2,rct1.curr_limbs);
    }
}
void LongDivisionChebyshev(double*& quotient,uint32_t & quotient_len,
                           double*& remainder,uint32_t & remainder_len,
                           double* f, uint32_t f_len,
                           double* g, uint32_t g_len
                           )
{
    uint32_t n = Degree(f,f_len);
    uint32_t k = Degree(g,g_len);

    if (n != f_len - 1) {//if (n != f.size() - 1) {
        cout<<"error LongDivisionChebyshev: The dominant coefficient of the divident is zero."<<endl;//OPENFHE_THROW(math_error, "LongDivisionChebyshev: The dominant coefficient of the divident is zero.");
    }

    if (k != g_len - 1) {//if (k != g.size() - 1) {
        cout<<"\"LongDivisionChebyshev: The dominant coefficient of the divisor is zero."<<endl;//OPENFHE_THROW(math_error, "LongDivisionChebyshev: The dominant coefficient of the divisor is zero.");
    }

    double *q;//std::vector<double> q;
    uint32_t q_len;
    uint32_t r_len=f_len;
    double *r=new double[r_len]; //std::vector<double> r = f;
    for (int i = 0; i < r_len; ++i) {
        r[i]=f[i];
    }

    if (int32_t(n - k) >= 0) {
        uint32_t q2_len= n-k+1;
        double *q2=new double[q2_len];//std::vector<double> q2(n - k + 1, 0.0);
        q_len=q2_len;
        q = new double[q_len];  //q = q2;
        for (int i = 0; i < q2_len; ++i) {
            q2[i]=0.0;
            q[i]=0.0;
        }

        while (int32_t(n - k) > 0) {
            q[n - k] = 2 * r[r_len-1];//q[n - k] = 2 * r.back();
            // 实现 IsNotEqualOne
            if (IsNotEqualOne(g[k])) {
                q[n - k] /= g[g_len-1];//q[n - k] /= g.back();
            }

            uint32_t d_len=n+1;
            double *d=new double [d_len]();//std::vector<double> d(n + 1, 0.0);
            for (int i = 0; i < d_len; ++i) {
                d[i]=0.0;
            }

            if (int32_t(k) == int32_t(n - k)) {
                d[0] = 2 * g[n - k];//d.front() = 2 * g[n - k];

                for (uint32_t i = 1; i < 2 * k + 1; i++) {
                    d[i] = g[abs(int32_t(n - k - i))];
                }
            }
            else {
                if (int32_t(k) > int32_t(n - k)) {
                    d[0] = 2 * g[n - k];//d.front() = 2 * g[n - k];
                    for (uint32_t i = 1; i < k - (n - k) + 1; i++) {
                        d[i] = g[abs(int32_t(n - k - i))] + g[int32_t(n - k + i)];
                    }

                    for (uint32_t i = k - (n - k) + 1; i < n + 1; i++) {
                        d[i] = g[abs(int32_t(i - n + k))];
                    }
                }
                else {
                    d[n - k] = g[0];//d[n - k] = g.front();
                    for (uint32_t i = n - 2 * k; i < n + 1; i++) {
                        if (i != n - k) {
                            d[i] = g[abs(int32_t(i - n + k))];
                        }
                    }
                }
            }

            if (IsNotEqualOne(r[r_len-1])) {//if (IsNotEqualOne(r.back())) {
                // d *= f[n]
//                std::transform(d.begin(), d.end(), d.begin(),
//                               std::bind(std::multiplies<double>(), std::placeholders::_1, r.back()));
                //注：初始化，vector r=f，因此 d*=f[n] 也写作 d*=r[n]，因此下面乘以 r[r_len-1]
                for (int i = 0; i < d_len; ++i) {
                    d[i]=d[i]*r[r_len-1];
                }
            }
            if (IsNotEqualOne(g[g_len-1])) {//if (IsNotEqualOne(g.back())) {
                // d /= g[k]
//                std::transform(d.begin(), d.end(), d.begin(),
//                               std::bind(std::divides<double>(), std::placeholders::_1, g.back()));
                for (int i = 0; i < d_len; ++i) {
                    d[i]=d[i]/g[g_len-1];
                }
            }

            // f-=d
//            std::transform(r.begin(), r.end(), d.begin(), r.begin(), std::minus<double>());
            if (r_len<d_len)
                cout<<"error: r_len<d_len!"<<endl;
            for (int i = 0; i < r_len; ++i) {
                r[i]=r[i]-d[i];
            }
            if (r_len > 1) {//if (r.size() > 1) {
                n = Degree(r,r_len);//n = Degree(r);
                r_len=n+1; //r.resize(n + 1);
            }
        }

        if (n == k) {
            q[0] = r[r_len-1];//q.front() = r.back();
            if (IsNotEqualOne(g[g_len-1])) {//if (IsNotEqualOne(g.back())) {
                // q[0] /= g[k]
                q[0]/=g[g_len-1];//q.front() /= g.back();
            }
            uint32_t d_len=g_len;
            double *d=new double [d_len];//std::vector<double> d = g;
            for (int i = 0; i < d_len; ++i) {
                d[i]=g[i];
            }
            if (IsNotEqualOne(r[r_len-1])) {//if (IsNotEqualOne(r.back())) {
                // d *= f[n]
//                std::transform(d.begin(), d.end(), d.begin(),
//                               std::bind(std::multiplies<double>(), std::placeholders::_1, r.back()));
                //注：初始化，vector r=f，因此 d*=f[n] 也写作 d*=r[n]，因此下面乘以 r[r_len-1]
                for (int i = 0; i < d_len; ++i) {
                    d[i]=d[i]*r[r_len-1];
                }

            }
            if (IsNotEqualOne(g[g_len-1])) {//if (IsNotEqualOne(g.back())) {
                // d /= g[k]
//                std::transform(d.begin(), d.end(), d.begin(),
//                               std::bind(std::divides<double>(), std::placeholders::_1, g.back()));
                for (int i = 0; i < d_len; ++i) {
                    d[i]=d[i]/g[g_len-1];
                }
            }
            // f-=d
//            std::transform(r.begin(), r.end(), d.begin(), r.begin(), std::minus<double>());
            if (r_len<d_len)
                cout<<"error: r_len<d_len!"<<endl;
            for (int i = 0; i < r_len; ++i) {
                r[i]=r[i]-d[i];
            }
            if (r_len > 1) {//if (r.size() > 1) {
                n = Degree(r,r_len);//n = Degree(r);
                r_len=n+1; //r.resize(n + 1);
            }
        }
        // Because we want to have [c0] in the last spot, not [c0/2]
        q[0]*=2;//q.front() *= 2;
    }
    else {
        uint32_t q2_len=1;
        double* q2 = new double[q2_len];//std::vector<double> q2(1, 0.0);
        q_len=q2_len;  //q = q2;
        q = new double[q_len];
        for (int i = 0; i < q2_len; ++i) {
            q2[i]=0.0;
            q[i]=q2[i];
        }

        if (r_len<f_len)
            cout<<"error: r_len<f_len"<<endl;
        for (int i = 0; i < f_len; ++i) {   //r = f;
            r[i]=f[i];
        }

    }

    quotient_len=q_len;
    quotient=q;
    remainder=r;
    remainder_len=r_len;
//    return std::make_shared<longDiv>(q, r);

}

double PREC = std::pow(2, -20);
inline bool IsNotEqualOne(double val) {
    if (1 - PREC >= val) {
        return true;
    }
    if (1 + PREC <= val) {
        return true;
    }
    return false;
}

void EvalLinearWSumMutable(Ciphertext & wsum,
                           Ciphertext * ciphertexts, uint32_t ciphertexts_num, double* constants,
                           Scheme scheme)
{

//    if (cryptoParams->GetScalingTechnique() != FIXEDMANUAL) {
//        // Check to see if input ciphertexts are of same level
//        // and adjust if needed to the max level among them
//        //...省略一些 高阶 scaling technique 带来的复杂的 adjust level and depth 操作
//    }

    // Check to see if input ciphertexts are of same level
    // and adjust if needed to the max level among them
    // NOTE: yx: 注： 本工程中 limbs 是做一次乘法减一的，所以应该所有密文的 limbs 统一到最小的那个值上，与 openfhe 中统一到最大值是相反的。
    uint32_t minLevel = ciphertexts[0].curr_limbs;
    uint32_t minIdx   = 0;
    for (uint32_t i= 1; i < ciphertexts_num; ++i) {
        if (ciphertexts[i].curr_limbs < minLevel){
            minLevel = ciphertexts[i].curr_limbs;
            minIdx   = i;
        }
    }

    for (uint32_t i = 0; i < minIdx; ++i) {
        if(ciphertexts[i].curr_limbs > minLevel)
            scheme.modDownToAndEqual(ciphertexts[i],minLevel);
    }

    for (uint32_t i = minIdx+1; i < ciphertexts_num; ++i) {
        if(ciphertexts[i].curr_limbs > minLevel)
            scheme.modDownToAndEqual(ciphertexts[i],minLevel);
    }

    //weightedSum Core computation
//    wsum = scheme.multByConst(ciphertexts[0], constants[0]);//Ciphertext<DCRTPoly> weightedSum = cc->EvalMult(ciphertexts[0], constants[0]);
    wsum = EvalMult(ciphertexts[0], constants[0],scheme);//Ciphertext<DCRTPoly> weightedSum = cc->EvalMult(ciphertexts[0], constants[0]);
//    wsum = scheme.Lattigo_MultByConst(ciphertexts[0], constants[0]);//Ciphertext<DCRTPoly> weightedSum = cc->EvalMult(ciphertexts[0], constants[0]);
    Ciphertext tmp;//Ciphertext<DCRTPoly> tmp;
    for (uint32_t i = 1; i < ciphertexts_num; i++) {//for (uint32_t i = 1; i < ciphertexts.size(); i++) {
//        tmp = scheme.multByConst(ciphertexts[i], constants[i]);//tmp = cc->EvalMult(ciphertexts[i], constants[i]);
        tmp = EvalMult(ciphertexts[i], constants[i],scheme);//tmp = cc->EvalMult(ciphertexts[i], constants[i]);
//        tmp = scheme.Lattigo_MultByConst(ciphertexts[i], constants[i]);//tmp = cc->EvalMult(ciphertexts[i], constants[i]);
        scheme.addAndEqual(wsum,tmp);//cc->EvalAddInPlace(weightedSum, tmp);
    }
//    scheme.reScaleByAndEqual(wsum, 1);//cc->ModReduceInPlace(weightedSum);
    ModReduceInternalInPlace(wsum,1,scheme);


}

void InnerEvalChebyshevPS(Ciphertext& result, const Ciphertext x,
                                double* coefficients, uint32_t coefficients_len,
                                uint32_t k, uint32_t m,
                                Ciphertext* T,Ciphertext*T2,
                                Scheme scheme)
{
//    auto cc = x->GetCryptoContext();

    // Compute k*2^{m-1}-k because we use it a lot
    uint32_t k2m2k = k * (1 << (m - 1)) - k;

    // Divide coefficients by T^{k*2^{m-1}}
    //std::vector<double> Tkm(int32_t(k2m2k + k) + 1, 0.0);
    uint32_t Tkm_len=int32_t(k2m2k + k) + 1;
    double *Tkm=new double[Tkm_len]();
    for (int i = 0; i < Tkm_len; ++i) {
        Tkm[i]=0.0;
    }
    Tkm[Tkm_len-1]=1;//Tkm.back() = 1;
//    auto divqr = LongDivisionChebyshev(coefficients, Tkm);
    double *divqr_q,*divqr_r;
    uint32_t divqr_q_len, divqr_r_len;
    LongDivisionChebyshev(divqr_q,divqr_q_len,divqr_r,divqr_r_len,coefficients,coefficients_len,Tkm,Tkm_len);
    delete[]Tkm;

    // Subtract x^{k(2^{m-1} - 1)} from r
//    std::vector<double> r2 = divqr->r;
//    if (int32_t(k2m2k - deg_divqr_r_len) <= 0) {//if (int32_t(k2m2k - Degree(divqr->r)) <= 0) {
//        r2[int32_t(k2m2k)] -= 1;
//        r2.resize(Degree(r2) + 1);
//    }
//    else {
//        r2.resize(int32_t(k2m2k + 1), 0.0);
//        r2.back() = -1;
//    }
    uint32_t r2_len;
    double *r2;
    if (int32_t(k2m2k - Degree(divqr_r,divqr_r_len)) <= 0){
        divqr_r[int32_t(k2m2k)]-=1;
        r2_len = Degree(divqr_r,divqr_r_len)+1;
        r2 = new double [r2_len];
        for (int i = 0;(i < divqr_r_len) && (i<r2_len); ++i) {
            r2[i]=divqr_r[i];
        }
        divqr_r[int32_t(k2m2k)]+=1; //处理 r2[int32_t(k2m2k)] -= 1;
    }
    else{
        r2_len=int32_t(k2m2k + 1);
        r2 =new double [r2_len];
        for (int i = 0; (i < divqr_r_len) && (i<r2_len); ++i) {
            r2[i]=divqr_r[i];
        }
        for (int i = divqr_r_len; i < r2_len-1; ++i) {
            r2[i]=0.0;
        }
        r2[r2_len-1]=-1;
    }

    // Divide r2 by q
    //auto divcs = LongDivisionChebyshev(r2, divqr->q);
    double *divcs_q, *divcs_r;
    uint32_t divcs_q_len, divcs_r_len;
    LongDivisionChebyshev(divcs_q, divcs_q_len, divcs_r, divcs_r_len, r2, r2_len, divqr_q, divqr_q_len);
    delete[] r2;

    // Add x^{k(2^{m-1} - 1)} to s
    double *s2;//std::vector<double> s2 = divcs->r;
    uint32_t s2_len=int32_t(k2m2k + 1);//s2.resize(int32_t(k2m2k + 1), 0.0);
    s2=new double[s2_len];
    for (uint32_t i = 0; (i <divcs_r_len)&&(i<s2_len); ++i) {
        s2[i]=divcs_r[i];
    }
    for (uint32_t i = divcs_r_len; i <s2_len-1; ++i) {
        s2[i]=0.0;
    }
    s2[s2_len-1]=1;//s2.back() = 1;

    // Evaluate c at u
    Ciphertext cu;//Ciphertext<DCRTPoly> cu;
    uint32_t dc = Degree(divcs_q,divcs_q_len);//uint32_t dc = Degree(divcs->q);
    bool flag_c = false;
    if (dc >= 1) {
        if (dc == 1) {
            if (divcs_q[1] != 1) {//if (divcs->q[1] != 1) {
//                cu = scheme.multByConst(T[0], divcs_q[1]);//cu = cc->EvalMult(T.front(), divcs->q[1]);
                cu = EvalMult(T[0], divcs_q[1],scheme);//cu = cc->EvalMult(T.front(), divcs->q[1]);
//                cu = scheme.Lattigo_MultByConst(T[0], divcs_q[1]);//cu = cc->EvalMult(T.front(), divcs->q[1]);
//                scheme.reScaleByAndEqual(cu, 1);//cc->ModReduceInPlace(cu);
                ModReduceInternalInPlace(cu,1,scheme);
            } else {
                cu = T[0];//cu = T.front();
            }
        }
        else {
            Ciphertext *ctxs = new Ciphertext[dc](); //std::vector<Ciphertext<DCRTPoly>> ctxs(dc);
            double *weights = new double[dc](); //std::vector<double> weights(dc);

            for (uint32_t i = 0; i < dc; i++) {
                ctxs[i] = T[i];
                weights[i] = divcs_q[i + 1];//weights[i] = divcs->q[i + 1];
            }
            //cu = cc->EvalLinearWSumMutable(ctxs, weights);
            EvalLinearWSumMutable(cu, ctxs, dc, weights, scheme);

            delete[] ctxs;
            delete[] weights;
        }

        // adds the free term (at x^0)
        scheme.addConstAndEqual(cu, divcs_q[0] / 2);//cc->EvalAddInPlace(cu, divcs->q.front() / 2);

        // Need to reduce levels up to the level of T2[m-1].
        //usint levelDiff = T2[m - 1]->GetLevel() - cu->GetLevel();
        //cc->LevelReduceInPlace(cu, nullptr, levelDiff);
        // 把两个cipher的limbs对齐到更小的那个值
        scheme.modDownToAndEqual(cu,T2[m-1].curr_limbs);

        flag_c = true;
    }

    // Evaluate q and s2 at u. If their degrees are larger than k, then recursively apply the Paterson-Stockmeyer algorithm.
    Ciphertext qu;//Ciphertext<DCRTPoly> qu;
    if (Degree(divqr_q,divqr_q_len) > k) {//if (Degree(divqr->q) > k) {
        // Question: yx: InnerEvalChebyshevPS 的作用和意义
//        qu = InnerEvalChebyshevPS(x, divqr->q, k, m - 1, T, T2);
        InnerEvalChebyshevPS(qu, x, divqr_q,divqr_q_len, k, m - 1, T, T2,scheme);
    }
    else {
        // dq = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        double *qcopy; //auto qcopy = divqr->q;
        uint32_t qcopy_len=k; //qcopy.resize(k);
        qcopy=new double [qcopy_len]();
        for (int i = 0; (i < divqr_q_len)&&(i < qcopy_len); ++i) {
            qcopy[i]=divqr_q[i];
        }
        for (uint32_t i = divqr_q_len; i < qcopy_len; ++i) {
            qcopy[i]=0.0;
        }

        uint32_t deg_qcopy=Degree(qcopy,qcopy_len);
        delete []qcopy;
        if (deg_qcopy > 0) {//if (Degree(qcopy) > 0) {
            Ciphertext* ctxs=new Ciphertext[deg_qcopy]();//std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(qcopy));
            double* weights = new double[deg_qcopy](); //std::vector<double> weights(Degree(qcopy));

            for (uint32_t i = 0; i < deg_qcopy; i++) {//for (uint32_t i = 0; i < Degree(qcopy); i++) {
                ctxs[i]    = T[i];
                weights[i] = divqr_q[i + 1];//weights[i] = divqr->q[i + 1];
            }

            //qu = cc->EvalLinearWSumMutable(ctxs, weights);
            EvalLinearWSumMutable(qu,ctxs,deg_qcopy,weights,scheme);

            // the highest order coefficient will always be a power of two up to 2^{m-1} because q is "monic" but the Chebyshev rule adds a factor of 2
            // we don't need to increase the depth by multiplying the highest order coefficient, but instead checking and summing, since we work with m <= 4.
            Ciphertext sum=T[k-1]; //Ciphertext<DCRTPoly> sum = T[k - 1];
            for (uint32_t i = 0; i < log2(divqr_q[divqr_q_len-1]); i++) {//for (uint32_t i = 0; i < log2(divqr->q.back()); i++) {
                sum = scheme.add(sum, sum);//sum = cc->EvalAdd(sum, sum);
            }

            CheckAndAdjusLevel(qu,sum,qu,sum,scheme);
            scheme.addAndEqual(qu, sum);//cc->EvalAddInPlace(qu, sum);

            delete[] ctxs;
            delete[] weights;
        }
        else{
            Ciphertext sum=T[k-1];//Ciphertext<DCRTPoly> sum = T[k - 1];
            for (uint32_t i = 0; i < log2(divqr_q[divqr_q_len-1]); i++) {//for (uint32_t i = 0; i < log2(divqr->q.back()); i++) {
                sum = scheme.add(sum, sum);//sum = cc->EvalAdd(sum, sum);
            }
            qu = sum;
        }

        // adds the free term (at x^0)
        scheme.addConstAndEqual(qu,divqr_q[0]/2); //cc->EvalAddInPlace(qu, divqr->q.front() / 2);
        // The number of levels of qu is the same as the number of levels of T[k-1] or T[k-1] + 1.
        // No need to reduce it to T2[m-1] because it only reaches here when m = 2.
    }

    Ciphertext su;//Ciphertext<DCRTPoly> su;
    uint32_t deg_s2= Degree(s2,s2_len);
    if (deg_s2 > k) {//if (Degree(s2) > k) {
//        su = InnerEvalChebyshevPS(x, s2, k, m - 1, T, T2);
        InnerEvalChebyshevPS(su, x, s2,s2_len, k, m - 1, T, T2,scheme);
    }
    else {
        // ds = k from construction
        // perform scalar multiplication for all other terms and sum them up if there are non-zero coefficients
        double *scopy; //auto scopy = s2;
        uint32_t scopy_len=k; //scopy.resize(k);
        scopy=new double [scopy_len]();
        for (int i = 0; (i < s2_len)&&(i < scopy_len); ++i) {
            scopy[i]=s2[i];
        }
        for (uint32_t i = s2_len; i < scopy_len; ++i) {
            scopy[i]=0.0;
        }

        uint32_t deg_scopy=Degree(scopy,scopy_len);
        delete []scopy;
        if (deg_scopy > 0) {
            Ciphertext * ctxs=new Ciphertext[deg_scopy](); //std::vector<Ciphertext<DCRTPoly>> ctxs(Degree(scopy));
            double * weights=new double [deg_scopy](); //std::vector<double> weights(Degree(scopy));

            for (uint32_t i = 0; i < deg_scopy; i++) {//for (uint32_t i = 0; i < Degree(scopy); i++) {
                ctxs[i]    = T[i];
                weights[i] = s2[i + 1];
            }

            // su = cc->EvalLinearWSumMutable(ctxs, weights);
            EvalLinearWSumMutable(su,ctxs,deg_scopy,weights,scheme);

            // the highest order coefficient will always be 1 because s2 is monic.
            //cc->EvalAddInPlace(su, T[k - 1]);
            Ciphertext tmp_T;
            if (T[k-1].curr_limbs > su.curr_limbs){
                CheckAndAdjusLevel(su,tmp_T,su,T[k-1],scheme);
                scheme.addAndEqual(su, tmp_T);
            }
            else{
                CheckAndAdjusLevel(su,T[k-1],su,T[k-1],scheme);
                scheme.addAndEqual(su, T[k - 1]);
            }

            delete[] ctxs;
            delete[] weights;
        }
        else {
            su = T[k - 1];
        }

        // adds the free term (at x^0)
        scheme.addConstAndEqual(su,s2[0]/2);//cc->EvalAddInPlace(su, s2.front() / 2);
        // The number of levels of su is the same as the number of levels of T[k-1] or T[k-1] + 1. Need to reduce it to T2[m-1] + 1.
        // su = cc->LevelReduce(su, nullptr, su->GetElements()[0].GetNumOfElements() - Lm + 1) ;
        //注：上面两行为原始代码中的注释. FIXEDMANUAL 模式需要执行：cc->LevelReduceInPlace(su, nullptr);
        scheme.modDownByAndEqual(su,1);//cc->LevelReduceInPlace(su, nullptr);
    }

    delete[]s2;

//    Ciphertext result;//Ciphertext<DCRTPoly> result;

    if (flag_c) {
        CheckAndAdjusLevel(T2[m - 1],cu,T2[m - 1],cu,scheme);
        result=scheme.add(T2[m - 1],cu);//result = cc->EvalAdd(T2[m - 1], cu);
    }
    else {
        result =scheme.addConst(T2[m - 1], divcs_q[0] / 2);//result = cc->EvalAdd(T2[m - 1], divcs->q.front() / 2);
    }

    CheckAndAdjusLevel(result,qu,result,qu,scheme);
    myMultAndEqual(result,qu, scheme);//result = cc->EvalMult(result, qu);
//    scheme.reScaleByAndEqual(result,1);//cc->ModReduceInPlace(result);
    ModReduceInternalInPlace(result,1,scheme);

    CheckAndAdjusLevel(result,su,result,su,scheme);
    scheme.addAndEqual(result,su);//cc->EvalAddInPlace(result, su);

//    return result;
}