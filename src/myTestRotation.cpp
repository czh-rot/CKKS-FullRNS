//
// Created by EYx on 2023/4/25.
//

#include "TestScheme.h"
#include "Context.h"
#include "SecretKey.h"
#include "Scheme.h"
#include "EvaluatorUtils.h"
#include "StringUtils.h"
#include "TimeUtils.h"

#include "myTestRotation.h"
#include "../data/testConstValue.h"
#include "myRotation.h"
#include "myKeySwitch.h"

//rotation
void testAutomorphismTransform()
{
    uint64_t tmp_len=4096;
    uint64_t* myinput = new uint64_t[tmp_len]();
    uint64_t* myoutput = new uint64_t[tmp_len]();
    for (int i = 0; i < tmp_len; ++i) {
        myinput[i]=input_automorph[i];
    }

    int N=512;
    AutomorphismTransform(myoutput, myinput, 8, N, 25, automap_VEC);

    int correct=1;
    for (int i = 0; i < tmp_len; ++i) {
        if (myoutput[i]!=output_automorph[i]){
            cout<<i<<" error in automorphism"<<endl;
            correct=0;
        }
    }
    if(correct)
        cout<<"testAutomorphismTransform correct"<<endl;
}

void testPrecomputeAutoMap()
{
    uint32_t N=512;
    uint32_t autoIndex=25;

    uint32_t vec_len=N;
    uint32_t* vec = new uint32_t[vec_len]();
    PrecomputeAutoMap(N, autoIndex, vec);

    int correct=1;
    for (int i = 0; i < vec_len; ++i) {
        if (vec[i]!=automap_VEC[i]){
            cout<<i<<" error in PrecomputeAutoMap"<<endl;
            correct=0;
        }
    }
    if(correct)
        cout<<"testPrecomputeAutoMap correct"<<endl;

}



void testHRotate(long logN, long L, long logp, long dnum, long logSlots) {
    cout << "!!! START TEST ROTATE BATCH with hybrid KS !!!" << endl;
    //-----------------------------------------
    long K = L/dnum;
    Context context(logN, logp, L, K, dnum);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = (1 << logSlots);
    long rot_list[]={1,4,3,2,-1, }; //todo: 当前仅支持通过正值表明左轮转，尚未实现 通过负值表达右轮转。
    long rot_list_len=4;
    complex<double>* mvec = new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        mvec[i]=i*0.1;
    }
    complex<double>** golden_value = new complex<double>* [rot_list_len];
    for (int i = 0; i < rot_list_len; ++i) {
        golden_value[i] = new complex<double> [slots];
        for (int k = 0; k < slots; ++k) {
            golden_value[i][k]=mvec[k];
        }
    }

    for (int i = 0; i < rot_list_len; ++i) {
        long rotSlots=rot_list[i];
        if(rotSlots>0) {
            EvaluatorUtils::leftRotateAndEqual(golden_value[i], slots, rotSlots);
        } else {
            long idx=-rotSlots;
            EvaluatorUtils::rightRotateAndEqual(golden_value[i], slots, idx);
        }
    }

    Ciphertext cipher1 = scheme.encrypt(mvec, slots, L);

    //test HRotate (thesis: Does Fully Homomorphic Encryption Need Compute Acceleration? Algo4)
    Ciphertext * rot_ciphers=new Ciphertext [rot_list_len];
    //fixme: HRotate需要设计对于 负值表示的右轮转的 相应evk的存储
    HRotate_KeyGen(secretKey,rot_list,rot_list_len,scheme);
    HRotate(rot_ciphers,cipher1,rot_list,rot_list_len,scheme);
    for (int i = 0; i < rot_list_len; ++i) {
        complex<double>* dvec5 = scheme.decrypt(secretKey, rot_ciphers[i]);
        StringUtils::showcompare(golden_value[i], dvec5, slots, "HROTATE_"+ to_string(rot_list[i]));
    }

    cout << "!!! END TEST HROTATE !!!" << endl<<endl;
}

void testFastRotate_demo(long logN, long L, long logp, long dnum, long rotSlots, long logSlots) {
    cout << "!!! START TEST FastRotate_demo with hybrid KS !!!" << endl;
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
    if (rotSlots >= 0){
        EvaluatorUtils::leftRotateAndEqual(mvec, slots, rotSlots);
    }else{
        long idx=-rotSlots;
        EvaluatorUtils::rightRotateAndEqual(mvec, slots, idx);
    }

    FastRotate_KeyGen(secretKey,rotSlots,scheme);
    FastRotate_demo(cipher,rotSlots,scheme);

    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    StringUtils::showcompare(mvec, dvec, slots, "FastRotate_demo");
}

void testEvalFastRotation(long logN, long L, long logp, long dnum, long rotSlots, long logSlots) {
    cout << "!!! START TEST EvalFastRotation(OpenFHE) !!!" << endl;
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
    if (rotSlots >= 0){
        EvaluatorUtils::leftRotateAndEqual(mvec, slots, rotSlots);
    }else{
        long idx=-rotSlots;
        EvaluatorUtils::rightRotateAndEqual(mvec, slots, idx);
    }

    ////------------------------------------------------
    ////EvalFastRotation with
    /// 1)FastRotate_KeyGen,
    /// 2)EvalFastRotationPrecompute,
    /// 3)EvalFastRotation
    ////------------------------------------------------

    FastRotate_KeyGen(secretKey, rotSlots, scheme);

    long curr_limbs = cipher.curr_limbs;
    long beta = std::ceil((curr_limbs * 1.0 / scheme.context.alpha));

    //preallocate KS memory
    //total beta groups, each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys, the size of each poly is 1<<context.logN
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *digits = new uint64_t[beta * expand_length];//allocate at where it calls
    EvalFastRotationPrecompute(digits, cipher, scheme);

    uint32_t tmp_len = (curr_limbs << logN);;
    uint64_t *tmp_ax = new uint64_t[tmp_len]();
    uint64_t *tmp_bx = new uint64_t[tmp_len]();
    Ciphertext result(tmp_ax, tmp_bx, cipher.N, cipher.slots, cipher.curr_limbs);
    EvalFastRotation(result, cipher, digits, rotSlots, scheme);
    delete[]digits;

    complex<double>* dvec = scheme.decrypt(secretKey, result);
    StringUtils::showcompare(mvec, dvec, slots, "EvalFastRotation");

}

void testEvalFastRotationExt(long logN, long L, long logp, long dnum, long rotSlots, long logSlots)
{
    cout << "!!! START TEST EvalFastRotationExt(OpenFHE) !!!" << endl;
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
    if (rotSlots >= 0){
        EvaluatorUtils::leftRotateAndEqual(mvec, slots, rotSlots);
    }else{
        long idx=-rotSlots;
        EvaluatorUtils::rightRotateAndEqual(mvec, slots, idx);
    }

    ////------------------------------------------------
    ////EvalFastRotation with
    /// 1)FastRotate_KeyGen,
    /// 2)EvalFastRotationPrecompute,
    /// 3)EvalFastRotationExt,
    /// 4)KeySwitchDown
    ////------------------------------------------------

    FastRotate_KeyGen(secretKey, rotSlots, scheme);

    long curr_limbs = cipher.curr_limbs;
    long beta = std::ceil((curr_limbs * 1.0 / scheme.context.alpha));

    //preallocate KS memory
    //total beta groups, each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys, the size of each poly is 1<<context.logN
    long expand_length = ((curr_limbs + K) << logN);
    uint64_t *digits = new uint64_t[beta * expand_length];//allocate at where it calls
    EvalFastRotationPrecompute(digits, cipher, scheme);

    uint32_t tmpExt_limbs = (curr_limbs + K);
    uint32_t tmpExt_len = (tmpExt_limbs << logN);
    uint64_t *tmpExt_ax = new uint64_t[tmpExt_len]();
    uint64_t *tmpExt_bx = new uint64_t[tmpExt_len]();
    Ciphertext tmpExt(tmpExt_ax, tmpExt_bx, cipher.N, cipher.slots, tmpExt_limbs);
    EvalFastRotationExt(tmpExt, cipher, rotSlots, digits, true, scheme);
    delete[]digits;

    uint32_t tmp_len = (curr_limbs << logN);;
    uint64_t *tmp_ax = new uint64_t[tmp_len]();
    uint64_t *tmp_bx = new uint64_t[tmp_len]();
    Ciphertext result(tmp_ax, tmp_bx, cipher.N, cipher.slots, cipher.curr_limbs);
    KeySwitchDown(result.ax, result.bx, tmpExt.ax, tmpExt.bx, curr_limbs, scheme);


    complex<double> *dvec = scheme.decrypt(secretKey, result);
    StringUtils::showcompare(mvec, dvec, slots, "EvalFastRotationExt");

}

void testPartialSum(long logN, long L, long logp, long dnum, long logSlots) {
    cout << "!!! START TEST PartialSum !!!" << endl;
    //-----------------------------------------
    long K = L/dnum;
    Context context(logN, logp, L, K,dnum);
    long N=context.N;
    cout<<"N: "<<N<<endl
        <<"L: "<<L<<endl
        <<"dnum: "<<dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = (1 << logSlots);
    complex<double>* mvec = new complex<double> [slots];
    complex<double>* golden_value = new complex<double> [slots];
    for (int i = 0; i < slots; ++i) {
        mvec[i]=i*0.1;
        golden_value[i]=mvec[i];
    }

    Ciphertext cipher = scheme.encrypt(mvec, slots, L);

    for (int j = 1; j < N / (2 * slots); j <<= 1) {
        Ciphertext temp=cipher;
        int rotSlots=j*slots;
        FastRotate_KeyGen(secretKey,rotSlots,scheme); //todo: 将密钥生成挪到预计算中去
        FastRotate_demo(temp,rotSlots,scheme);

        scheme.addAndEqual(cipher, temp);//cc->EvalAddInPlace(raised, temp);

        for (int i = 0; i < slots; ++i) {
            golden_value[i]*=2;
        }
    }
    complex<double>* dvec = scheme.decrypt(secretKey, cipher);
    StringUtils::showcompare(golden_value, dvec, slots, "PartialSum");
}

void testConjugate_demo(long logN, long L, long logp, long dnum, long logSlots)
{
    cout << "!!! START TEST Conjugate_demo with hybrid KS !!!" << endl;
    //-----------------------------------------
    long K = L/dnum;
    Context context(logN, logp, L, K,dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots = (1 << logSlots);
    complex<double>* mvec = EvaluatorUtils::randomComplexArray(slots);
    complex<double>* mvecconj = new complex<double>[slots];
    for (long i = 0; i < slots; ++i) {
//        cout<<mvec[i]<<",\t";
        mvecconj[i] = conj(mvec[i]);
    }
//    cout<<endl;
    Ciphertext cipher = scheme.encrypt(mvec, slots, L);

    Conjugate_KeyGen(secretKey,scheme);
    Conjugate_demo(cipher,scheme);

    complex<double>* dvecconj = scheme.decrypt(secretKey, cipher);
    StringUtils::showcompare(mvecconj, dvecconj, slots, "testConjugate_demo");
}
