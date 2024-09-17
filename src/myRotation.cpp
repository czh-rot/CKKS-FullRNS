//
// Created by EYx on 2023/4/18.
//

#include "myRotation.h"
#include "myKeySwitch.h"
#include "Scheme.h"
#include "../data/testConstValue.h"

void HRotate_KeyGen(SecretKey secretKey, long* rot_list, long rot_list_len, Scheme& scheme)
{
    for (int i = 0; i < rot_list_len; ++i) {
        if (rot_list[i]>=0){
            scheme.addLeftRotKey(secretKey,rot_list[i]);
        }
        else{
            //todo：yx: 右轮转，to be implemented
            long idx = scheme.context.Nh - rot_list[i];
            scheme.addLeftRotKey(secretKey,idx);
        }

    }
}

//todo: yx: 注：该函数尚未支持 右轮转，也没有想好右轮转是通过负值表示还是增加一个 bool 表明左右，HRotate_KeyGen也需要相应处理。
void HRotate(Ciphertext*& out_ct, Ciphertext in_ct, long* rot_list, long rot_list_len, Scheme scheme)
{

    uint32_t curr_limbs= in_ct.curr_limbs;
    uint32_t N=in_ct.N;
    uint32_t M=N<<1;
    uint32_t alpha=scheme.context.alpha;
    uint32_t logN=scheme.context.logN;
    uint32_t K=scheme.context.K;

    ////------------------------------------------------
    ////KeySwitchCore in_ct.ax, now [bx=rot(in_ct.ax)]!
    ////------------------------------------------------

    //KeySwitchCore 操作对象
    uint64_t *KS_input = in_ct.ax;

    //precompute KS params
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //preallocate KS memory
    //total beta groups, each group expand (alpha) RNS polys to (curr_limbs+K) RNS polys, the size of each poly is N
    long expand_limbs = (curr_limbs + K);
    long expand_length = (expand_limbs << logN);
    uint64_t *aHat = new uint64_t[beta * expand_length];

    //Algo4 Line2-3
    EvalKeySwitchPrecomputeCore(aHat, KS_input, curr_limbs, beta, scheme);

    for (int i = 0; i < rot_list_len; ++i) {
        long curr_rotSlots = 0;
        if (rot_list[i] == 0) continue;
        else if (rot_list[i] > 0) {
            curr_rotSlots = rot_list[i];
        } else {
            cout<<"Right Rotation to be implemented!"<<endl;
//            return;
            curr_rotSlots = scheme.context.Nh-rot_list[i];//fixme: yx: 这么写不对。需要结合HRotate_KeyGen修改。
        }

        //Algo4 Line5
        uint64_t* aHat_rot = new uint64_t[beta * expand_length];
        long autoIndex = FindAutomorphismIndex2nComplex(curr_rotSlots, M); //该函数实现、接口皆参考OpenFHE
        uint32_t vec_len=N;
        uint32_t * vec  = new uint32_t [vec_len]();
        PrecomputeAutoMap(N, autoIndex, vec); //该函数实现、接口皆参考OpenFHE
        //automorph beta digits
        for (int i = 0; i < beta; ++i) {
            uint64_t *aHat_rotj = aHat_rot + expand_length * i;
            uint64_t *aHatj = aHat + expand_length * i;
            AutomorphismTransform(aHat_rotj, aHatj, expand_limbs, N, autoIndex, vec); //该函数实现、接口皆参考OpenFHE
        }

        //Algo4 Line6-7
        //here res_ax is directly replaced by in_ct.ax
        long res_len=curr_limbs << logN;
        uint64_t* res_ax = new uint64_t[res_len]();
        uint64_t* res_bx = new uint64_t[res_len]();
        EvalFastKeySwitchCore(res_ax, res_bx, aHat_rot, curr_rotSlots, scheme.leftRotKeyMap,
                              expand_length, curr_limbs, beta, scheme);
        delete[] aHat_rot;

        //Algo4 Line8
        uint64_t* bxrot = new uint64_t[curr_limbs << logN]();
        AutomorphismTransform(bxrot, in_ct.bx, curr_limbs, N, autoIndex, vec);
        delete[] vec;

        //Algo4 Line9
        scheme.context.addAndEqual(res_bx, bxrot, curr_limbs);
        delete[] bxrot;
        out_ct[i]=Ciphertext(res_ax,res_bx,in_ct.N,in_ct.slots,in_ct.curr_limbs);
    }
    delete[] aHat;
}

void FastRotate_KeyGen_List(SecretKey secretKey, int* index, long index_len, Scheme& scheme)
{
    long L=scheme.context.L;
    long K=scheme.context.K;
    long logN=scheme.context.logN;
    long N=scheme.context.N;
    long M=N<<1;

    uint64_t* sxrot = new uint64_t[(L + K) << logN]();
    uint32_t vec_key_len=N;
    uint32_t * vec_key  = new uint32_t [vec_key_len]();
    for (long i = 0; i < index_len; ++i) {
        long autoIndex= FindAutomorphismIndex2nComplex(index[i], M);//该函数实现、接口皆参考OpenFHE
        uint64_t autoIndexInv = invMod(autoIndex,2 * N);//注：这里 invMod函数要取ind的最小正模逆！所以这个invMod函数是从原始库里改过的！
        autoIndexInv = (autoIndexInv*autoIndex)%(2*N);//这一步现在已经添加到上述 invMod函数中

        PrecomputeAutoMap(N, autoIndexInv, vec_key);

        AutomorphismTransform(sxrot, secretKey.sx, L+K, N, autoIndexInv, vec_key);
        // FastRotate 的 evk通过 autoIndex=automorph(index) 进行索引
        scheme.KS_KeyGen(sxrot, secretKey.sx, autoIndex, scheme.BootstrapFastRotationKeyMap);//对scheme类中的leftRotKeyMap进行写操作
    }
    delete []vec_key;
    delete []sxrot;
}

void FastRotate_KeyGen(SecretKey secretKey, int index, Scheme& scheme)
{
    long L=scheme.context.L;
    long K=scheme.context.K;
    long logN=scheme.context.logN;
    long N=scheme.context.N;
    long M=N<<1;

    long autoIndex= FindAutomorphismIndex2nComplex(index, M);//该函数实现、接口皆参考OpenFHE
    uint64_t autoIndexInv = invMod(autoIndex,2 * N);//注：这里 invMod函数要取ind的最小正模逆！所以这个invMod函数是从原始库里改过的！
    autoIndexInv = (autoIndexInv*autoIndex)%(2*N);//这一步现在已经添加到上述 invMod函数中

    uint32_t vec_key_len=N;
    uint32_t * vec_key  = new uint32_t [vec_key_len]();
    PrecomputeAutoMap(N, autoIndexInv, vec_key);

    uint64_t* sxrot = new uint64_t[(L + K) << logN]();
    AutomorphismTransform(sxrot, secretKey.sx, L+K, N, autoIndexInv, vec_key);
    // FastRotate 的evk通过 autoIndex=automorph(index) 进行索引
    scheme.KS_KeyGen(sxrot, secretKey.sx, autoIndex, scheme.leftRotKeyMap);//对scheme类中的leftRotKeyMap进行写操作
    delete []sxrot;
    delete []vec_key;
}

void FastRotate_demo(Ciphertext& cipher, int index, Scheme scheme)
{
    uint32_t curr_limbs= cipher.curr_limbs;
    uint32_t N=cipher.N;
    uint32_t M=N<<1;
    uint32_t logN=scheme.context.logN;

    long autoIndex = FindAutomorphismIndex2nComplex(index, M); //该函数实现、接口皆参考OpenFHE

    ////------------------------------------------------
    ////KeySwitchCore cipher.ax, now [bx=rot(cipher.ax)]!
    ////------------------------------------------------

    //KeySwitchCore 操作对象
    uint64_t* KS_input=cipher.ax;

    long res_len=curr_limbs << logN;
    uint64_t* res_ax = new uint64_t[res_len]();
    uint64_t* res_bx = new uint64_t[res_len]();
    // FastRotate 的evk通过 autoIndex=automorph(index) 进行索引
    KeySwitchCore(res_ax, res_bx, KS_input, autoIndex, scheme.leftRotKeyMap, curr_limbs, scheme);

    uint64_t* bxrot = new uint64_t[curr_limbs << logN]();
    scheme.context.add( bxrot, cipher.bx, res_bx,  curr_limbs);

    uint32_t vec_len=N;
    uint32_t * vec  = new uint32_t [vec_len]();
    PrecomputeAutoMap(N, autoIndex, vec); //该函数实现、接口皆参考OpenFHE

    AutomorphismTransform(cipher.ax, res_ax, curr_limbs, N, autoIndex, vec);
    AutomorphismTransform(cipher.bx,bxrot, curr_limbs, N, autoIndex, vec);
    delete [] vec;
}

void EvalFastRotationPrecompute(uint64_t *digits, const Ciphertext &cipher, Scheme scheme) {
    long alpha = scheme.context.alpha;
    long curr_limbs = cipher.curr_limbs;
    long beta = std::ceil((curr_limbs * 1.0 / alpha));

    //openfhe 的 cv[1] 多项式对应这里的 cipher.ax
    //base-leveledshe.cpp
    // algo->EvalKeySwitchPrecomputeCore(cv[1], ciphertext->GetCryptoParameters());
    EvalKeySwitchPrecomputeCore(digits, cipher.ax, curr_limbs, beta, scheme);
}

void EvalFastRotation(Ciphertext& result, Ciphertext& cipher,
                      uint64_t* digits, int index,
                      Scheme scheme)
{
    if (index == 0) {
        result=cipher; //Ciphertext<Element> result = ciphertext->Clone(); return result;
    }

    //KS params
    long curr_limbs=cipher.curr_limbs;
    uint32_t alpha = scheme.context.alpha;
    uint32_t logN = scheme.context.logN;
    uint32_t N = 1 << logN;
    uint32_t M = N << 1;
    uint32_t K = scheme.context.K;

    long autoIndex = FindAutomorphismIndex2nComplex(index, M);

    long beta = std::ceil(
            (curr_limbs * 1.0 / alpha));
    long expand_limbs = (curr_limbs + K);
    long expand_length = (expand_limbs << logN);
    uint32_t res_len = (curr_limbs << logN);
    uint64_t *res_ax = new uint64_t[res_len]();
    uint64_t *res_bx = new uint64_t[res_len]();
    // FastRotate 的evk通过 autoIndex=automorph(index) 进行索引
    EvalFastKeySwitchCore(res_ax, res_bx, digits, autoIndex, scheme.leftRotKeyMap,
                          expand_length, curr_limbs, beta, scheme);

    scheme.context.addAndEqual(res_bx, cipher.bx, curr_limbs);

    uint32_t vec_len=N;
    uint32_t * vec  = new uint32_t [vec_len]();
    PrecomputeAutoMap(N, autoIndex, vec);//  PrecomputeAutoMap(N, autoIndex, &vec);

    AutomorphismTransform(result.ax,res_ax,curr_limbs,N,autoIndex,vec);
    AutomorphismTransform(result.bx,res_bx,curr_limbs,N,autoIndex,vec);

    delete [] vec;
    delete []res_ax;
    delete []res_bx;
}

void EvalFastRotationExt(Ciphertext& resultExt, Ciphertext cipher, int index,
                         uint64_t* digits, bool addFirst, Scheme scheme )
{
    if (index == 0) {
        resultExt=cipher; //Ciphertext<Element> resultExt = ciphertext->Clone(); return resultExt;
    }

    uint32_t curr_limbs= cipher.curr_limbs;
    uint32_t N = cipher.N;
    uint32_t M = N << 1;
    uint32_t alpha = scheme.context.alpha;
    uint32_t logN = scheme.context.logN;
    uint32_t K = scheme.context.K;
    long beta = std::ceil((curr_limbs * 1.0 / alpha));//current level=curr_limbs, therefore there's curr_limbs modulus

    // Find the automorphism index that corresponds to rotation index.
    long autoIndex = FindAutomorphismIndex2nComplex(index, M); //该函数实现、接口皆参考OpenFHE

    long expand_limbs = (curr_limbs + K);
    long expand_length = (expand_limbs << logN);
    uint64_t *sumaxmult = new uint64_t[expand_length]();
    uint64_t *sumbxmult = new uint64_t[expand_length]();

    // InnerProduct
    // FastRotate 的evk通过 autoIndex=automorph(index) 进行索引
    EvalFastKeySwitchCoreExt(sumaxmult, sumbxmult, digits, autoIndex, scheme.leftRotKeyMap,
                             expand_length, beta, curr_limbs, scheme);

    if (addFirst) {
        // 计算 cMult = cipher.bx * PModQ 实现对cipher.bx的升模，升模的结果为cMult，与 sumbxmult 对齐 模数 以进行下一步计算
        uint64_t *cMult = new uint64_t[curr_limbs << logN];
        for (int i = 0; i < curr_limbs; i++) {
            uint64_t *cMultj = cMult + (i << logN);
            uint64_t *cipher_bxj = cipher.bx + (i << logN);
            uint64_t PModqi = scheme.context.PModq[i];
            for (int j = 0; j < N; ++j) {
                //compute cMultj[i]=cipher_bxj[i]*PModqi;
                mulModBarrett(cMultj[j], cipher_bxj[j], PModqi,
                              scheme.context.qVec[i], scheme.context.qrVec[i], scheme.context.qTwok[i]);
            }
        }
        // sumbxmult+=cipher.bx
        scheme.context.addAndEqual(sumbxmult, cMult, curr_limbs);
    }

    uint32_t vec_len = N;
    uint32_t *vec = new uint32_t[vec_len]();
    PrecomputeAutoMap(N, autoIndex, vec);//  PrecomputeAutoMap(N, autoIndex, &vec);

    AutomorphismTransform(resultExt.ax, sumaxmult, expand_limbs, N, autoIndex, vec);
    AutomorphismTransform(resultExt.bx, sumbxmult, expand_limbs, N, autoIndex, vec);

    delete[] sumaxmult;
    delete[] sumbxmult;
}


////------------------------------------------------
////rotation utils
////------------------------------------------------
//todo: yx: set as another precompute map?
long FindAutomorphismIndex2nComplex(int32_t i, uint32_t m) {
    if (i == 0) {
        return 1;
    }

    // conjugation automorphism
    if (i == int32_t(m - 1)) {
        return long(i);
    }
    else {
        // generator
        int32_t g0;

        if (i < 0) {
//            g0 = NativeInteger(5).ModInverse(m).ConvertToInt();
            g0 = invMod(5, m);
            g0 = (g0*5)%m;
        }
        else {
            g0 = 5;
        }
        uint32_t i_unsigned = (uint32_t)std::abs(i);

        int32_t g = g0;
        for (size_t j = 1; j < i_unsigned; j++) {
            g = (g * g0) % m;
        }
        return long(g);
    }
}

void PrecomputeAutoMap(uint32_t n, uint32_t k, uint32_t* precomp) {
    uint32_t m    = n << 1;  // cyclOrder
    uint32_t logm = std::round(log2(m));
    uint32_t logn = std::round(log2(n));
    for (uint32_t j = 0; j < n; j++) {
        uint32_t jTmp    = ((j << 1) + 1);
        uint32_t idx        = ((jTmp * k) - (((jTmp * k) >> logm) << logm)) >> 1;
        uint32_t jrev       = ReverseBits(j, logn);
        uint32_t idxrev     = ReverseBits(idx, logn);
        precomp[jrev] = idxrev;
    }
}

inline uint64_t ReverseBits(uint64_t num, uint64_t msb) {
    uint64_t msbb = (msb >> 3) + (msb & 0x7 ? 1 : 0);

    switch (msbb) {
        case 1:
            return (reverse_byte((num)&0xff) >> shift_trick[msb & 0x7]);

        case 2:
            return (reverse_byte((num)&0xff) << 8 | reverse_byte((num >> 8) & 0xff)) >> shift_trick[msb & 0x7];

        case 3:
            return (reverse_byte((num)&0xff) << 16 | reverse_byte((num >> 8) & 0xff) << 8 |
                    reverse_byte((num >> 16) & 0xff)) >>
                                                      shift_trick[msb & 0x7];

        case 4:
            return (reverse_byte((num)&0xff) << 24 | reverse_byte((num >> 8) & 0xff) << 16 |
                    reverse_byte((num >> 16) & 0xff) << 8 | reverse_byte((num >> 24) & 0xff)) >>
                                                                                              shift_trick[msb & 0x7];
        default:
            return -1;
            // OPENFHE_THROW(math_error, "msbb value not handled:" +
            // std::to_string(msbb));
    }
}

inline static unsigned char reverse_byte(unsigned char x) {
    static const unsigned char table[] = {
            0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0, 0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0, 0x70, 0xf0, 0x08, 0x88,
            0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8, 0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8, 0x04, 0x84, 0x44, 0xc4,
            0x24, 0xa4, 0x64, 0xe4, 0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4, 0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac,
            0x6c, 0xec, 0x1c, 0x9c, 0x5c, 0xdc, 0x3c, 0xbc, 0x7c, 0xfc, 0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2,
            0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2, 0x0a, 0x8a, 0x4a, 0xca, 0x2a, 0xaa, 0x6a, 0xea, 0x1a, 0x9a,
            0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa, 0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6, 0x16, 0x96, 0x56, 0xd6,
            0x36, 0xb6, 0x76, 0xf6, 0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee, 0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe,
            0x7e, 0xfe, 0x01, 0x81, 0x41, 0xc1, 0x21, 0xa1, 0x61, 0xe1, 0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
            0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9, 0x19, 0x99, 0x59, 0xd9, 0x39, 0xb9, 0x79, 0xf9, 0x05, 0x85,
            0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5, 0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5, 0x0d, 0x8d, 0x4d, 0xcd,
            0x2d, 0xad, 0x6d, 0xed, 0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd, 0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3,
            0x63, 0xe3, 0x13, 0x93, 0x53, 0xd3, 0x33, 0xb3, 0x73, 0xf3, 0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb,
            0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb, 0x07, 0x87, 0x47, 0xc7, 0x27, 0xa7, 0x67, 0xe7, 0x17, 0x97,
            0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7, 0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef, 0x1f, 0x9f, 0x5f, 0xdf,
            0x3f, 0xbf, 0x7f, 0xff,
    };
    return table[x];
}

void AutomorphismTransform(uint64_t* ra, uint64_t* a,
                           uint32_t l, uint32_t N,
                           uint32_t i, uint32_t* precomp_vec)
{
//    DCRTPolyType result(*this);
    for (uint32_t k = 0; k < l; k++) { //for (usint k = 0; k < m_vectors.size(); k++) {
        uint64_t* raj=ra+(k * N);//todo: yx: change to logn with right shift
        uint64_t* aj=a+(k * N);//todo: yx: change to logn with right shift
        // result.m_vectors[k] = m_vectors[k].AutomorphismTransform(i, vec);展平为以下代码
        if (i % 2 == 0) {
//            OPENFHE_THROW(math_error, "automorphism index should be odd\n");
            cout<<"automorphism index should be odd\n"<<endl;
        }
        for (uint32_t j = 0; j < N; ++j) {
                raj[j] = aj[precomp_vec[j]]; //(*result.m_values)[j] = (*m_values)[precomp[j]];
        }
    }
//    return result;
}

void Conjugate_KeyGen(SecretKey secretKey, Scheme& scheme)
{
    long L=scheme.context.L;
    long K=scheme.context.K;
    long logN=scheme.context.logN;
    long N=scheme.context.N;
    long M=N<<1;

    long autoIndex= 2 * N - 1;//该函数实现、接口皆参考OpenFHE
    uint64_t autoIndexInv = invMod(autoIndex,2 * N);//注：这里 invMod函数要取ind的最小正模逆！所以这个invMod函数是从原始库里改过的！
    autoIndexInv = (autoIndexInv*autoIndex)%(2*N);//这一步现在已经添加到上述 invMod函数中

    uint32_t vec_key_len=N;
    uint32_t * vec_key  = new uint32_t [vec_key_len]();
    PrecomputeAutoMap(N, autoIndexInv, vec_key);

    uint64_t* sxrot = new uint64_t[(L + K) << logN]();
    AutomorphismTransform(sxrot, secretKey.sx, L+K, N, autoIndexInv, vec_key);
    // FastRotate 的evk通过 autoIndex=automorph(index) 进行索引
    scheme.KS_KeyGen(sxrot, secretKey.sx, autoIndex, scheme.leftRotKeyMap);//对scheme类中的leftRotKeyMap进行写操作
    delete []sxrot;
    delete []vec_key;
}

void Conjugate_demo(Ciphertext& cipher, Scheme scheme)
{
    uint32_t curr_limbs= cipher.curr_limbs;
    uint32_t N=cipher.N;
    uint32_t M=N<<1;
    uint32_t logN=scheme.context.logN;

    long autoIndex = 2 * N - 1; //该函数实现、接口皆参考OpenFHE

    ////------------------------------------------------
    ////KeySwitchCore cipher.ax, note that now [bx=rot(cipher.ax)]!
    ////------------------------------------------------

    //KeySwitchCore 操作对象
    uint64_t* KS_input=cipher.ax;

    long res_len=curr_limbs << logN;
    uint64_t* res_ax = new uint64_t[res_len]();
    uint64_t* res_bx = new uint64_t[res_len]();
    // Conjugate 的evk通过 autoIndex 进行索引
    KeySwitchCore(res_ax, res_bx, KS_input, autoIndex, scheme.leftRotKeyMap, curr_limbs, scheme);

    uint64_t* bxrot = new uint64_t[curr_limbs << logN]();
    scheme.context.add( bxrot, cipher.bx, res_bx,  curr_limbs);

    uint32_t vec_len=N;
    uint32_t * vec  = new uint32_t [vec_len]();
    PrecomputeAutoMap(N, autoIndex, vec); //该函数实现、接口皆参考OpenFHE

    AutomorphismTransform(cipher.ax, res_ax, curr_limbs, N, autoIndex, vec);
    AutomorphismTransform(cipher.bx,bxrot, curr_limbs, N, autoIndex, vec);
    delete [] vec;
}
