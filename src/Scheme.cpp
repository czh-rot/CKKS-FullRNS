/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "../data/BSConstValue.h"
#include "../data/testConstValue.h"
#include "Scheme.h"
#include "Common.h"
#include "myUtils.h"
//Scheme::Scheme(Context& context) : context(context) {}

Scheme::Scheme(SecretKey &secretKey, Context &context) : secretKey(secretKey), context(context) {
    addEncKey(secretKey);
    addMultKey(secretKey);
}

Scheme::Scheme(SecretKey& secretKey, Context& context,
               ScalingTechnique rescaleTech, SecretKeyDist secretKeyDist,uint64_t L0) :
secretKey(secretKey),context(context), rescaleTech(rescaleTech), secretKeyDist(secretKeyDist),L0(L0) {
	addEncKey(secretKey);
	addMultKey(secretKey);
}


void Scheme::addEncKey(SecretKey& secretKey) {
	uint64_t* ex = new uint64_t[context.L << context.logN]();
	uint64_t* ax = new uint64_t[context.L << context.logN]();
	uint64_t* bx = new uint64_t[context.L << context.logN]();

	context.sampleUniform(ax, context.L);

	context.sampleGauss(ex, context.L);
	context.NTTAndEqual(ex, context.L);

	context.mul(bx, ax, secretKey.sx, context.L);
	context.sub2AndEqual(ex, bx, context.L);

	delete[] ex;

	keyMap.insert(pair<long, Key>(ENCRYPTION, Key(ax, bx)));
}

/*
 * @param beforeKS_sk: KS操作前解密需要的新的解密密钥
 * @param afterKS_sk: KS操作后解密密文使用的密钥
 */
void Scheme::KS_KeyGen(uint64_t* afterKS_sk, uint64_t* beforeKS_sk, long type, map<long, Key>& keyMap){
    for (int j = 0; j < context.dnum; ++j) {
        uint64_t* ax = new uint64_t[(context.L + context.K) << context.logN]();
        uint64_t* bx = new uint64_t[(context.L + context.K) << context.logN]();
        uint64_t* ex = new uint64_t[(context.L + context.K) << context.logN]();
        uint64_t* fac_x_newKey = new uint64_t[(context.L + context.K) << context.logN]();

        // factor * beforeKS_sk (mod Q)
        context.myEvalAndEqual(fac_x_newKey, beforeKS_sk, context.L, j);
        context.sampleGauss(ex, context.L, context.K);
        //fixme: ex 若不置零 hybrid KS 计算过程无法正确解密。ex置零影响安全性，不影响正确性。猜测为噪声溢出。
        for (int i = 0; i < (context.L + context.K) << context.logN; ++i) {
            ex[i] = 0;
        }
        context.NTTAndEqual(ex, context.L, context.K);
        context.addAndEqual(ex, fac_x_newKey, context.L);
        if(SET_SWK==1){
            //SET swk
            for (int i = 0; i < context.L + context.K; ++i) {
                uint64_t *axj = ax + (i << context.logN);
                for (int k = 0; k < context.N; ++k) {
                    axj[k] = BSswk_ax[i][k];
                }
            }
        }
        else{
            context.sampleUniform(ax, context.L, context.K);
        }
        context.mul(bx, ax, afterKS_sk, context.L, context.K);
        context.sub2AndEqual(ex, bx, context.L, context.K);
        keyMap.insert(pair<long, Key>(type * context.dnum + j, Key(ax, bx)));


        delete[] ex;
        delete[] fac_x_newKey;
    }
}

void Scheme::addMultKey(SecretKey& secretKey) {
    uint64_t* sxsx = new uint64_t[(context.L + context.K) << context.logN]();
    context.mul(sxsx, secretKey.sx, secretKey.sx, context.L);

    KS_KeyGen(secretKey.sx,sxsx,MULTIPLICATION,keyMap);

//    uint64_t* ax = new uint64_t[(context.L + context.K) << context.logN]();
//    uint64_t* bx = new uint64_t[(context.L + context.K) << context.logN]();
//    uint64_t* ex = new uint64_t[(context.L + context.K) << context.logN]();
//    context.evalAndEqual(sxsx, context.L);
//
//	context.sampleGauss(ex, context.L, context.K);
//	context.NTTAndEqual(ex, context.L, context.K);
//
//	context.addAndEqual(ex, sxsx, context.L);
//
//	context.sampleUniform(ax, context.L, context.K);
//	context.mul(bx, ax, secretKey.sx, context.L, context.K);
//	context.sub2AndEqual(ex, bx, context.L, context.K);
//	delete[] ex;
//    keyMap.insert(pair<long, Key>(MULTIPLICATION, Key(ax, bx)));

    delete[] sxsx;

}

void Scheme::addConjKey(SecretKey& secretKey) {
	uint64_t* sxconj = new uint64_t[(context.L + context.K) << context.logN]();

    context.conjugate(sxconj, secretKey.sx, context.L);

    KS_KeyGen(secretKey.sx,sxconj,CONJUGATION,keyMap);

//    uint64_t* ex = new uint64_t[(context.L + context.K) << context.logN]();
//    uint64_t* ax = new uint64_t[(context.L + context.K) << context.logN]();
//    uint64_t* bx = new uint64_t[(context.L + context.K) << context.logN]();
//    context.evalAndEqual(sxconj, context.L);
//
//	context.sampleGauss(ex, context.L, context.K);
//	context.NTTAndEqual(ex, context.L, context.K);
//
//	context.addAndEqual(ex, sxconj, context.L);
//
//	context.sampleUniform(ax, context.L, context.K);
//	context.mul(bx, ax, secretKey.sx, context.L, context.K);
//	context.sub2AndEqual(ex, bx, context.L, context.K);
//    delete[] ex;
//    keyMap.insert(pair<long, Key>(CONJUGATION, Key(ax, bx)));

    delete[] sxconj;

}

void Scheme::addLeftRotKey(SecretKey& secretKey, long rot) {
	uint64_t* sxrot = new uint64_t[(context.L + context.K) << context.logN]();

    context.leftRot(sxrot, secretKey.sx, context.L, rot);
    KS_KeyGen(secretKey.sx,sxrot,rot,leftRotKeyMap);

//    uint64_t* ex = new uint64_t[(context.L + context.K) << context.logN]();
//    uint64_t* ax = new uint64_t[(context.L + context.K) << context.logN]();
//    uint64_t* bx = new uint64_t[(context.L + context.K) << context.logN]();
//    context.evalAndEqual(sxrot, context.L);
//
//	context.sampleGauss(ex, context.L, context.K);
//	context.NTTAndEqual(ex, context.L, context.K);
//
//	context.addAndEqual(ex, sxrot, context.L);
//
//	context.sampleUniform(ax, context.L, context.K);
//	context.mul(bx, ax, secretKey.sx, context.L, context.K);
//    context.sub2AndEqual(ex, bx, context.L, context.K);
//    delete[] ex;
//    leftRotKeyMap.insert(pair<long, Key>(rot, Key(ax, bx)));

    delete[] sxrot;

}

void Scheme::addLeftRotKeys(SecretKey& secretKey) {
	for (long i = 0; i < context.logNh; ++i) {
		long idx = 1 << i;
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx);
		}
	}
}

void Scheme::addRightRotKeys(SecretKey& secretKey) {
	for (long i = 0; i < context.logNh; ++i) {
		long idx = context.Nh - (1 << i);
		if(leftRotKeyMap.find(idx) == leftRotKeyMap.end()) {
			addLeftRotKey(secretKey, idx);
		}
	}
}

Plaintext Scheme::encode(double* v, long slots, long l) {
	uint64_t* m = new uint64_t[l << context.logN]();
	context.encode(m, v, slots, l);
	return Plaintext(m, context.N, slots, l);
}

Plaintext Scheme::encode(complex<double>* v, long slots, long l) {
	uint64_t* m = new uint64_t[l << context.logN]();
	context.encode(m, v, slots, l);
	return Plaintext(m, context.N, slots, l);
}

Plaintext Scheme::encodeSingle(complex<double> val, long l) {
	uint64_t* m = new uint64_t[l << context.logN]();
	context.encodeSingle(m, val, l);
	return Plaintext(m, context.N, 1, l);
}

complex<double>* Scheme::decode(Plaintext& msg) {
	complex<double>* res = new complex<double>[msg.slots]();
	context.decode(msg.mx, res, msg.slots, msg.l);
	return res;
}

complex<double> Scheme::decodeSingle(Plaintext& msg) {
	complex<double> res;
	context.decodeSingle(msg.mx, res, msg.l);
	return res;
}

Ciphertext Scheme::encryptMsg(SecretKey& secretkey, Plaintext& message) {
	Ciphertext res;
	return res;
}


Ciphertext Scheme::encryptMsg(Plaintext& message) {
	Key key = keyMap.at(ENCRYPTION);

	uint64_t* ax = new uint64_t[message.l << context.logN]();
	uint64_t* bx = new uint64_t[message.l << context.logN]();
	uint64_t* vx = new uint64_t[message.l << context.logN]();
	uint64_t* ex = new uint64_t[message.l << context.logN]();

	context.sampleZO(vx, context.Nh, message.l);
	context.NTTAndEqual(vx, message.l);

	context.mul(ax, vx, key.ax, message.l);

	context.sampleGauss(ex, message.l);
	context.NTTAndEqual(ex, message.l);

    //Debug Only: set ex in encryption to zero
    for (int i = 0; i < message.l << context.logN; ++i) {
        ex[i]=0;
    }
	context.addAndEqual(ax, ex, message.l);

	context.mul(bx, vx, key.bx, message.l);

	context.sampleGauss(ex, message.l);
	context.NTTAndEqual(ex, message.l);

	context.addAndEqual(bx, ex, message.l);
	context.addAndEqual(bx, message.mx, message.l);

	return Ciphertext(ax, bx, context.N, message.slots, message.l);
}

Plaintext Scheme::decryptMsg(SecretKey& secretKey, Ciphertext& cipher) {
	uint64_t* mx = new uint64_t[context.N]();
	context.mul(mx, cipher.ax, secretKey.sx, 1);
	context.addAndEqual(mx, cipher.bx, 1);

	return Plaintext(mx, context.N, cipher.slots, 1);
}

Ciphertext Scheme::encrypt(double* vals, long slots, long l) {
	Plaintext msg = encode(vals, slots, l);
	return encryptMsg(msg);
}

Ciphertext Scheme::encrypt(complex<double>* vals, long slots, long l) {
	Plaintext msg = encode(vals, slots, l);
	return encryptMsg(msg);
}

Ciphertext Scheme::encryptSingle(complex<double> val, long l) {
	Plaintext msg = encodeSingle(val, l);
	return encryptMsg(msg);
}

complex<double>* Scheme::decrypt(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg = decryptMsg(secretKey, cipher);
	return decode(msg);
}

complex<double> Scheme::decryptSingle(SecretKey& secretKey, Ciphertext& cipher) {
	Plaintext msg = decryptMsg(secretKey, cipher);
	return decodeSingle(msg);
}

Ciphertext Scheme::negate(Ciphertext& cipher) {
	uint64_t* axres = new uint64_t[cipher.curr_limbs << context.logN];
	uint64_t* bxres = new uint64_t[cipher.curr_limbs << context.logN];

	context.negate(axres, cipher.ax, cipher.curr_limbs);
	context.negate(bxres, cipher.bx, cipher.curr_limbs);

	return Ciphertext(axres, bxres, context.N, cipher.slots, cipher.curr_limbs);
}

void Scheme::negateAndEqual(Ciphertext& cipher) {
	long shift = 0;
	for (long i = 0; i < cipher.curr_limbs; ++i) {
		context.qiNegateAndEqual(cipher.ax + shift, i);
		context.qiNegateAndEqual(cipher.bx + shift, i);
		shift += context.N;
	}
}

Ciphertext Scheme::add(Ciphertext& cipher1, Ciphertext& cipher2) {
	uint64_t* axres = new uint64_t[cipher1.curr_limbs << context.logN];
	uint64_t* bxres = new uint64_t[cipher1.curr_limbs << context.logN];

	context.add(axres, cipher1.ax, cipher2.ax, cipher1.curr_limbs);
	context.add(bxres, cipher1.bx, cipher2.bx, cipher1.curr_limbs);

	return Ciphertext(axres, bxres, context.N, cipher1.slots, cipher1.curr_limbs);
}

void Scheme::addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	if(cipher1.curr_limbs != cipher2.curr_limbs) {
		throw invalid_argument("Ciphertexts are not comparable");
	}

	context.addAndEqual(cipher1.ax, cipher2.ax, cipher1.curr_limbs);
	context.addAndEqual(cipher1.bx, cipher2.bx, cipher1.curr_limbs);
}

Ciphertext Scheme::sub(Ciphertext& cipher1, Ciphertext& cipher2) {
	if(cipher1.curr_limbs != cipher2.curr_limbs) {
		throw invalid_argument("Ciphertexts are not comparable");
	}

	uint64_t* axres = new uint64_t[cipher1.curr_limbs << context.logN];
	uint64_t* bxres = new uint64_t[cipher1.curr_limbs << context.logN];

	context.sub(axres, cipher1.ax, cipher2.ax, cipher1.curr_limbs);
	context.sub(bxres, cipher1.bx, cipher2.bx, cipher1.curr_limbs);

	return Ciphertext(axres, bxres, context.N, cipher1.slots, cipher1.curr_limbs);
}

void Scheme::subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	if(cipher1.curr_limbs != cipher2.curr_limbs) {
		throw invalid_argument("Ciphertexts are not comparable");
	}

	context.subAndEqual(cipher1.ax, cipher2.ax, cipher1.curr_limbs);
	context.subAndEqual(cipher1.bx, cipher2.bx, cipher1.curr_limbs);
}

void Scheme::sub2AndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {
	if(cipher1.curr_limbs != cipher2.curr_limbs) {
		throw invalid_argument("Ciphertexts are not comparable");
	}

	context.sub2AndEqual(cipher1.ax, cipher2.ax, cipher1.curr_limbs);
	context.sub2AndEqual(cipher1.bx, cipher2.bx, cipher1.curr_limbs);
}

Ciphertext Scheme::mult(Ciphertext& cipher1, Ciphertext& cipher2) {
    uint64_t* axbx1 = new uint64_t[cipher1.curr_limbs << context.logN]();
    uint64_t* axbx2 = new uint64_t[cipher1.curr_limbs << context.logN]();

    uint64_t* axax = new uint64_t[cipher1.curr_limbs << context.logN]();
    uint64_t* bxbx = new uint64_t[cipher1.curr_limbs << context.logN]();

    context.add(axbx1, cipher1.ax, cipher1.bx, cipher1.curr_limbs);
    context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.curr_limbs);
    context.mulAndEqual(axbx1, axbx2, cipher1.curr_limbs);
    context.mul(bxbx, cipher1.bx, cipher2.bx, cipher1.curr_limbs);
    context.mul(axax, cipher1.ax, cipher2.ax, cipher1.curr_limbs);
    //只关心输入的level curr_limbs，默认转到 L+K个special modulus上
    //ApproxModUp
    context.raiseAndEqual(axax, cipher1.curr_limbs);

    //InnerProduct
    Key key = keyMap.at(MULTIPLICATION);

    uint64_t* axmult = new uint64_t[(cipher1.curr_limbs + context.K) << context.logN]();
    uint64_t* bxmult = new uint64_t[(cipher1.curr_limbs + context.K) << context.logN]();
    //evk是一个二维向量，在FullRNS-Mult第三步中，evk的两个分量分别和升模后的d2(j)相乘
    context.mulKey(axmult, axax, key.ax, cipher1.curr_limbs);
    context.mulKey(bxmult, axax, key.bx, cipher1.curr_limbs);

    //Mod Down
    context.backAndEqual(axmult, cipher1.curr_limbs);
    context.backAndEqual(bxmult, cipher1.curr_limbs);

    //凑数
    context.addAndEqual(axmult, axbx1, cipher1.curr_limbs);
    context.subAndEqual(axmult, bxbx, cipher1.curr_limbs);
    context.subAndEqual(axmult, axax, cipher1.curr_limbs);
    context.addAndEqual(bxmult, bxbx, cipher1.curr_limbs);

    delete[] axax;
    delete[] bxbx;
    delete[] axbx1;
    delete[] axbx2;

    return Ciphertext(axmult, bxmult, context.N, cipher1.slots, cipher1.curr_limbs);
}

void Scheme::multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2) {

	uint64_t* axbx1 = new uint64_t[cipher1.curr_limbs << context.logN]();
	uint64_t* axbx2 = new uint64_t[cipher1.curr_limbs << context.logN]();

	uint64_t* axax = new uint64_t[cipher1.curr_limbs << context.logN]();
	uint64_t* bxbx = new uint64_t[cipher1.curr_limbs << context.logN]();

	context.add(axbx1, cipher1.ax, cipher1.bx, cipher1.curr_limbs); // ax1 + bx1 mod P, 0 mod Q
	context.add(axbx2, cipher2.ax, cipher2.bx, cipher2.curr_limbs); // ax2 + bx2 mod P, 0 mod Q

	context.mulAndEqual(axbx1, axbx2, cipher1.curr_limbs); // (ax1 + bx1) * (ax2 + bx2) mod P, 0 mod Q

	context.mul(bxbx, cipher1.bx, cipher2.bx, cipher1.curr_limbs); // bx1 * bx2 mod P, 0 mod Q
	context.mul(axax, cipher1.ax, cipher2.ax, cipher1.curr_limbs); // ax1 * ax2 mod P, 0 mod Q
	context.raiseAndEqual(axax, cipher1.curr_limbs); // ax1 * ax2 mod P, ax1 * ax2 + e * P mod Q//Mod up

	Key key = keyMap.at(MULTIPLICATION); // kbx - kax * sx = (sxsx * Q + ex mod P, ex mod Q)

	delete[] cipher1.ax;
	delete[] cipher1.bx;

	cipher1.ax = new uint64_t[(cipher1.curr_limbs + context.K) << context.logN]();
	cipher1.bx = new uint64_t[(cipher1.curr_limbs + context.K) << context.logN]();

	context.mulKey(cipher1.ax, axax, key.ax, cipher1.curr_limbs);
	context.mulKey(cipher1.bx, axax, key.bx, cipher1.curr_limbs);

	context.backAndEqual(cipher1.ax, cipher1.curr_limbs);    //Mod Down
	context.backAndEqual(cipher1.bx, cipher1.curr_limbs);    //Mod Down

	context.addAndEqual(cipher1.ax, axbx1, cipher1.curr_limbs);  //凑数字
	context.subAndEqual(cipher1.ax, bxbx, cipher1.curr_limbs);
	context.subAndEqual(cipher1.ax, axax, cipher1.curr_limbs);
	context.addAndEqual(cipher1.bx, bxbx, cipher1.curr_limbs);

	delete[] axax;
	delete[] bxbx;
	delete[] axbx1;
	delete[] axbx2;
}


Ciphertext Scheme::square(Ciphertext& cipher) {
    uint64_t* axbx = new uint64_t[cipher.curr_limbs << context.logN]();

	uint64_t* axax = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bxbx = new uint64_t[cipher.curr_limbs << context.logN]();

	uint64_t* axmult = new uint64_t[(cipher.curr_limbs + context.K) << context.logN]();
	uint64_t* bxmult = new uint64_t[(cipher.curr_limbs + context.K) << context.logN]();

	context.add(axbx, cipher.ax, cipher.bx, cipher.curr_limbs); // ax1 + bx1 mod P, 0 mod Q

	context.squareAndEqual(axbx, cipher.curr_limbs); // (ax1 + bx1) * (ax2 + bx2) mod P, 0 mod Q
	context.square(bxbx, cipher.bx, cipher.curr_limbs);
	context.square(axax, cipher.ax, cipher.curr_limbs);

	context.raiseAndEqual(axax, cipher.curr_limbs);

	Key key = keyMap.at(MULTIPLICATION);

	context.mulKey(axmult, axax, key.ax, cipher.curr_limbs);
	context.mulKey(bxmult, axax, key.bx, cipher.curr_limbs);

	context.backAndEqual(axmult, cipher.curr_limbs);
	context.backAndEqual(bxmult, cipher.curr_limbs);

	context.addAndEqual(axmult, axbx, cipher.curr_limbs);
	context.subAndEqual(axmult, bxbx, cipher.curr_limbs);
	context.subAndEqual(axmult, axax, cipher.curr_limbs);
	context.addAndEqual(bxmult, bxbx, cipher.curr_limbs);

	delete[] axax;
	delete[] bxbx;
	delete[] axbx;

	return Ciphertext(axmult, bxmult, context.N, cipher.slots, cipher.curr_limbs);
}

void Scheme::squareAndEqual(Ciphertext& cipher) {
	uint64_t* axbx = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* axax = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bxbx = new uint64_t[cipher.curr_limbs << context.logN]();

	context.add(axbx, cipher.ax, cipher.bx, cipher.curr_limbs); // ax1 + bx1 mod P, 0 mod Q

	context.squareAndEqual(axbx, cipher.curr_limbs); // (ax1 + bx1) * (ax2 + bx2) mod P, 0 mod Q
	context.square(bxbx, cipher.bx, cipher.curr_limbs);
	context.square(axax, cipher.ax, cipher.curr_limbs);

	context.raiseAndEqual(axax, cipher.curr_limbs);

	delete[] cipher.ax;
	delete[] cipher.bx;

	cipher.ax = new uint64_t[(cipher.curr_limbs + context.K) << context.logN]();
	cipher.bx = new uint64_t[(cipher.curr_limbs + context.K) << context.logN]();

	Key key = keyMap.at(MULTIPLICATION);

	context.mulKey(cipher.ax, axax, key.ax, cipher.curr_limbs);
	context.mulKey(cipher.bx, axax, key.bx, cipher.curr_limbs);

	context.backAndEqual(cipher.ax, cipher.curr_limbs);
	context.backAndEqual(cipher.bx, cipher.curr_limbs);

	context.addAndEqual(cipher.ax, axbx, cipher.curr_limbs);
	context.subAndEqual(cipher.ax, bxbx, cipher.curr_limbs);
	context.subAndEqual(cipher.ax, axax, cipher.curr_limbs);
	context.addAndEqual(cipher.bx, bxbx, cipher.curr_limbs);

	delete[] axax;
	delete[] bxbx;
	delete[] axbx;
}

Ciphertext Scheme::imult(Ciphertext& cipher) {
    //TODO implement method
}

void Scheme::imultAndEqual(Ciphertext& cipher) {
	//TODO implement method
}

Ciphertext Scheme::idiv(Ciphertext& cipher) {
	//TODO implement method
	Ciphertext res;
	return res;

}

void Scheme::idivAndEqual(Ciphertext& cipher) {
	//TODO implement method
}

Ciphertext Scheme::addConst(Ciphertext& cipher, double cnst) {
	uint64_t tmpr = abs(cnst) * context.p;
	uint64_t* ax = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();
	copy(cipher.ax, cipher.ax + (cipher.curr_limbs << context.logN), ax);

	if(cnst >= 0) {
		context.addConst(bx, cipher.bx, tmpr, cipher.curr_limbs);
	} else {
		context.subConst(bx, cipher.bx, tmpr, cipher.curr_limbs);
	}

	return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);

}

Ciphertext Scheme::addConst(Ciphertext& cipher, complex<double> cnst) {
	//TODO implement method
	Ciphertext res;
	return res;

}

void Scheme::addConstAndEqual(Ciphertext& cipher, double cnst) {
	uint64_t tmpr = abs(cnst) * context.p;
	if(cnst >= 0) {
		context.addConstAndEqual(cipher.bx, tmpr, cipher.curr_limbs);
	} else {
		context.subConstAndEqual(cipher.bx, tmpr, cipher.curr_limbs);
	}
}

void Scheme::addConstAndEqual(Ciphertext& cipher, complex<double> cnst) {
	//TODO implement method
}

void Scheme::addPcAndEqual(Ciphertext& cipher) {
	context.addAndEqual(cipher.bx, context.pccoeff, cipher.curr_limbs);
}

void Scheme::addP2AndEqual(Ciphertext& cipher) {
	context.addAndEqual(cipher.bx, context.p2coeff, cipher.curr_limbs);
}

void Scheme::addP2hAndEqual(Ciphertext& cipher) {
	context.addAndEqual(cipher.bx, context.p2hcoeff, cipher.curr_limbs);
}

Ciphertext Scheme::multByConst(Ciphertext& cipher, double cnst) {
    if (DEBUG_SK=="NULL"){
        uint64_t tmpr = abs(cnst) * context.p;

        uint64_t* ax = new uint64_t[cipher.curr_limbs << context.logN]();
        uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();

        context.mulConst(ax, cipher.ax, tmpr, cipher.curr_limbs);
        context.mulConst(bx, cipher.bx, tmpr, cipher.curr_limbs);

        if(cnst < 0) {
            context.negateAndEqual(ax, cipher.curr_limbs);
            context.negateAndEqual(bx, cipher.curr_limbs);
        }
        return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);
    }
    else if(DEBUG_SK=="CURR_SK"){
        //fixme: //Debug Only 暂时忽视对 负整数、正整数的不支持，直接用解密再加密保证计算正确性来完成上层debug
        // "Check: Ciphertext Scheme::multByConst(Ciphertext& cipher, double cnst), if computation failed!"
        complex<double> *dvecCMult = decrypt(secretKey, cipher);
        for (int i = 0; i < cipher.slots; ++i) {
            dvecCMult[i]*=cnst;
        }

        Ciphertext res = encrypt(dvecCMult, cipher.slots, cipher.curr_limbs);
        context.mulConstAndEqual( res.ax, context.p, cipher.curr_limbs);
        context.mulConstAndEqual( res.bx, context.p, cipher.curr_limbs);
        return res;
    }
}

Ciphertext Scheme::multByConst(Ciphertext& cipher, complex<double> cnst) {
    uint64_t tmpr = abs(cnst.real()) * context.p;
    uint64_t tmpi = abs(cnst.imag()) * context.p;

    uint64_t* ax = new uint64_t[cipher.curr_limbs << context.logN]();
    uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();
    uint64_t* axi = new uint64_t[cipher.curr_limbs << context.logN]();
    uint64_t* bxi = new uint64_t[cipher.curr_limbs << context.logN]();

    //法一：//[Coeff]下实现ct_ax * X^{N/2}, ct_bx * X^{N/2}
//    context.INTTAndEqual(cipher.ax, cipher.curr_limbs);
//    context.INTTAndEqual(cipher.bx, cipher.curr_limbs);
//	  context.mulByMonomial(axi, cipher.ax, cipher.curr_limbs, context.Nh);
//	  context.mulByMonomial(bxi, cipher.bx, cipher.curr_limbs, context.Nh);
//    context.NTTAndEqual(cipher.ax, cipher.curr_limbs);
//    context.NTTAndEqual(cipher.bx, cipher.curr_limbs);
//    context.NTTAndEqual(axi, cipher.curr_limbs);
//    context.NTTAndEqual(bxi, cipher.curr_limbs);

    //法二：//[Eval]下实现ct_ax * X^{N/2}, ct_bx * X^{N/2}
    SetZero(axi,0,cipher.curr_limbs << context.logN);
    SetZero(bxi,0,cipher.curr_limbs << context.logN);
    for (int i = 0; i < cipher.curr_limbs; ++i) {
        int index=context.Nh+(i<< context.logN);
        axi[index]=1;
        bxi[index]=1;
    }
    context.NTTAndEqual(axi, cipher.curr_limbs);
    context.NTTAndEqual(bxi, cipher.curr_limbs);
    context.mulAndEqual(axi, cipher.ax,  cipher.curr_limbs);
    context.mulAndEqual(bxi, cipher.bx,  cipher.curr_limbs);

    //Δ * abs(cnst.real) * ct_ax, Δ * abs(cnst.real) * ct_bx
    context.mulConst(ax, cipher.ax, tmpr, cipher.curr_limbs);
    context.mulConst(bx, cipher.bx, tmpr, cipher.curr_limbs);

    //Δ * abs(cnst.imag) * ct_bx, Δ * abs(cnst.imag) * ct_bx
    context.mulConstAndEqual(axi, tmpi, cipher.curr_limbs);
    context.mulConstAndEqual(bxi, tmpi, cipher.curr_limbs);

    if (cnst.real()<0){
        context.negateAndEqual(ax, cipher.curr_limbs);
        context.negateAndEqual(bx, cipher.curr_limbs);
    }
    if (cnst.imag()<0){
        context.negateAndEqual(axi, cipher.curr_limbs);
        context.negateAndEqual(bxi, cipher.curr_limbs);
    }

    context.addAndEqual(ax, axi, cipher.curr_limbs);
    context.addAndEqual(bx, bxi, cipher.curr_limbs);

    return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);
}

Ciphertext Scheme::multByConstVec(Ciphertext& cipher, double* cnstVec, long slots) {
    cout<<"multByConstVec needs to encode a double vector, which is not implemented!"<<endl;
    uint64_t* ax = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* mx = new uint64_t[cipher.curr_limbs << context.logN]();

    //FIXME: Encode a double vector `cnstVec` is not implemented!
	context.encode(mx, cnstVec, slots, cipher.curr_limbs);
	context.mul(ax, cipher.ax, mx, cipher.curr_limbs);
	context.mul(bx, cipher.bx, mx, cipher.curr_limbs);
	delete[] mx;
	return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);
}

Ciphertext Scheme::multByConstVec(Ciphertext& cipher, complex<double>* cnstVec, long slots) {
	uint64_t* ax = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* mx = new uint64_t[cipher.curr_limbs << context.logN]();

	context.encode(mx, cnstVec, slots, cipher.curr_limbs);
	context.mul(ax, cipher.ax, mx, cipher.curr_limbs);
	context.mul(bx, cipher.bx, mx, cipher.curr_limbs);
	delete[] mx;
	return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);
}

void Scheme::multByConstAndEqual(Ciphertext& cipher, double cnst) {
	uint64_t tmpr = abs(cnst) * context.p;

	context.mulConstAndEqual(cipher.ax, tmpr, cipher.curr_limbs);
	context.mulConstAndEqual(cipher.bx, tmpr, cipher.curr_limbs);

	if(cnst < 0) {
		context.negateAndEqual(cipher.ax, cipher.curr_limbs);
		context.negateAndEqual(cipher.bx, cipher.curr_limbs);
	}
}

void Scheme::multByConstAndEqual(Ciphertext& cipher, complex<double> cnst) {
}

void Scheme::multByPolyAndEqual(Ciphertext& cipher, uint64_t* poly) {
	context.mulAndEqual(cipher.ax, poly, cipher.curr_limbs);
	context.mulAndEqual(cipher.bx, poly, cipher.curr_limbs);
}

void Scheme::multByMonomial(Ciphertext &result, Ciphertext &cipher, long monomialDegree) {

    long l=cipher.curr_limbs;

    context.INTTAndEqual(cipher.ax, cipher.curr_limbs);
    context.INTTAndEqual(cipher.bx, cipher.curr_limbs);

    //mulByMonomial should be done in coeff mode!
    context.mulByMonomial(result.ax, cipher.ax, l, monomialDegree);
    context.mulByMonomial(result.bx, cipher.bx, l, monomialDegree);

    context.NTTAndEqual(cipher.ax, cipher.curr_limbs);
    context.NTTAndEqual(cipher.bx, cipher.curr_limbs);
    context.NTTAndEqual(result.ax, cipher.curr_limbs);
    context.NTTAndEqual(result.bx, cipher.curr_limbs);
}

void Scheme::multByMonomialAndEqual(Ciphertext& cipher, long monomialDegree) {

    long l=cipher.curr_limbs;

    context.INTTAndEqual(cipher.ax, cipher.curr_limbs);
    context.INTTAndEqual(cipher.bx, cipher.curr_limbs);

    //mulByMonomial should be done in coeff mode!
    context.mulByMonomialAndEqual(cipher.ax, l, monomialDegree);
    context.mulByMonomialAndEqual(cipher.bx, l, monomialDegree);

    context.NTTAndEqual(cipher.ax, cipher.curr_limbs);
    context.NTTAndEqual(cipher.bx, cipher.curr_limbs);
}


Ciphertext Scheme::reScaleBy(Ciphertext& cipher, long dl) {
	//TODO implement method
	Ciphertext res;
	return res;
}

void Scheme::reScaleByAndEqual(Ciphertext& cipher, long dl) {
	for (long i = 0; i < dl; ++i) {
		context.reScaleAndEqual(cipher.ax, cipher.curr_limbs);
		context.reScaleAndEqual(cipher.bx, cipher.curr_limbs);
		cipher.curr_limbs -= 1;
	}
}

Ciphertext Scheme::reScaleTo(Ciphertext& cipher, long l) {
	long dl = cipher.curr_limbs - l;
	return reScaleBy(cipher, dl);
}

void Scheme::reScaleToAndEqual(Ciphertext& cipher, long l) {
	long dl = cipher.curr_limbs - l;
	reScaleByAndEqual(cipher, dl);
}

Ciphertext Scheme::modDownBy(Ciphertext& cipher, long dl) {
	uint64_t* ax = context.modDown(cipher.ax, cipher.curr_limbs, dl);
	uint64_t* bx = context.modDown(cipher.bx, cipher.curr_limbs, dl);
	return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs - dl);
}

void Scheme::modDownByAndEqual(Ciphertext& cipher, long dl) {
	context.modDownAndEqual(cipher.ax, cipher.curr_limbs, dl);
	context.modDownAndEqual(cipher.bx, cipher.curr_limbs, dl);
	cipher.curr_limbs -= dl;
}

Ciphertext Scheme::modDownTo(Ciphertext& cipher, long l) {
	long dl = cipher.curr_limbs - l;
	return modDownBy(cipher, dl);
}

void Scheme::modDownToAndEqual(Ciphertext& cipher, long l) {
	long dl = cipher.curr_limbs - l;
	modDownByAndEqual(cipher, dl);
}

Ciphertext Scheme::leftRotateFast(Ciphertext& cipher, long rotSlots) {
	uint64_t* bxrot = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();

	uint64_t* ax = new uint64_t[(cipher.curr_limbs + context.K) << context.logN]();

	Key key = leftRotKeyMap.at(rotSlots);

	context.leftRot(bxrot, cipher.bx, cipher.curr_limbs, rotSlots);
	context.leftRot(bx, cipher.ax, cipher.curr_limbs, rotSlots);

	context.raiseAndEqual(bx, cipher.curr_limbs);

	context.mulKey(ax, bx, key.ax, cipher.curr_limbs);
	context.mulKey(bx, bx, key.bx, cipher.curr_limbs);

	context.backAndEqual(ax, cipher.curr_limbs);
	context.backAndEqual(bx, cipher.curr_limbs);

	context.addAndEqual(bx, bxrot, cipher.curr_limbs);

	return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);
}

void Scheme::leftRotateAndEqualFast(Ciphertext& cipher, long rotSlots) {
	uint64_t* bxrot = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* bx = new uint64_t[cipher.curr_limbs << context.logN]();
	uint64_t* ax = new uint64_t[(cipher.curr_limbs + context.K) << context.logN]();

	context.leftRot(bxrot, cipher.bx, cipher.curr_limbs, rotSlots);
	context.leftRot(bx, cipher.ax, cipher.curr_limbs, rotSlots);

	Key key = leftRotKeyMap.at(rotSlots);

	context.raiseAndEqual(bx, cipher.curr_limbs);

	context.mulKey(ax, bx, key.ax, cipher.curr_limbs);
	context.mulKey(bx, bx, key.bx, cipher.curr_limbs);

	context.back(cipher.ax, ax, cipher.curr_limbs);
	context.back(cipher.bx, bx, cipher.curr_limbs);

	context.addAndEqual(cipher.bx, bxrot, cipher.curr_limbs);
}

Ciphertext Scheme::leftRotateByPo2(Ciphertext& cipher, long logRotSlots) {
	long rotSlots = (1 << logRotSlots);
	return leftRotateFast(cipher, rotSlots);
}

void Scheme::leftRotateByPo2AndEqual(Ciphertext& cipher, long logRotSlots) {
	long rotSlots = (1 << logRotSlots);
	leftRotateAndEqualFast(cipher, rotSlots);
}

Ciphertext Scheme::rightRotateByPo2(Ciphertext& cipher, long logRotSlots) {
	long rotSlots = context.Nh - (1 << logRotSlots);
	return leftRotateFast(cipher, rotSlots);
}

void Scheme::rightRotateByPo2AndEqual(Ciphertext& cipher, long logRotSlots) {
	long rotSlots = context.Nh - (1 << logRotSlots);
	leftRotateAndEqualFast(cipher, rotSlots);
}

Ciphertext Scheme::leftRotate(Ciphertext& cipher, long rotSlots) {
	Ciphertext res = cipher;
	leftRotateAndEqual(res, rotSlots);
	return res;
}

void Scheme::leftRotateAndEqual(Ciphertext& cipher, long rotSlots) {
	long remrotSlots = rotSlots % cipher.slots;
	long logrotSlots = log2((double)remrotSlots) + 1;
	for (long i = 0; i < logrotSlots; ++i) {
		if(remrotSlots & 1 << i) {
			leftRotateByPo2AndEqual(cipher, i);
		}
	}
}

Ciphertext Scheme::rightRotate(Ciphertext& cipher, long rotSlots) {
	Ciphertext res = cipher;
	rightRotateAndEqual(res, rotSlots);
	return res;
}

void Scheme::rightRotateAndEqual(Ciphertext& cipher, long rotSlots) {
	long remrotSlots = rotSlots % cipher.slots;
	long logrotSlots = log2((double)remrotSlots) + 1;
	for (long i = 0; i < logrotSlots; ++i) {
		if(remrotSlots & 1 << i) {
			rightRotateByPo2AndEqual(cipher, i);
		}
	}
}

Ciphertext Scheme::conjugate(Ciphertext& cipher) {
	uint64_t dnum=context.dnum;
    uint64_t* bxconj = new uint64_t[context.N * cipher.curr_limbs];
	uint64_t* bx = new uint64_t[context.N * cipher.curr_limbs];
	uint64_t* ax = new uint64_t[context.N * (cipher.curr_limbs + context.K)];

	context.conjugate(bxconj, cipher.bx, cipher.curr_limbs);
	context.conjugate(bx, cipher.ax, cipher.curr_limbs);

	Key key = keyMap.at(CONJUGATION*dnum);

	context.raiseAndEqual(bx, cipher.curr_limbs);

	long shift = 0;
	for (long i = 0; i < cipher.curr_limbs; ++i) {
		context.qiMul(ax + shift, bx + shift, key.ax + shift, i);
		context.qiMulAndEqual(bx + shift, key.bx + shift, i);
		shift += context.N;
	}

	long msshift = context.N * context.L;
	for (long i = 0; i < context.K; ++i) {
		context.piMul(ax + shift, bx + shift, key.ax + msshift, i);
		context.piMulAndEqual(bx + shift, key.bx + msshift, i);
		shift += context.N;
		msshift += context.N;
	}

	context.backAndEqual(ax, cipher.curr_limbs);
	context.backAndEqual(bx, cipher.curr_limbs);

	context.addAndEqual(bx, bxconj, cipher.curr_limbs);

	delete[] bxconj;

	return Ciphertext(ax, bx, context.N, cipher.slots, cipher.curr_limbs);
}

void Scheme::conjugateAndEqual(Ciphertext& cipher) {
	uint64_t* bxconj = new uint64_t[context.N * cipher.curr_limbs];
	uint64_t* bx = new uint64_t[context.N * cipher.curr_limbs];
	uint64_t* ax = new uint64_t[context.N * (cipher.curr_limbs + context.K)];

	context.conjugate(bxconj, cipher.bx, cipher.curr_limbs);
	context.conjugate(bx, cipher.ax, cipher.curr_limbs);

	Key key = keyMap.at(CONJUGATION);

	context.raiseAndEqual(bx, cipher.curr_limbs);

	long shift = 0;
	for (long i = 0; i < cipher.curr_limbs; ++i) {
		context.qiMul(ax + shift, bx + shift, key.ax + shift, i);
		context.qiMulAndEqual(bx + shift, key.bx + shift, i);
		shift += context.N;
	}

	long msshift = context.N * context.L;
	for (long i = 0; i < context.K; ++i) {
		context.piMul(ax + shift, bx + shift, key.ax + msshift, i);
		context.piMulAndEqual(bx + shift, key.bx + msshift, i);
		shift += context.N;
		msshift += context.N;
	}

	context.back(cipher.ax, ax, cipher.curr_limbs);
	context.back(cipher.bx, bx, cipher.curr_limbs);

	context.addAndEqual(cipher.bx, bxconj, cipher.curr_limbs);

	delete[] bxconj;
}

// The advanced method to compute the Ciphertext in mult a const
Ciphertext Scheme::Lattigo_MultByConst(Ciphertext& in, complex<double> cnst)
{
    auto l = in.curr_limbs - 1;

    double cReal, cImag;
    long scale;
    scale = 1;
    cReal = cnst.real();
    cImag = cnst.imag();
    if (cReal != 0)
    {
        auto valueInt = int(cReal);
        auto valueFloat = cReal - double(valueInt);

        if (valueFloat != 0)
        {
            scale = context.p;
        }
    }
    if (cImag != 0)
    {
        auto valueInt = int(cImag);
        auto valueFloat = cImag - double(valueInt);

        if (valueFloat != 0)
        {
            scale = context.p;
        }
    }

    //auto context = Scheme::context;
    uint64_t scaledConst, scaledConstReal, scaledConstImag;
    uint64_t* ax = new uint64_t[in.curr_limbs << context.logN]();
    uint64_t* bx = new uint64_t[in.curr_limbs << context.logN]();
    for (uint64_t i = 0; i < l+1; i++)
    {
        uint64_t qi = context.qVec[i];
        uint64_t bredParams[] = {0, 0};
        Context::ComputeBRedParameters(qi, bredParams[0], bredParams[1]);
        uint64_t mredParams = context.qInvVec[i];

        scaledConstReal = 0, scaledConstImag = 0, scaledConst = 0;

        if (cReal != 0)
        {
            scaledConstReal = Context::ScaleUpExact(cReal, scale, qi);
            scaledConst = scaledConstReal;
        }

        if (cImag != 0)
        {
            scaledConstImag = Context::ScaleUpExact(cImag, scale, qi);
            scaledConstImag = Context::MRed(scaledConstImag, context.nttPsi[i][1], qi, mredParams);
            scaledConst = Context::CRed(scaledConst + scaledConstImag, qi);
        }

        scaledConst = Context::MForm(scaledConst, qi, bredParams);

        for (long j = 0; j < (1<<(context.logN-1)); ++j)
        {
            uint64_t* ai = in.ax + (i << context.logN);
            uint64_t* bi = in.bx + (i << context.logN);
            uint64_t* resa = ax + (i <<context.logN);
            uint64_t* resb = bx + (i <<context.logN);
            resa[j] = Context::MRed(ai[j], scaledConst, qi, mredParams);
            resb[j] = Context::MRed(bi[j], scaledConst, qi, mredParams);
        }

        if (cImag != 0)
        {
            scaledConst = Context::CRed(scaledConstReal + (qi - scaledConstImag), qi);
            scaledConst = Context::MForm(scaledConst, qi, bredParams);
        }

        for (long j = (1<<(context.logN-1)); j < (1<<context.logN); ++j)
        {
            uint64_t* ai = in.ax + (i << context.logN);
            uint64_t* bi = in.bx + (i << context.logN);
            uint64_t* resa = ax + (i <<context.logN);
            uint64_t* resb = bx + (i <<context.logN);
            resa[j] = Context::MRed(ai[j], scaledConst, qi, mredParams);
            resb[j] = Context::MRed(bi[j], scaledConst, qi, mredParams);
        }
    }
    return Ciphertext(ax, bx, context.N, in.slots, in.curr_limbs);
}

complex<double>* Scheme::Lattigo_decrypt(SecretKey& secretKey, Ciphertext& cipher) {
    uint64_t* ax = new uint64_t[cipher.curr_limbs << context.logN]();
    auto level = cipher.curr_limbs - 1;
    auto scale = context.p;
    for (int i = 0; i < level + 1; ++i)
    {
        for (uint64_t j = 0; j < (1<<(context.logN)); ++j)
        {
            uint64_t* ai = ax + (i << context.logN);
            uint64_t* bi = cipher.bx + (i << context.logN);
            ai[j] = bi[j];
        }
    }
    // 1 = Cipher.Degree()
    for (int m = 1; m > 0; m--)
    {
        for (uint64_t n = 0; n < level + 1; n++)
        {   // bx * s
            auto qi = context.qVec[n];
            auto p1 = ax + (n << context.logN);
            auto p2 = secretKey.sx + (n << context.logN);
            auto p3 = ax + (n << context.logN);
            auto mredParams = context.qInvVec[n];
            for (uint64_t k = 0; k < (1<<(context.logN)); ++k)
            {
                p3[k] = Context::MRed(p1[k], p2[k], qi, mredParams);
            }
        }

        for (uint64_t n = 0; n < level + 1; n++)
        {   // bx * s + ax
            auto qi = context.qVec[n];
            auto p1 = ax + (n << context.logN);
            auto p2 = cipher.ax + (n << context.logN);
            auto p3 = ax + (n << context.logN);
            for (uint64_t k = 0; k < (1<<(context.logN)); ++k)
            {
                p3[k] = Context::CRed(p1[k] + p2[k], qi);
            }
        }

        if ((m&7) == 7)
        {
            //ReduceLvl();...
            Scheme::ReduceLvl(level, ax, ax);
        }
    }
    //----Decryption Done!

    //return decode(msg);
}

void Scheme::ReduceLvl(uint64_t level, uint64_t* p1, uint64_t* p2)
{
    uint64_t qi;
    for (uint64_t n = 0; n < level + 1; n++)
    {   // bx * s
        qi = context.qVec[n];
        auto p1tmp = p1 + (n << context.logN);
        auto p2tmp = p2 + (n << context.logN);
        uint64_t bredParams[] = {0, 0};
        Context::ComputeBRedParameters(qi, bredParams[0], bredParams[1]);
        for (uint64_t k = 0; k < (1<<(context.logN)); ++k)
        {
            p2tmp[k] = Context::BRedAdd(p1[k], qi, bredParams);
        }
    }
}

complex<double>* Scheme::fulldecrypt(SecretKey& secretKey, Ciphertext& cipher) {
    Plaintext msg = fulldecryptMsg(secretKey, cipher);
    return fulldecode(msg);
}

Plaintext Scheme::fulldecryptMsg(SecretKey& secretKey, Ciphertext& cipher) {
    uint64_t* mx = new uint64_t[cipher.curr_limbs << context.logN]();
    context.mul(mx, cipher.ax, secretKey.sx, cipher.curr_limbs);
    context.addAndEqual(mx, cipher.bx, cipher.curr_limbs);

    return Plaintext(mx, context.N, cipher.slots, cipher.curr_limbs);
}

complex<double>* Scheme::fulldecode(Plaintext& msg) {
    auto res = new complex<double>[msg.slots]();
    context.fulldecode(msg.mx, res, msg.slots, msg.l);
    return res;
}

Plaintext Scheme::fullencode(complex<double>* v, long slots, long l) {
    uint64_t* m = new uint64_t[l << context.logN]();
    context.fullencode(m, v, slots, l);
    return Plaintext(m, context.N, slots, l);
}

Ciphertext Scheme::fullencrypt(complex<double>* vals, long slots, long l) {
    Plaintext msg = fullencode(vals, slots, l);
    return encryptMsg(msg);
}





