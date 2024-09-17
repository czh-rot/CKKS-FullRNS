/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <stdint.h>

#include "Context.h"
#include "EvaluatorUtils.h"
#include "../data/testConstValue.h"
#include <boost/multiprecision/cpp_int.hpp>

Context::Context(long logN, long logp, long L, long K, long dnum, long h, double sigma) :
		logN(logN), logp(logp), L(L), K(K), dnum(dnum), h(h), sigma(sigma) {

	N = 1L << logN;
	M = N << 1;
	logNh = logN - 1;
	Nh = N >> 1;
	p = 1L << logp;


	qVec = new uint64_t[L]();
	qrVec = new uint64_t[L]();
	qTwok = new long[L]();
	qkVec = new uint64_t[L]();
	qdVec = new uint64_t[L]();
	qInvVec = new uint64_t[L]();
	qRoots = new uint64_t[L]();
	qRootsInv = new uint64_t[L]();
	qRootPows = new uint64_t*[L];
	qRootScalePows = new uint64_t*[L];
	qRootScalePowsOverq = new uint64_t*[L];
	qRootScalePowsInv = new uint64_t*[L];
	qRootPowsInv = new uint64_t*[L];
	NInvModq = new uint64_t[L]();
	NScaleInvModq = new uint64_t[L]();

	// Generate Primes //
	long bnd = 1;
	long cnt = 1;

	bnd = 1;
	while(1) {
			uint64_t prime = (1ULL << Q0_BIT_SIZE) + bnd * M + 1;
			if(primeTest(prime)) {
				qVec[0] = prime;
				break;
			}
			bnd++;
	}

	bnd = 1;
	while(cnt < L) {
		uint64_t prime1 = (1ULL << logp) + bnd * M + 1;
		if(primeTest(prime1)) {
			qVec[cnt] = prime1;
			cnt++;
		}
		uint64_t prime2 = (1ULL << logp) - bnd * M + 1;
		if(primeTest(prime2)) {
			qVec[cnt] = prime2;
			cnt++;
		}
		bnd++;
	}

	if(logp - logN - 1 - ceil(log2(bnd)) < 10) {
		cerr << "ERROR: too small number of precision" << endl;
		cerr << "TRY to use larger logp or smaller depth" << endl;
	}

	for (long i = 0; i < L; ++i) {
		qTwok[i] = (2 * ((long)log2(qVec[i]) + 1)); //ModBarrett 算法里的 k，详见笔记
		qrVec[i] = (static_cast<unsigned __int128>(1) << qTwok[i]) / qVec[i];   //ModBarret算法里的 m，详见笔记
		qkVec[i] = static_cast<uint64_t>(((static_cast<unsigned __int128>(invMod(((uint64_t)(1) << 62), qVec[i])) << 62) - 1) / qVec[i]);
		qdVec[i] = qVec[i] << 1;
		qRoots[i] = findMthRootOfUnity(M, qVec[i]);
		qRootsInv[i] = invMod(qRoots[i], qVec[i]);
		NInvModq[i] = invMod(N, qVec[i]);
		mulMod(NScaleInvModq[i], NInvModq[i], (static_cast<uint64_t>(1) << 32), qVec[i]);
		mulMod(NScaleInvModq[i], NScaleInvModq[i], (static_cast<uint64_t>(1) << 32), qVec[i]);
		qInvVec[i] = inv(qVec[i]);
		qRootPows[i] = new uint64_t[N]();
		qRootPowsInv[i] = new uint64_t[N]();
		qRootScalePows[i] = new uint64_t[N]();
		qRootScalePowsOverq[i] = new uint64_t[N]();
		qRootScalePowsInv[i] = new uint64_t[N]();
		uint64_t power = static_cast<uint64_t>(1);
		uint64_t powerInv = static_cast<uint64_t>(1);

		for (long j = 0; j < N; ++j) {
			uint64_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN);
			qRootPows[i][jprime] = power;
			unsigned __int128 tmp = (static_cast<unsigned __int128>(power) << 64);
			qRootScalePowsOverq[i][jprime] = static_cast<uint64_t>(tmp / qVec[i]);
			mulMod(qRootScalePows[i][jprime], qRootPows[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);
			mulMod(qRootScalePows[i][jprime], qRootScalePows[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);
			qRootPowsInv[i][jprime] = powerInv;
			mulMod(qRootScalePowsInv[i][jprime], qRootPowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);
			mulMod(qRootScalePowsInv[i][jprime], qRootScalePowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), qVec[i]);

			if (j < N - 1) {
				mulMod(power, power, qRoots[i], qVec[i]);
				mulMod(powerInv, powerInv, qRootsInv[i], qVec[i]);
			}
		}
	}

	pVec = new uint64_t[K]();
	prVec = new uint64_t[K]();
	pTwok = new long[K]();
	pkVec = new uint64_t[K]();
	pdVec = new uint64_t[K]();
	pInvVec = new uint64_t[K]();
	pRoots = new uint64_t[K]();
	pRootsInv = new uint64_t[K]();
	pRootPows = new uint64_t*[K];
	pRootPowsInv = new uint64_t*[K];
	pRootScalePows = new uint64_t*[K];
	pRootScalePowsOverp = new uint64_t*[K];
	pRootScalePowsInv = new uint64_t*[K];
	NInvModp = new uint64_t[K]();
	NScaleInvModp = new uint64_t[K]();

	// Generate Special Primes //
	cnt = 0;
	while(cnt < K) {
		uint64_t prime1 = (1ULL << logp) + bnd * M + 1;
		if(primeTest(prime1)) {
			pVec[cnt] = prime1;
			cnt++;
		}
		if(cnt == K) break;
		uint64_t prime2 = (1ULL << logp) - bnd * M + 1;
		if(primeTest(prime2)) {
			pVec[cnt] = prime2;
			cnt++;
		}
		bnd++;
	}

	for (long i = 0; i < K; ++i) {
		pTwok[i] = (2 * ((long)log2(pVec[i]) + 1));
		prVec[i] = (static_cast<unsigned __int128>(1) << pTwok[i]) / pVec[i];
		pkVec[i] = static_cast<uint64_t>(((static_cast<unsigned __int128>(invMod(((uint64_t)(1) << 62), pVec[i])) << 62) - 1) / pVec[i]);
		pdVec[i] = pVec[i] << 1;
		pRoots[i] = findMthRootOfUnity(M, pVec[i]);
		pRootsInv[i] = invMod(pRoots[i], pVec[i]);
		NInvModp[i] = invMod(N, pVec[i]);
		mulMod(NScaleInvModp[i], NInvModp[i], (static_cast<uint64_t>(1) << 32), pVec[i]);
		mulMod(NScaleInvModp[i], NScaleInvModp[i], (static_cast<uint64_t>(1) << 32), pVec[i]);
		pRootPows[i] = new uint64_t[N]();
		pRootScalePows[i] = new uint64_t[N]();
		pRootScalePowsOverp[i] = new uint64_t[N]();
		pRootScalePowsInv[i] = new uint64_t[N]();
		pRootPowsInv[i] = new uint64_t[N]();
		pInvVec[i] = inv(pVec[i]);
		uint64_t power = static_cast<uint64_t>(1);
		uint64_t powerInv = static_cast<uint64_t>(1);
		for (long j = 0; j < N; ++j) {
			uint64_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN);
			pRootPows[i][jprime] = power;
			unsigned __int128 tmp = (static_cast<unsigned __int128>(power) << 64);
			mulMod(pRootScalePows[i][jprime], pRootPows[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
			mulMod(pRootScalePows[i][jprime], pRootScalePows[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
			pRootPowsInv[i][jprime] = powerInv;
			mulMod(pRootScalePowsInv[i][jprime], pRootPowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
			mulMod(pRootScalePowsInv[i][jprime], pRootScalePowsInv[i][jprime], (static_cast<uint64_t>(1) << 32), pVec[i]);
			if (j < N - 1) {
				mulMod(power, power, pRoots[i], pVec[i]);
				mulMod(powerInv, powerInv, pRootsInv[i], pVec[i]);
			}
		}
	}

	qHatModq = new uint64_t*[L]; // [curr_limbs][i] (phat_i)_l mod p_i
	for (long l = 0; l < L; ++l) {
		qHatModq[l] = new uint64_t[l + 1]();
		for (long i = 0; i < l + 1; ++i) {
			qHatModq[l][i] = 1;
			for (long j = 0; j < i; ++j) {
				uint64_t temp = qVec[j] % qVec[i];
				mulMod(qHatModq[l][i], qHatModq[l][i], temp, qVec[i]);
			}
			for (long j = i + 1; j < l + 1; ++j) {
				uint64_t temp = qVec[j] % qVec[i];
				mulMod(qHatModq[l][i], qHatModq[l][i], temp, qVec[i]);
			}
		}
	}

	qHatInvModq = new uint64_t*[L]; // [curr_limbs][i] (phat_i)_l^-1 mod p_i
	for (long l = 0; l < L; ++l) {
		qHatInvModq[l] = new uint64_t[l + 1]();
		for (long i = 0; i < l + 1; ++i) {
			qHatInvModq[l][i] = invMod(qHatModq[l][i], qVec[i]);
		}
	}

	pHatModp = new uint64_t[K](); // [k] qhat_k mod q_k

	for (long k = 0; k < K; ++k) {
		pHatModp[k] = 1;
		for (long j = 0; j < k; ++j) {
			uint64_t temp = pVec[j] % pVec[k];
			mulMod(pHatModp[k], pHatModp[k], temp, pVec[k]);
		}
		for (long j = k + 1; j < K; ++j) {
			uint64_t temp = pVec[j] % pVec[k];
			mulMod(pHatModp[k], pHatModp[k], temp, pVec[k]);
		}
	}

	pHatInvModp = new uint64_t[K](); // [k] qhat_k^-1 mod q_k
	for (long k = 0; k < K; ++k) {
		pHatInvModp[k] = invMod(pHatModp[k], pVec[k]);
	}

	qHatModp = new uint64_t**[L];  // [curr_limbs] [i] [k]  (phat_i)_l mod q_k
	for (long l = 0; l < L; ++l) {
		qHatModp[l] = new uint64_t*[l + 1];
		for (long i = 0; i < l + 1; ++i) {
			qHatModp[l][i] = new uint64_t[K]();
			for (long k = 0; k < K; ++k) {
				qHatModp[l][i][k] = 1;
				for (long j = 0; j < i; ++j) {
					uint64_t temp = qVec[j] % pVec[k];
					mulMod(qHatModp[l][i][k], qHatModp[l][i][k], temp, pVec[k]);
				}
				for (long j = i + 1; j < l + 1; ++j) {
					uint64_t temp = qVec[j] % pVec[k];
					mulMod(qHatModp[l][i][k], qHatModp[l][i][k], temp, pVec[k]);
				}
			}
		}
	}

	pHatModq = new uint64_t*[K]; // [k][i] qhat_k mod p_i
	for (long k = 0; k < K; ++k) {
		pHatModq[k] = new uint64_t[L]();
		for (long i = 0; i < L; ++i) {
			pHatModq[k][i] = 1;
			for (long s = 0; s < k; ++s) {
				uint64_t temp = pVec[s] % qVec[i];
				mulMod(pHatModq[k][i], pHatModq[k][i], temp, qVec[i]);
			}
			for (long s = k + 1; s < K; ++s) {
				uint64_t temp = pVec[s] % qVec[i];
				mulMod(pHatModq[k][i], pHatModq[k][i], temp, qVec[i]);
			}
		}
	}



	PModq = new uint64_t[L](); // [i] qprod mod p_i
	for (long i = 0; i < L; ++i) {
		PModq[i] = 1;
		for (long k = 0; k < K; ++k) {
			uint64_t temp = pVec[k] % qVec[i];
			mulMod(PModq[i], PModq[i], temp, qVec[i]);
		}
	}

	PInvModq = new uint64_t[L](); // [i] qprod^-1 mod p_i
	for (long i = 0; i < L; ++i) {
		PInvModq[i] = invMod(PModq[i], qVec[i]);
	}

	QModp = new uint64_t*[L];
	for (long i = 0; i < L; ++i) {
		QModp[i] = new uint64_t[K]();
		for (long k = 0; k < K; ++k) {
			QModp[i][k] = 1;
			for (long j = 0; j < i + 1; ++j) {
				uint64_t temp = qVec[j] % pVec[k];
				mulMod(QModp[i][k], QModp[i][k], temp, pVec[k]);
			}
		}
	}
	QInvModp = new uint64_t*[L];
	for (long i = 0; i < L; ++i) {
		QInvModp[i] = new uint64_t[K]();
		for (long k = 0; k < K; ++k) {
			QInvModp[i][k] = invMod(QModp[i][k], pVec[k]);
		}
	}

	qInvModq = new uint64_t*[L]; // [i][j] p_i^-1 mod p_j
	for (long i = 0; i < L; ++i) {
		qInvModq[i] = new uint64_t[L]();
		for (long j = 0; j < i; ++j) {
			qInvModq[i][j] = invMod(qVec[i], qVec[j]);
		}
		for (long j = i + 1; j < L; ++j) {
			qInvModq[i][j] = invMod(qVec[i], qVec[j]);
		}
	}

	rotGroup = new long[Nh]();
	long fivePows = 1;
	for (long i = 0; i < Nh; ++i) {
		rotGroup[i] = fivePows;
		fivePows *= 5;
		fivePows %= M;
	}

	ksiPows = new complex<double>[M + 1];
	for (long j = 0; j < M; ++j) {
		double angle = 2.0 * M_PI * j / M;
		ksiPows[j].real(cos(angle));
		ksiPows[j].imag(sin(angle));
	}

	ksiPows[M] = ksiPows[0];

	p2coeff = new uint64_t[L << logN];

	for (long i = 0; i < L; ++i) {
		for (long n = 0; n < N; ++n) {
			mulModBarrett(p2coeff[n + (i << logN)], p, p, qVec[i], qrVec[i], qTwok[i]);
		}
	}
	p2hcoeff = new uint64_t[L << logN];

	for (long i = 0; i < L; ++i) {
		for (long n = 0; n < N; ++n) {
			mulModBarrett(p2hcoeff[n + (i << logN)], (p >> 1), p, qVec[i], qrVec[i], qTwok[i]);
		}
	}

	pccoeff = new uint64_t[L << logN]();
	for (long i = 0; i < L; ++i) {
		for (long n = 0; n < N; ++n) {
			mulModBarrett(pccoeff[n + (i << logN)], (94.2372881) * (p >> 20), (1L << 20), qVec[i], qrVec[i], qTwok[i]);
		}
	}
	negateAndEqual(pccoeff, L);


	taylorCoeffsMap.insert(pair<string, double*>(LOGARITHM, new double[11]{0,1,-0.5,1./3,-1./4,1./5,-1./6,1./7,-1./8,1./9,-1./10}));
	taylorCoeffsMap.insert(pair<string, double*>(EXPONENT, new double[11]{1,1,0.5,1./6,1./24,1./120,1./720,1./5040, 1./40320,1./362880,1./3628800}));
	taylorCoeffsMap.insert(pair<string, double*>(SIGMOID, new double[11]{1./2,1./4,0,-1./48,0,1./480,0,-17./80640,0,31./1451520,0}));

    myPrecomputation();

}

uint64_t Context::BitReverse64(uint64_t index, uint64_t bitLen)
{
    return Reverse64(index) >> (64 -bitLen);
}

uint64_t Context::Reverse64(uint64_t x)
{
    x = ((x & 0x5555555555555555) << 1) | ((x >> 1) & 0x5555555555555555);
    x = ((x & 0x3333333333333333) << 2) | ((x >> 2) & 0x3333333333333333);
    x = ((x & 0x0F0F0F0F0F0F0F0F) << 4) | ((x >> 4) & 0x0F0F0F0F0F0F0F0F);
    x = ((x & 0x00FF00FF00FF00FF) << 8) | ((x >> 8) & 0x00FF00FF00FF00FF);
    x = ((x & 0x0000FFFF0000FFFF) << 16) | ((x >> 16) & 0x0000FFFF0000FFFF);
    x = (x << 32) | (x >> 32);
    return x;
}


void Context::arrayBitReverse(complex<double>* vals, const long size) {
	for (long i = 1, j = 0; i < size; ++i) {
		long bit = size >> 1;
		for (; j >= bit; bit >>= 1) {
			j -= bit;
		}
		j += bit;
		if (i < j) {
			swap(vals[i], vals[j]);
		}
	}
}

void Context::arrayBitReverse(uint64_t* vals, const long size) {
	for (long i = 1, j = 0; i < size; ++i) {
		long bit = size >> 1;
		for (; j >= bit; bit >>= 1) {
			j -= bit;
		}
		j += bit;
		if (i < j) {
			swap(vals[i], vals[j]);
		}
	}
}

void Context::fft(complex<double>* vals, const long size) {
	arrayBitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		long lenh = len >> 1;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = j * MoverLen;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiPows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftInvLazy(complex<double>* vals, const long size) {
	arrayBitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		long MoverLen = M / len;
		long lenh = len >> 1;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < lenh; ++j) {
				long idx = (len - j) * MoverLen;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiPows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftInv(complex<double>* vals, const long size) {
	fftInvLazy(vals, size);
	for (long i = 0; i < size; ++i) {
		vals[i] /= size;
	}
}

void Context::fftSpecial(complex<double>* vals, const long size) {
	arrayBitReverse(vals, size);
	for (long len = 2; len <= size; len <<= 1) {
		for (long i = 0; i < size; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			for (long j = 0; j < lenh; ++j) {
				long idx = ((rotGroup[j] % lenq)) * M / lenq;
				complex<double> u = vals[i + j];
				complex<double> v = vals[i + j + lenh];
				v *= ksiPows[idx];
				vals[i + j] = u + v;
				vals[i + j + lenh] = u - v;
			}
		}
	}
}

void Context::fftSpecialInvLazy(complex<double>* vals, const long size) {
	for (long len = size; len >= 1; len >>= 1) {
		for (long i = 0; i < size; i += len) {
			long lenh = len >> 1;
			long lenq = len << 2;
			for (long j = 0; j < lenh; ++j) {
				long idx = (lenq - (rotGroup[j] % lenq)) * M / lenq;
				complex<double> u = vals[i + j] + vals[i + j + lenh];
				complex<double> v = vals[i + j] - vals[i + j + lenh];
				v *= ksiPows[idx];
				vals[i + j] = u;
				vals[i + j + lenh] = v;
			}
		}
	}
	arrayBitReverse(vals, size);
}

void Context::fftSpecialInv(complex<double>* vals, const long size) {
	fftSpecialInvLazy(vals, size);
	for (long i = 0; i < size; ++i) {
		vals[i] /= size;
	}
}

void Context::encode(uint64_t* a, complex<double>* v, long slots, long l) {
	complex<double>* uvals = new complex<double> [slots]();
	copy(v, v + slots, uvals);

	long gap = Nh / slots;

	fftSpecialInv(uvals, slots);

	for (long j = 0; j < slots; ++j) {
		uvals[j] *= p;
	}

	for (long i = 0; i < l; ++i) {
		uint64_t* mi = a + i * N;
		for (long j = 0, jdx = Nh, idx = 0; j < slots; ++j, jdx += gap, idx += gap) {
			long mir = uvals[j].real();
			long mii = uvals[j].imag();
			mi[idx] = mir >= 0 ? (uint64_t) mir : (uint64_t) (qVec[i] + mir);
			mi[jdx] = mii >= 0 ? (uint64_t) mii : (uint64_t) (qVec[i] + mii);
		}
		qiNTTAndEqual(mi, i);
	}
	delete[] uvals;
}

void Context::encode(uint64_t* ax, double* vals, long slots, long l) {
    //TODO implement method
}

void Context::encodeSingle(uint64_t* ax, complex<double>& val, long l) {
	long vr = val.real() * p;
	long vi = val.imag() * p;

	for (long i = 0; i < l; ++i) {
		uint64_t* ai = ax + i * N;
		ai[0] = vr >= 0 ? vr : qVec[i] + vr;
		ai[Nh] = vi >= 0 ? vi : qVec[i] + vi;
		qiNTTAndEqual(ai, i);
	}
}

void Context::encodeSingle(uint64_t* ax, double val, long l) {

}

void Context::decode(uint64_t* a, complex<double>* v, long slots, long l) {
	uint64_t* tmp = new uint64_t[N]();
	copy(a, a + N, tmp);
	long gap = Nh / slots;
	qiINTTAndEqual(tmp, 0);

	uint64_t pr = qVec[0];
	uint64_t pr_2 = qVec[0] / 2;

	for (long j = 0, jdx = Nh, idx = 0; j < slots; ++j, jdx += gap, idx += gap) {
		double mir = tmp[idx] <= pr_2 ? ((double) (tmp[idx]) / p) : (((double) (tmp[idx]) - (double) (pr)) / p);
		double mii = tmp[jdx] <= pr_2 ? ((double) (tmp[jdx]) / p) : (((double) (tmp[jdx]) - (double) (pr)) / p);
		v[j].real(mir);
		v[j].imag(mii);
	}
	fftSpecial(v, slots);
}

void Context::decodeSingle(uint64_t* ax, complex<double>& val, long l) {
	uint64_t* tmp = new uint64_t[N]();
	copy(ax, ax + N, tmp);
	qiINTTAndEqual(tmp, 0);

	uint64_t pr = qVec[0];
	uint64_t pr_2 = qVec[0] / 2;

	double vr = tmp[0] <= pr_2 ? ((double) tmp[0]) / p : (((double) tmp[0]) - ((double) pr)) / (double) p;
	double vi = tmp[Nh] <= pr_2 ? ((double) tmp[Nh]) / p : (((double) tmp[Nh]) - ((double) pr)) / (double) p;

	val.real(vr);
	val.imag(vi);
}

void Context::decodeSingle(uint64_t* ax, double val, long l) {
	//TODO implement method

}

void Context::qiNTT(uint64_t* res, uint64_t* a, long index) {
	copy(a, a + N, res);
	qiNTTAndEqual(res, index);
}

void Context::piNTT(uint64_t* res, uint64_t* a, long index) {
	copy(a, a + N, res);
	piNTTAndEqual(res, index);
}

void Context::NTT(uint64_t* res, uint64_t* a, long l, long k) {
	for (long index = 0; index < l; ++index) {
		uint64_t* ai = a + (index << logN);
		uint64_t* resi = a + (index << logN);
		qiNTT(resi, ai, index);
	}

	for (long index = l; index < l + k; ++index) {
		uint64_t* ai = a + (index << logN);
		uint64_t* resi = a + (index << logN);
		piNTT(resi, ai, index - l);
	}
}

void Context::qiNTTAndEqual(uint64_t* a, long index) {
	long t = N;
	long logt1 = logN + 1;
	uint64_t q = qVec[index];
//	uint64_t qd = qdVec[index];
	uint64_t qInv = qInvVec[index];
	for (long m = 1; m < N; m <<= 1) {
		t >>= 1;
		logt1 -= 1;
		for (long i = 0; i < m; i++) {
			long j1 = i << logt1;
			long j2 = j1 + t - 1;
			uint64_t W = qRootScalePows[index][m + i];
//			uint64_t W = qRootScalePowsOverq[index][m + i];
//			uint64_t w = qRootPows[index][m + i];
			for (long j = j1; j <= j2; j++) {

				uint64_t T = a[j + t];
				unsigned __int128 U = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(U);
				uint64_t U1 = static_cast<uint64_t>(U >> 64);
				uint64_t Q = U0 * qInv;
				unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * q;
				uint64_t H = static_cast<uint64_t>(Hx >> 64);
				uint64_t V = U1 < H ? U1 + q - H : U1 - H;
				a[j + t] = a[j] < V ? a[j] + q - V: a[j] - V;
				a[j] += V;
				if(a[j] > q) a[j] -= q;

//				if(a[j] >= qd) a[j] -= qd;
//				uint64_t T = a[j + t];
//				unsigned __int128 U = static_cast<unsigned __int128>(T) * W;
//				uint64_t Q = static_cast<uint64_t>(U >> 64);
//				T *= w;
//				uint64_t T1 = Q * q;
//				T -= T1;
//				a[j + t] = a[j] + qd - T;
//				a[j] += T;
			}
		}
	}
//	for(long i = 0; i < N; i++) {
//		if(a[i] >= qd) a[i] -= qd;
//		if(a[i] >= q) a[i] -= q;
//	}
}

void Context::piNTTAndEqual(uint64_t* a, long index) {
	long t = N;
	long logt1 = logN + 1;
	uint64_t pi = pVec[index];
//	uint64_t pd = pdVec[index];
	uint64_t pInv = pInvVec[index];
	for (long m = 1; m < N; m <<= 1) {
		t >>= 1;
		logt1 -= 1;
		for (long i = 0; i < m; i++) {
			long j1 = i << logt1;
			long j2 = j1 + t - 1;
			uint64_t W = pRootScalePows[index][m + i];
//			uint64_t W = pRootScalePowsOverp[index][m + i];
//			uint64_t w = pRootPows[index][m + i];
			for (long j = j1; j <= j2; j++) {

				uint64_t T = a[j + t];
				unsigned __int128 U = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(U);
				uint64_t U1 = static_cast<uint64_t>(U >> 64);
				uint64_t Q = U0 * pInv;
				unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * pi;
				uint64_t H = static_cast<uint64_t>(Hx >> 64);
				uint64_t V = U1 < H ? U1 + pi - H : U1 - H;
				a[j + t] = a[j] < V ? a[j] + pi - V: a[j] - V;
				a[j] += V;
				if(a[j] > pi) a[j] -= pi;

//				if(a[j] >= pd) a[j] -= pd;
//				uint64_t T = a[j + t];
//				unsigned __int128 U = static_cast<unsigned __int128>(T) * W;
//				uint64_t Q = static_cast<uint64_t>(U >> 64);
//				T *= w;
//				uint64_t T1 = Q * pi;
//				T -= T1;
//				a[j + t] = a[j] + pd - T;
//				a[j] += T;
			}
		}
	}
//	for(long i = 0; i < N; i++) {
//		if(a[i] >= pd) a[i] -= pd;
//		if(a[i] >= p) a[i] -= pi;
//	}
}

void Context::NTTAndEqual(uint64_t* a, long l, long k) {
	for (long index = 0; index < l; ++index) {
		uint64_t* ai = a + (index << logN);
		qiNTTAndEqual(ai, index);
	}

	for (long index = l; index < l + k; ++index) {
		uint64_t* ai = a + (index << logN);
		piNTTAndEqual(ai, index - l);
	}
}

void Context::qiINTT(uint64_t* res, uint64_t* a, long index) {
	copy(a, a + N, res);
	qiINTTAndEqual(res, index);
}

void Context::qiINTTAndEqual(uint64_t* a, long index) {
	uint64_t q = qVec[index];
	uint64_t qd = qdVec[index];
	uint64_t qInv = qInvVec[index];
	long t = 1;
	for (long m = N; m > 1; m >>= 1) {
		long j1 = 0;
		long h = m >> 1;
		for (long i = 0; i < h; i++) {
			long j2 = j1 + t - 1;
			uint64_t W = qRootScalePowsInv[index][h + i];
			for (long j = j1; j <= j2; j++) {
				uint64_t T = a[j] + qd;
				T -= a[j + t];
				a[j] += a[j + t];
				if(a[j] >= qd) a[j] -= qd;
				unsigned __int128 UU = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(UU);
				uint64_t U1 = static_cast<uint64_t>(UU >> 64);
				uint64_t Q = U0 * qInv;
				unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * q;
				uint64_t H = static_cast<uint64_t>(Hx >> 64);
				a[j + t] = U1 + q;
				a[j + t] -= H;
			}
			j1 += (t << 1);
		}
		t <<= 1;
	}

	uint64_t NScale = NScaleInvModq[index];
	for (long i = 0; i < N; i++) {
		uint64_t T = (a[i] < q) ? a[i] : a[i] - q;
		unsigned __int128 U = static_cast<unsigned __int128>(T) * NScale;
		uint64_t U0 = static_cast<uint64_t>(U);
		uint64_t U1 = static_cast<uint64_t>(U >> 64);
		uint64_t Q = U0 * qInv;
		unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * q;
		uint64_t H = static_cast<uint64_t>(Hx >> 64);
		a[i] = (U1 < H) ? U1 + q - H : U1 - H;
	}
}

void Context::piINTTAndEqual(uint64_t* a, long index) {
	uint64_t pi = pVec[index];
	uint64_t pd = pdVec[index];
	uint64_t pInv = pInvVec[index];
	long t = 1;
	for (long m = N; m > 1; m >>= 1) {
		long j1 = 0;
		long h = m >> 1;
		for (long i = 0; i < h; i++) {
			long j2 = j1 + t - 1;
			uint64_t W = pRootScalePowsInv[index][h + i];
			for (long j = j1; j <= j2; j++) {
				uint64_t T = a[j] + pd;
				T -= a[j + t];
				a[j] += a[j + t];
				if(a[j] >= pd) a[j] -= pd;
				unsigned __int128 UU = static_cast<unsigned __int128>(T) * W;
				uint64_t U0 = static_cast<uint64_t>(UU);
				uint64_t U1 = static_cast<uint64_t>(UU >> 64);
				uint64_t Q = U0 * pInv;
				unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * pi;
				uint64_t H = static_cast<uint64_t>(Hx >> 64);
				a[j + t] = U1 + pi;
				a[j + t] -= H;
			}
			j1 += (t << 1);
		}
		t <<= 1;
	}

	uint64_t NScale = NScaleInvModp[index];
	for (long i = 0; i < N; i++) {
		uint64_t T = (a[i] < pi) ? a[i] : a[i] - pi;
		unsigned __int128 U = static_cast<unsigned __int128>(T) * NScale;
		uint64_t U0 = static_cast<uint64_t>(U);
		uint64_t U1 = static_cast<uint64_t>(U >> 64);
		uint64_t Q = U0 * pInv;
		unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * pi;
		uint64_t H = static_cast<uint64_t>(Hx >> 64);
		a[i] = (U1 < H) ? U1 + pi - H : U1 - H;
	}
}

void Context::INTTAndEqual(uint64_t* a, long l, long k) {
	for (long index = 0; index < l; ++index) {
		uint64_t* ai = a + (index << logN);
		qiINTTAndEqual(ai, index);
	}

	for (long index = l; index < l + k; ++index) {
		uint64_t* ai = a + (index << logN);
		piINTTAndEqual(ai, index - l);
	}
}

void Context::qiNegate(uint64_t* res, uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		res[i] = (a[i] == 0) ? 0 : qVec[index] - a[i];
	}
}

void Context::piNegate(uint64_t* res, uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		res[i] = (a[i] == 0) ? 0 : pVec[index] - a[i];
	}
}

void Context::negate(uint64_t* res, uint64_t* a, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiNegate(resi, ai, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		piNegate(resi, ai, i - l);
	}
}

void Context::qiNegateAndEqual(uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		a[i] = (a[i] == 0) ? 0 : qVec[index] - a[i];
	}
}

void Context::piNegateAndEqual(uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		a[i] = (a[i] == 0) ? 0 : pVec[index] - a[i];
	}
}

void Context::negateAndEqual(uint64_t* a, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		qiNegateAndEqual(ai, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		piNegateAndEqual(ai, i - l);
	}
}

void Context::qiAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(res[i], a[i], c, qVec[index]);
	}
}

void Context::piAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(res[i], a[i], c, pVec[index]);
	}
}

void Context::addConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiAddConst(resi, ai, c, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		piAddConst(resi, ai, c, i - l);
	}
}

void Context::qiAddConstAndEqual(uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(a[i], a[i], c, qVec[index]);
	}
}

void Context::piAddConstAndEqual(uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(a[i], a[i], c, pVec[index]);
	}
}

void Context::addConstAndEqual(uint64_t* a, uint64_t c, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		qiAddConstAndEqual(ai, c, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		piAddConstAndEqual(ai, c, i - l);
	}
}

void Context::qiSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(res[i], a[i], c, qVec[index]);
	}
}

void Context::piSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(res[i], a[i], c, pVec[index]);
	}
}

void Context::subConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiSubConst(resi, ai, c, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		piSubConst(resi, ai, c, i - l);
	}
}

void Context::qiSubConstAndEqual(uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(a[i], a[i], c, qVec[index]);
	}
}

void Context::piSubConstAndEqual(uint64_t* a, uint64_t c, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(a[i], a[i], c, pVec[index]);
	}
}

void Context::subConstAndEqual(uint64_t* a, uint64_t c, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		qiSubConstAndEqual(ai, c, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		piSubConstAndEqual(ai, c, i - l);
	}
}

void Context::qiAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(res[i], a[i], b[i], qVec[index]);
	}
}

void Context::piAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(res[i], a[i], b[i], pVec[index]);
	}
}

void Context::add(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiAdd(resi, ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		uint64_t* resi = res + (i << logN);
		piAdd(resi, ai, bi, i - l);
	}
}

void Context::qiAddAndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(a[i], a[i], b[i], qVec[index]);
	}
}

void Context::piAddAndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		addMod(a[i], a[i], b[i], pVec[index]);
	}
}

void Context::addAndEqual(uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		qiAddAndEqual(ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		piAddAndEqual(ai, bi, i - l);
	}
}

void Context::qiSub(uint64_t* res, uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(res[i], a[i], b[i], qVec[index]);
	}
}

void Context::piSub(uint64_t* res, uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(res[i], a[i], b[i], pVec[index]);
	}
}

void Context::sub(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiSub(resi, ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		uint64_t* resi = res + (i << logN);
		piSub(resi, ai, bi, i - l);
	}
}

void Context::qiSubAndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(a[i], a[i], b[i], qVec[index]);
	}
}

void Context::piSubAndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(a[i], a[i], b[i], pVec[index]);
	}
}

void Context::subAndEqual(uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		qiSubAndEqual(ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		piSubAndEqual(ai, bi, i - l);
	}
}

void Context::qiSub2AndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(b[i], a[i], b[i], qVec[index]);
	}
}

void Context::piSub2AndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		subMod(b[i], a[i], b[i], pVec[index]);
	}
}

void Context::sub2AndEqual(uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		qiSub2AndEqual(ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		piSub2AndEqual(ai, bi, i - l);
	}
}

void Context::qiMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], a[i], cnst, qVec[index], qrVec[index], qTwok[index]);
	}
}

void Context::piMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], a[i], cnst, pVec[index], prVec[index], pTwok[index]);
	}
}

void Context::mulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiMulConst(resi, ai, cnst, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		piMulConst(resi, ai, cnst, i - l);
	}
}

void Context::qiMulConstAndEqual(uint64_t* res, uint64_t cnst, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], res[i], cnst, qVec[index], qrVec[index], qTwok[index]);
	}
}

void Context::piMulConstAndEqual(uint64_t* res, uint64_t cnst, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], res[i], cnst, pVec[index], prVec[index], pTwok[index]);
	}
}

void Context::mulConstAndEqual(uint64_t* a, uint64_t cnst, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		qiMulConstAndEqual(ai, cnst, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		piMulConstAndEqual(ai, cnst, i - l);
	}
}

void Context::qiMul(uint64_t* res, uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], a[i], b[i], qVec[index], qrVec[index], qTwok[index]);
	}
}

void Context::piMul(uint64_t* res, uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], a[i], b[i], pVec[index], prVec[index], pTwok[index]);
	}
}

void Context::mul(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiMul(resi, ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		uint64_t* resi = res + (i << logN);
		piMul(resi, ai, bi, i - l);
	}
}

void Context::mulKey(uint64_t* res, uint64_t* a, uint64_t* key, long l) {
	for (long i = 0; i < l; ++i) {
        //compute the index of polys
		uint64_t* ai = a + (i << logN);
		uint64_t* keyi = key + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiMul(resi, ai, keyi, i);
	}

	for (long i = l; i < l + K; ++i) {
        //compute the index of polys
		uint64_t* ai = a + (i << logN);
		uint64_t* keyi = key + ((i - l + L) << logN);
		uint64_t* resi = res + (i << logN);
		piMul(resi, ai, keyi, i - l);
	}
}

void Context::qiMulAndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(a[i], a[i], b[i], qVec[index], qrVec[index], qTwok[index]);
	}
}

void Context::piMulAndEqual(uint64_t* a, uint64_t* b, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(a[i], a[i], b[i], pVec[index], prVec[index], pTwok[index]);
	}
}

void Context::mulAndEqual(uint64_t* a, uint64_t* b, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		qiMulAndEqual(ai, bi, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* bi = b + (i << logN);
		piMulAndEqual(ai, bi, i - l);
	}
}

void Context::qiSquare(uint64_t* res, uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], a[i], a[i], qVec[index], qrVec[index], qTwok[index]);
	}
}

void Context::piSquare(uint64_t* res, uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(res[i], a[i], a[i], pVec[index], prVec[index], pTwok[index]);
	}
}

void Context::square(uint64_t* res, uint64_t* a, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		qiSquare(resi, ai, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		uint64_t* resi = res + (i << logN);
		piSquare(resi, ai, i - l);
	}
}

void Context::qiSquareAndEqual(uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(a[i], a[i], a[i], qVec[index], qrVec[index], qTwok[index]);
	}
}

void Context::piSquareAndEqual(uint64_t* a, long index) {
	for (long i = 0; i < N; ++i) {
		mulModBarrett(a[i], a[i], a[i], pVec[index], prVec[index], pTwok[index]);
	}
}

void Context::squareAndEqual(uint64_t* a, long l, long k) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		qiSquare(ai, ai, i);
	}

	for (long i = l; i < l + k; ++i) {
		uint64_t* ai = a + (i << logN);
		piSquare(ai, ai, i - l);
	}
}

void Context::evalAndEqual(uint64_t* a, long l) {
	INTTAndEqual(a, l);
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		for (long n = 0; n < N; ++n) {
			mulModBarrett(ai[n], ai[n], PModq[i], qVec[i], qrVec[i], qTwok[i]);
		}
	}
	NTTAndEqual(a, l);
}

void Context::raise(uint64_t* res, uint64_t* a, long l) {
	//TODO implement method

}

void Context::raiseAndEqual(uint64_t*& a, long l) {
	uint64_t* ra = new uint64_t[(l + K) << logN]();
	copy(a, a + (l << logN), ra);

	INTTAndEqual(a, l);

	uint64_t* tmp3 = new uint64_t[l << logN];
	for(long i = 0; i < l; ++i) {
		uint64_t* tmp3i = tmp3 + (i << logN);
		uint64_t* ai = a + (i << logN);
		for(long n = 0; n < N; ++n) {
			mulModBarrett(tmp3i[n], ai[n], qHatInvModq[l - 1][i], qVec[i], qrVec[i], qTwok[i]);
		}
	}
	for (long k = 0; k < K; ++k) {
		uint64_t* rak = ra + ((l + k) << logN);
		for (long n = 0; n < N; ++n) {
			uint64_t tt = tmp3[n];
			unsigned __int128 sum = static_cast<unsigned __int128>(tt) * qHatModp[l - 1][0][k];
			for (long i = 1; i < l; ++i) {
				tt = tmp3[n + (i << logN)];
				sum += static_cast<unsigned __int128>(tt) * qHatModp[l - 1][i][k];
			}
			modBarrett(rak[n], sum, pVec[k], prVec[k], pTwok[k]);
		}
	}
	NTTAndEqual(ra + (l << logN), 0, K);

	delete[] a;
	a = ra;

}

void Context::back(uint64_t* res, uint64_t* a, long l) {
	uint64_t* tmp = new uint64_t[(l + K) << logN];
	copy(a, a + ((l + K) << logN), tmp);
	INTTAndEqual(tmp, l, K);
	uint64_t* tmp3 = new uint64_t[K << logN];
	for(long k = 0; k < K; k++) {
		uint64_t* tmpk = tmp + ((k + l) << logN);
		uint64_t* tmp3k = tmp3 + (k << logN);
		for(long n = 0; n < N; ++n) {
			mulModBarrett(tmp3k[n], tmpk[n], pHatInvModp[k], pVec[k], prVec[k], pTwok[k]);
		}
	}
	for (long i = 0; i < l; ++i) {
		uint64_t* resi = res + (i << logN);
		uint64_t* tmpi = tmp + (i << logN);

		for (long n = 0; n < N; ++n) {
			uint64_t tt = tmp3[n];
			unsigned __int128 sum = static_cast<unsigned __int128>(tt) * pHatModq[0][i];
			for (long k = 1; k < K; ++k) {
				tt = tmp3[n + (k << logN)];
				sum += static_cast<unsigned __int128>(tt) * pHatModq[k][i];
			}
			modBarrett(resi[n], sum, qVec[i], qrVec[i], qTwok[i]);
			subMod(resi[n], tmpi[n], resi[n], qVec[i]);
			mulModBarrett(resi[n], resi[n], PInvModq[i], qVec[i], qrVec[i], qTwok[i]);
		}
	}

	delete[] tmp;

	NTTAndEqual(res, l);
}

void Context::backAndEqual(uint64_t*& a, long l) {

	INTTAndEqual(a, l, K);

	uint64_t* ra = new uint64_t[l << logN]();
	uint64_t* tmp3 = new uint64_t[K << logN];

	for(long k = 0; k < K; k++) {
		uint64_t* tmp3k = tmp3 + (k << logN);
		uint64_t* ak = a + ((k + l) << logN);
		for(long n = 0; n < N; ++n) {
			mulModBarrett(tmp3k[n], ak[n], pHatInvModp[k], pVec[k], prVec[k], pTwok[k]);
		}
	}

	for (long i = 0; i < l; ++i) {
		uint64_t* rai = ra + (i << logN);
		uint64_t* ai = a + (i << logN);
		for (long n = 0; n < N; ++n) {
			uint64_t tt = tmp3[n];
			unsigned __int128 sum = static_cast<unsigned __int128>(tt) * pHatModq[0][i];
			for (long k = 1; k < K; ++k) {
				tt = tmp3[n + (k << logN)];
				sum += static_cast<unsigned __int128>(tt) * pHatModq[k][i];
			}
			modBarrett(rai[n], sum, qVec[i], qrVec[i], qTwok[i]);
			subMod(rai[n], ai[n], rai[n], qVec[i]);
			mulModBarrett(rai[n], rai[n], PInvModq[i], qVec[i], qrVec[i], qTwok[i]);
		}
	}
	NTTAndEqual(ra, l);
	delete[] a;
	a = ra;
}

void Context::reScale(uint64_t* res, uint64_t* a, long l) {
	//TODO implement method
}

void Context::reScaleAndEqual(uint64_t*& a, long l) {
    uint64_t *ra = new uint64_t[(l - 1) << logN]();
    uint64_t *al = a + ((l - 1) << logN);
    qiINTTAndEqual(al, l - 1);
    for (long i = 0; i < l - 1; ++i) {//openfhe: l-1=extra.m_vectors.size()
        uint64_t *rai = ra + (i << logN);
        uint64_t *ai = a + (i << logN);

        for (long n = 0; n < N; ++n) {
            modBarrett(rai[n], al[n], qVec[i], qrVec[i], qTwok[i]);//openfhe: temp.SwitchModulus()
        }
        qiNTTAndEqual(rai, i);


        for (long n = 0; n < N; ++n) {
            subMod(rai[n], ai[n], rai[n], qVec[i]);
            mulModBarrett(rai[n], rai[n], qInvModq[l - 1][i], qVec[i], qrVec[i], qTwok[i]);
        }
    }
	delete[] a;
	a = ra;
}

void Context::modDownAndEqual(uint64_t*& a, long l, long dl) {
	uint64_t* ra = new uint64_t[(l - dl) << logN]();
	copy(a, a + ((l - dl) << logN), ra);
	delete[] a;
	a = ra;
}

uint64_t* Context::modDown(uint64_t* a, long l, long dl) {
	uint64_t* ra = new uint64_t[(l - dl) << logN]();
	copy(a, a + ((l - dl) << logN), ra);
	return ra;
}

void Context::leftRot(uint64_t* res, uint64_t* a, long l, long rotSlots) {
//	long idx = rotSlots % Nh;
//	for (long n = 0; n < N; ++n) {
//		uint32_t reversed = bitReverse(static_cast<uint32_t>(n)) >> (32 - logN);
//		uint64_t index_raw = rotGroup[idx] * (2 * reversed + 1);
//		index_raw &= (M - 1);
//		long index = bitReverse((static_cast<uint32_t>(index_raw) - 1) >> 1) >> (32 - logN);
//		for (long i = 0; i < curr_limbs; ++i) {
//			res[n + (i << logN)] = a[index + (i << logN)];
//		}
//	}

	uint64_t* tmp = new uint64_t[l << logN]();
	copy(a, a + (l << logN), tmp);
	INTTAndEqual(tmp, l);
	long pow = rotGroup[rotSlots];
	for (long i = 0; i < l; ++i) {
		uint64_t* resi = res + (i << logN);
		uint64_t* tmpi = tmp + (i << logN);
		for (long n = 0; n < N; ++n) {
			long npow = n * pow;
			long shift = npow % M;
			if(shift < N) {
				resi[shift] = tmpi[n];
			} else {
				resi[shift - N] = qVec[i] - tmpi[n];
			}
		}
	}
	NTTAndEqual(res, l);
}

void Context::leftRotAndEqual(uint64_t* a, long l, long rotSlots) {
	uint64_t* tmp = new uint64_t[l << logN];
	copy(a, a + (l << logN), tmp);
	long idx = rotSlots % Nh;
	for (long n = 0; n < N; ++n) {
		uint32_t reversed = bitReverse(static_cast<uint32_t>(n)) >> (32 - logN);
		uint64_t index_raw = rotGroup[idx] * (2 * reversed + 1);
		index_raw &= M - 1;
		uint32_t index = bitReverse((static_cast<uint32_t>(index_raw) - 1) >> 1) >> (32 - logN);
		for (long i = 0; i < l; ++i) {
			a[n + (i << logN)] = tmp[index + (i << logN)];
		}
	}
//	INTTAndEqual(tmp, curr_limbs);
//	long pow = rotGroup[rotSlots];
//	for (long i = 0; i < curr_limbs; ++i) {
//		uint64_t* ai = a + (i << logN);
//		uint64_t* tmpi = tmp + (i << logN);
//		for (long n = 0; n < N; ++n) {
//			long jpow = n * pow;
//			long shift = jpow % M;
//			if(shift < N) {
//				ai[shift] = tmpi[n];
//			} else {
//				ai[shift - N] = qVec[i] - tmpi[n];
//			}
//		}
//	}
//	NTTAndEqual(a, curr_limbs);
}

void Context::conjugate(uint64_t* res, uint64_t* a, long l) {
	for (long i = 0; i < l; ++i) {
		uint64_t* resi = res + (i << logN);
		uint64_t* ai = a + (i << logN);
		for (int n = 0; n < N; ++n) {
			resi[n] = ai[N - 1 - n];
		}
	}
}

void Context::conjugateAndEqual(uint64_t* a, long l) {
	for (long i = 0; i < l; ++i) {
		uint64_t* ai = a + (i << logN);
		for (int n = 0; n < N; ++n) {
			swap(ai[n], ai[N - 1 - n]);
		}
	}
}

void Context::mulByMonomial(uint64_t* res, uint64_t* a, long l, long monomialDeg) {
    long shift=monomialDeg% M ; // M is the Cyclotomic order
    if (shift==0){
        //directly copy the input
        for (long i = 0; i < l; ++i) {
            uint64_t* resi = res + (i << logN);
            uint64_t* ai = a + (i << logN);
            for (long n = 0; n < N; ++n)
                resi[n]=ai[n];
        }
        return;
    }

    if (shift<N){
        for (long i = 0; i < l; ++i) {
            uint64_t* resi = res + (i << logN);
            uint64_t* ai = a + (i << logN);
            for (long n = 0; n < shift; ++n)
                resi[n]=qVec[i]-ai[N - shift + n];
            for (long n = shift; n < N; ++n)
                resi[n]=ai[n-shift];
        }
    }
    else{
        //merge the negate operation and the assignment
        shift%=N;
        for (long i = 0; i < l; ++i) {
            uint64_t* resi = res + (i << logN);
            uint64_t* ai = a + (i << logN);
            for (long n = 0; n < shift; ++n)
                resi[n]=ai[N - shift + n];
            for (long n = shift; n < N; ++n)
                resi[n]=qVec[i]-ai[n-shift];
        }
    }

}

void Context::mulByMonomialAndEqual(uint64_t* a, long l, long monomialDeg) {

    long shift=monomialDeg% M ; // M is the Cyclotomic order
    if (shift==0){
       return;
    }

    //todo: 可以借鉴HEAX ntt 不用辅助数组的实现方法省略此处的tmp数组？
    uint64_t * tmp=new uint64_t[l*N];
    if (shift<N){
        for (long i = 0; i < l; ++i) {
            uint64_t* tmpi = tmp + (i << logN);
            uint64_t* ai = a + (i << logN);
            for (long n = 0; n < N; ++n)
                tmpi[n]=ai[n];
        }
    }
    else{
        //negate
        for (long i = 0; i < l; ++i) {
            uint64_t* tmpi = tmp + (i << logN);
            uint64_t* ai = a + (i << logN);
            for (long n = 0; n < N; ++n)
                tmpi[n]=qVec[i]-ai[n];
        }
    }

    shift %=N;
    for (long i = 0; i < l; ++i) {
        uint64_t* tmpi = tmp + (i << logN);
        uint64_t* ai = a + (i << logN);
        for (long n = 0; n < shift; ++n)
            ai[n]=qVec[i]-tmpi[N - shift + n];
        for (long n  = shift; n < N; n++)
            ai[n]=tmpi[n-shift];
    }

    delete[] tmp;


}

void Context::sampleGauss(uint64_t* res, long l, long k) {
	static long const bignum = 0xfffffff;
	for (long i = 0; i < N; i += 2) {
		double r1 = (1 + (uint64_t)rand() % bignum) / ((double) bignum + 1);
		double r2 = (1 + (uint64_t)rand() % bignum) / ((double) bignum + 1);
		double theta = 2 * M_PI * r1;
		double rr = sqrt(-2.0 * log(r2)) * sigma;

		long g1 = floor(rr * cos(theta) + 0.5);
		long g2 = floor(rr * sin(theta) + 0.5);

		for (long j = 0; j < l; ++j) {
			uint64_t* resj = res + (j << logN);
			resj[i] = g1 >= 0 ? g1 : qVec[j] + g1;
			resj[i + 1] = g2 >= 0 ? g2 : qVec[j] + g2;
		}
		for (long j = 0; j < k; ++j) {
			uint64_t* resj = res + ((j + l) << logN);
			resj[i] = g1 >= 0 ? g1 : pVec[j] + g1;
			resj[i + 1] = g2 >= 0 ? g2 : pVec[j] + g2;
		}
	}
}

void Context::sampleZO(uint64_t* res, long s, long l, long k) {
	for (long i = 0; i < N; ++i) {
		long zo = (rand() % 2) == 0 ? 0 : (rand() % 2) ? 1 : -1;
		for (long j = 0; j < l; ++j) {
			uint64_t* resj = res + (j << logN);
			resj[i] = zo >= 0 ? zo : qVec[j] + zo;
		}
		for (long j = 0; j < k; ++j) {
			uint64_t* resj = res + ((j + l) << logN);
			resj[i] = zo >= 0 ? zo : pVec[j] + zo;
		}
	}
}

void Context::sampleUniform(uint64_t* res, long l, long k) {
	for (long j = 0; j < l; ++j) {
		uint64_t* resj = res + (j << logN);
		for (long n = 0; n < N; ++n) {
			resj[n] = floor(((double) rand() / (RAND_MAX)) * qVec[j]);
		}
	}
	for (long j = 0; j < k; ++j) {
		uint64_t* resj = res + ((j + l) << logN);
		for (long n = 0; n < N; ++n) {
			resj[n] = floor(((double) rand() / (RAND_MAX)) * pVec[j]);
		}
	}
}

void Context::sampleHWT(uint64_t* res, long l, long k) {
	long idx = 0;
	while (idx < h) {
		long i = ((double) rand() / (RAND_MAX)) * N;
		if (res[i] == 0) {
			long hwt = (rand() % 2) ? 1 : -1;
			for (long j = 0; j < l; ++j) {
				uint64_t* resj = res + (j << logN);
				resj[i] = hwt >= 0 ? hwt : qVec[j] + hwt;
			}
			for (long j = 0; j < k; ++j) {
				uint64_t* resj = res + ((j + l) << logN);
				resj[i] = hwt >= 0 ? hwt : pVec[j] + hwt;
			}
			idx++;
		}
	}
}

// By Czh, tested and right
uint64_t Context::ScaleUpExact(double value, long n, uint64_t q)
{
    bool isNegative = false;
    if (value < 0) {
        isNegative = true;
        value = -n * value;
    } else {
        value = n * value;
    }

    value += 0.5;
    uint64_t res = static_cast<uint64_t>(value) % q;

    if (isNegative) {
        res = q - res;
    }

    return res;
}

// By Czh, tested and right
uint64_t Mul64(uint64_t x,uint64_t y, uint64_t &lo)
{
    uint64_t mask32 = 4294967295; // 1 << 32 - 1
    uint64_t x0 = x & mask32;
    uint64_t x1 = x >> 32;
    uint64_t y0 = y & mask32;
    uint64_t y1 = y >> 32;
    uint64_t w0 = x0 * y0;
    uint64_t t = x1*y0 + (w0 >> 32);
    uint64_t w1 = t & mask32;
    uint64_t w2 = t >>32;
    w1 += x0 * y1;
    uint64_t hi = x1*y1 + w2 + (w1>>32);
    lo = x * y;
    return hi;
}

void Add64(uint64_t x, uint64_t y, uint64_t carry, uint64_t& sum, uint64_t& carryout)
{
    sum = x + y + carry;
    carryout = ((x & y) | ((x | y) &~ sum)) >> 63;
}

uint64_t BRed(uint64_t x, uint64_t y, uint64_t q, uint64_t u[])
{
    uint64_t lhi, mhi, mlo, s0, s1, carry;
    uint64_t ahi, alo, useless, useless2;
    ahi = Mul64(x, y, alo);
    lhi = Mul64(alo, u[1], useless);
    mhi = Mul64(alo, u[0], mlo);
    Add64(mlo, lhi, 0, s0, carry);
    s1 = mhi + carry;
    mhi = Mul64(ahi, u[1], mlo);
    Add64(mlo, s0, 0, useless2, carry);
    lhi = mhi + carry;
    s0 = ahi*u[0] + s1 + lhi;
    auto res = alo - s0*q;
    if (res >= q)
    {
        res -= q;
    }
    return res;
}

uint64_t Context::BRedAdd(uint64_t x, uint64_t q, uint64_t u[])
{
    uint64_t useless;
    auto s0 = Mul64(x, u[0], useless);
    auto r = x - s0*q;
    if (r >= q)
    {
        r -= q;
    }
    return r;
}


uint64_t Context::ModExp(uint64_t x, uint64_t e, uint64_t p) {
    uint64_t params[] = {0, 0};
    Context::ComputeBRedParameters(p, params[0], params[1]);
    uint64_t result = 1;
    for (long i = e; i > 0; i >>= 1)
    {
        if (i&1 == 1)
        {
            result = BRed(result, x, p, params);
        }
        x = BRed(x, x, p, params);
    }
    return result;
}


// By Czh, tested and right
uint64_t Context::MRed(uint64_t x, uint64_t y, uint64_t q, uint64_t qInv)
{
    uint64_t ahi, alo;
    ahi = Mul64(x, y, alo);
    auto R = alo * qInv;
    uint64_t H, Useless;
    H = Mul64(R, q, Useless);
    uint64_t r = ahi - H + q;

    if (r >= q)
    {
        r -= q;
    }
    return r;
}

void Context::MRed2(uint64_t *x, uint64_t y, uint64_t q, uint64_t qInv, uint64_t &r)
{
    uint64_t ahi, alo;
    ahi = Mul64(*x, y, alo);
    auto R = alo * qInv;
    uint64_t H, Useless;
    H = Mul64(R, q, Useless);
    r = ahi - H + q;

    if (r >= q)
    {
        r -= q;
    }
}

// By Czh, tested and right
uint64_t Context::MForm(uint64_t a, uint64_t q, const uint64_t u[])
{
    uint64_t mhi, mlo;
    mhi = Mul64(a, u[1], mlo);
    uint64_t r = -(a*u[0] + mhi) * q;
    if (r >= q)
    {
        r -= q;
    }
    return r;
}

uint64_t Context::CRed(uint64_t a, uint64_t q)
{
    if (a >= q)
    {
        return a - q;
    }
    return  a;
}

using namespace boost::multiprecision;

// convert a bred para into uint64{mhi, mlo}, this arr is meg para
void Context::ComputeBRedParameters(uint64_t q, uint64_t& mhi, uint64_t& mlo) {
    cpp_int bigR = cpp_int(1) << 128;
    bigR /= q;

    mhi = static_cast<uint64_t>(bigR >> 64);
    mlo = static_cast<uint64_t>(bigR);
}

int Context::Len64(uint64_t x) {
    if (x == 0)
        return 0;

    int len = 0;
    while (x != 0) {
        x >>= 1;
        len++;
    }

    return len;
}

void Context::fulldecode(uint64_t* a, complex<double>* v, long slots, long l) {
    auto tmp = new uint64_t[l << logN]();
    auto bigintCoeffs = new cpp_int[N]();
    copy(a, a + (l << logN), tmp);

    INTTAndEqual(tmp, l, 0);
    PolyToBigint(tmp, bigintCoeffs, l);

    auto Q = bigintChain[l-1];
    auto qHalf = Q >> 1;
    long gap = Nh / slots;

    int sign;

    for (uint64_t i = 0, idx = 0; i < slots; i = i + 1, idx = idx + gap)
    {
        bigintCoeffs[idx] %= Q;
        sign = cmp(bigintCoeffs[idx], qHalf);
        if (sign == 1 | sign == 0)
        {
            bigintCoeffs[idx] -= Q;
        }

        bigintCoeffs[idx+Nh] %= Q;
        sign = cmp(bigintCoeffs[idx+Nh], qHalf);
        if (sign == 1 | sign == 0)
        {
            bigintCoeffs[idx+Nh] -= Q;
        }
        v[i] = complex(scaleDown(bigintCoeffs[idx], double(p)), scaleDown(bigintCoeffs[idx+Nh], double(p)));
    }

    fftSpecial(v, slots);
}

double Context::scaleDown(cpp_int& coeff, double n)
{
    return coeff.convert_to<double>()/double(n);
//    cpp_dec_float_100 x(coeff.str());
//    x /= n;
//
//    return x.convert_to<double>();

}

int Context::cmp(cpp_int a, cpp_int b)
{
    if (a == b)
    {
        return 0;
    } else if (a > b) {
        return  1;
    } else {
        return -1;
    }
}

void Context::PolyToBigint(uint64_t* p1, cpp_int* coeffsBigint, long l)
{
    uint64_t  qi, level;
    level = l - 1;
    auto crtReconstruction = new cpp_int[level + 1]();

    cpp_int QiB;
    cpp_int tmp("1");
    cpp_int modulusBigint("1");


    for (uint64_t i = 0; i < level + 1; i++)
    {
        qi = qVec[i];
        QiB = cpp_int(qi);
        modulusBigint *= QiB;
        crtReconstruction[i] = cpp_int();
        crtReconstruction[i] = ModulusBigint / QiB;

        tmp = 1;
        for (uint64_t j = 0; j < L; ++j)
        {
            if (i != j)
            {
                auto imp = qInvModq[j][i];
                tmp *= cpp_int(imp);
                tmp %= qVec[i];
            }
        }

        tmp %= QiB;
        crtReconstruction[i] *= tmp;
    }
    for (uint64_t x = 0; x < N; x++)
    {
        tmp = cpp_int(uint64_t(0));
        for (uint64_t y = 0; y < level + 1; y++)
        {
            auto p1i = p1 + y * N;
            tmp = cpp_int(p1i[x]) * crtReconstruction[y];
            coeffsBigint[x] += tmp;
        }
        coeffsBigint[x] %= modulusBigint;
    }
}

// tested
cpp_int Context::ModInverse(cpp_int g, cpp_int n)
{
    cpp_int z;
    if (n < 0)
    {
        n = -n;
    }
    if (g < 0)
    {
        g = g % n;
    }
    cpp_int d, x, y;

    d = GCD(x, y, g, n);
    if (d != 1)
    {
        //printf("something wrong");
        return {};
    }

    if (x < 0)
    {
        z = x + n;
    } else {
        z = x;
    }

    return z;
}

// tested
cpp_int Context::GCD(cpp_int& x, cpp_int& y, cpp_int a, cpp_int b) {
    cpp_int z;
    if (a == 0 && b == 0) {
        x = y = 0;
        return 0;
    }

    // 处理特殊情况 a = 0，b != 0
    if (a == 0 && b != 0) {
        z = abs(b);
        x = 0;
        y = (b < 0) ? -1 : 1;
    }

    // 处理特殊情况 a != 0，b = 0
    if (a != 0 && b == 0) {
        z = abs(a);
        x = (a < 0) ? -1 : 1;
        y = 0;
    }

    // 计算最大公约数
    z = gcd(abs(a), abs(b));

    // 扩展欧几里得算法求解 x 和 y
    cpp_int old_x = 1;
    x = 0;
    cpp_int old_y = 0;
    y = 1;

    cpp_int temp;
    cpp_int quotient;

    while (b != 0) {
        quotient = a / b;

        temp = b;
        b = a - quotient * b;
        a = temp;

        temp = x;
        x = old_x - quotient * x;
        old_x = temp;

        temp = y;
        y = old_y - quotient * y;
        old_y = temp;
    }

    // 根据条件设置 x 和 y 的值
    if (x != 0 || y != 0) {
        x = old_x;
        y = old_y;
    }

    return z;
}

void Context::fullencode(uint64_t* a, complex<double>* v, long slots, long l) {
    fftSpecialInv(v, slots);
    long gap = Nh / slots;

//    if (l != L)
//    {
//        printf("Input Parameters Errors! I suggest you need to Guarantee l == L!");
//    }

    auto valuesfloat = new double[N]();
    for (uint64_t i = 0, jdx = Nh, idx = 0; i < slots; i = i + 1, jdx = jdx + gap, idx = idx + gap)
    {
        valuesfloat[idx] = v[i].real();
        valuesfloat[jdx] = v[i].imag();
    }

    auto subqVec = new uint64_t[l]();
    for (int i = 0; i < l; ++i)
    {
        subqVec[i] = qVec[i];
    }
    scaleUpVecExact(valuesfloat, p, subqVec, a, l);
    NTTAndEqual(a, l, 0);
    fftSpecial(v, slots);
}

void Context::scaleUpVecExact(double* values, double n, uint64_t *moduli, uint64_t* coeffs, long l) {
    bool isNegative;
    double xFlo;
    cpp_int xInt;
    cpp_int tmp;

    for (size_t i = 0; i < N; i++) {
        if (n * values[i] > 1.8446744073709552e+19) { // 2 << 64
            isNegative = false;
            if (values[i] < 0) {
                isNegative = true;
                xFlo = -n * values[i];
            } else {
                xFlo = n * values[i];
            }

            xFlo += 0.5;

            xInt = cpp_int(xFlo);

            for (size_t j = 0; j < l; j++) {
                tmp = xInt % moduli[j];
                auto coeffsi = coeffs + j * N;
                if (isNegative) {
                    coeffsi[i] = moduli[j] - tmp.convert_to<uint64_t>();
                } else {
                    coeffsi[i] = tmp.convert_to<uint64_t>();
                }
            }
        } else {
            if (values[i] < 0) {
                for (size_t j = 0; j < l; j++) {
                    auto coeffsi = coeffs + j * N;
                    coeffsi[i] = moduli[j] - (uint64_t)(-n * values[i] + double(0.5)) % moduli[j];
                }
            } else {
                for (size_t j = 0; j < l; j++) {
                    auto coeffsi = coeffs + j * N;
                    coeffsi[i] = (uint64_t)(n * values[i] + double(0.5)) % moduli[j];
                }
            }
        }
    }
}
