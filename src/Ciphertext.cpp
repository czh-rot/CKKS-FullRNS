/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"

Ciphertext::Ciphertext() : ax(nullptr), bx(nullptr), N(0), slots(0), curr_limbs(0) {}

Ciphertext::Ciphertext(uint64_t* ax, uint64_t* bx, long N, long slots, long l, long k)
: ax(ax), bx(bx), N(N), slots(slots), curr_limbs(l),special_limbs(k){}

Ciphertext::Ciphertext(long N, long slots, long l, long k)
        : ax(nullptr), bx(nullptr),N(N), slots(slots), curr_limbs(l), special_limbs(k)
        {
            long limbs=curr_limbs+special_limbs;
            ax = new uint64_t[N * limbs];
            bx = new uint64_t[N * limbs];
        }

Ciphertext::Ciphertext(const Ciphertext& cipher) : N(cipher.N), slots(cipher.slots),
curr_limbs(cipher.curr_limbs), special_limbs(cipher.special_limbs) {
    long limbs=curr_limbs+special_limbs;
	ax = new uint64_t[N * limbs];
	bx = new uint64_t[N * limbs];
	for (long i = 0; i < N * limbs; ++i) {
		ax[i] = cipher.ax[i];
		bx[i] = cipher.bx[i];
	}
}

Ciphertext& Ciphertext::operator=(const Ciphertext& o) {
	if(this == &o) return *this; // handling of self assignment, thanks for your advice, arul.
	delete[] ax;
	delete[] bx;
	N = o.N;
    curr_limbs = o.curr_limbs;
    special_limbs = o.special_limbs;
	slots = o.slots;
    long limbs=curr_limbs+special_limbs;
	ax = new uint64_t[N * limbs];
	bx = new uint64_t[N * limbs];
	for (long i = 0; i < N * limbs; ++i) {
		ax[i] = o.ax[i];
		bx[i] = o.bx[i];
	}
	return *this;
}
