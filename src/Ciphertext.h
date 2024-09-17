/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_CIPHERTEXT_H_
#define HEAANNTT_CIPHERTEXT_H_

#include "Common.h"

class Ciphertext {

public:

	uint64_t* bx;
	uint64_t* ax;

	long N; ///< Dimension of Ring

	long slots; ///< The length of plaintext vector

	long curr_limbs; ///< The level of this ciphertext

    long special_limbs; //The number of special limbs
    // Default constructor
    Ciphertext();

    // Constructor
    Ciphertext(uint64_t *ax, uint64_t *bx, long N, long slots, long l, long k=0);

    // My Constructor
    Ciphertext(long N, long slots, long l, long k=0);

    //My Destructor
    ~Ciphertext(void);

    // Copy constructor
    Ciphertext(const Ciphertext &cipher);

    Ciphertext &operator=(const Ciphertext &o);

};

#endif
