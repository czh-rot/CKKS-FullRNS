/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_SECRETKEY_H_
#define HEAANNTT_SECRETKEY_H_

#include "Common.h"
#include "Context.h"
#include "constants.h"

class SecretKey {
public:

	uint64_t* sx;

	SecretKey(Context& context);
	SecretKey(Context& context, string STRING);


};


#endif
