/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_TESTSCHEME_H_
#define HEAANNTT_TESTSCHEME_H_

#include <iostream>

using namespace std;

class TestScheme {
public:
	static void testEncodeSingle(long logN, long L, long logp);

	static void testEncodeBatch(long logN, long L, long logp, long logSlots);

	static void testBasic(long logN, long L, long logp, long logSlots);

	static void testConjugateBatch(long logN, long L, long logp, long logSlots);

	static void testimultBatch(long logN, long L, long logp, long logSlots);

	static void testRotateByPo2Batch(long logN, long L, long logp, long logRotSlots, long logSlots, bool isLeft);

	static void testRotateBatch(long logN, long L, long logp, long rotSlots, long logSlots, bool isLeft);

	static void testSlotsSum(long logN, long L, long logp, long logSlots);

    //test ready-made func
    static void testConstAdd(long logN, long L, long logp, long logSlots);
    static void testComplexConstMult(long logN, long L, long logp, long logSlots);
    static void testComplexConstVecMult(long logN, long L, long logp, long logSlots);

    //test my basic func
    static void testmySquare(long logN, long L, long logp, long logSlots, long dnum=1);
    static void testmyMultAndEqual(long logN, long L, long logp, long logSlots, long dnum=1);
    static void testmyMult(long logN, long L, long logp, long logSlots, long dnum=1);
    static void testmyMult2(long logN, long L, long logp, long logSlots, long dnum=1);
    static void testmyMult3(long logN, long L, long logp, long logSlots, long dnum=1);
    static void DebugmyMult(long logN, long L, long logp, long logSlots, long dnum=1);
    static void testModUp(long logN, long L, long logp, long logSlots);
    static void testModDown(long logN, long L, long logp, long logSlots);
    static void testModUpAndDown(long logN, long L, long logp, long logSlots);

    static void KeySwitchingExtTrue(long logN, long L, long logp, long dnum, long logSlots);
    static void KeySwitchingExtFalse(long logN, long L, long logp, long dnum, long logSlots);

    static void testMultByIntegerInPlace();
    static void testConstMult(long logN, long L, long logp, long logSlots);

    static void testmyRescale(long logN, long L, long logp, long logSlots, long dnum);


    //----------------------------------------------------------------------------------
	//   POWER & PRODUCT TESTS
	//----------------------------------------------------------------------------------

	static void testPowerOf2Batch(long logN, long L, long logp, long logDegree, long logSlots);

	static void testPowerBatch(long logN, long L, long logp, long degree, long logSlots);

	static void testProdOfPo2Batch(long logN, long L, long logp, long logDegree, long logSlots);

	static void testProdBatch(long logN, long L, long logp, long degree, long logSlots);

    //----------------------------------------------------------------------------------
	//   FUNCTION TESTS
	//----------------------------------------------------------------------------------

	static void testInverseBatch(long logN, long L, long logp, long invSteps, long logSlots);

	static void testLogarithmBatch(long logN, long L, long logp, long degree, long logSlots);

	static void testExponentBatch(long logN, long L, long logp, long degree, long logSlots);

	static void testSigmoidBatch(long logN, long L, long logp, long degree, long logSlots);




};

#endif /* TESTSCHEME_H_ */
