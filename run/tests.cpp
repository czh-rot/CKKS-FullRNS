//
// Created by EYx on 2023/5/28.
//
/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "../src/TestScheme.h"
#include <time.h>       /* time_t, struct tm, time, localtime, asctime */

// 引入文件流库
#include <fstream>
#include <iostream>
using namespace std;
#include <sstream>// 引入字符串流库

#include "../src/myTestRotation.h"
#include "../src/myTestChebyshevFunc.h"
#include "../src/myBootstrapExample.h"
#include "../src/myTestBootstrap.h"

int main() {
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    printf ("\n\n"
            "The current date/time is: %s\n\n", asctime (timeinfo) );

    //Debug 区
    {
//todo: swk 生成有问题
//    TestScheme::DebugmyMult(6, 2, 55, 1,1);
//    TestScheme::DebugmyMult(6, 4, 55, 1,2);

//todo: 需要mvec和cnstvec各自取正负值，正负值的绝对值小于1和大于1的情况的正确性测试
//        TestScheme::testComplexConstVecMult(15, 11, 55, 1);

    }

//测试dnum的影响
    {
//Level=24
//    TestScheme::testmyMult(16, 24, 45, 3,24);
//    TestScheme::testmyMult(16, 24, 45, 3,12);
//    TestScheme::testmyMult(16, 24, 45, 3,6);
//    TestScheme::testmyMult(16, 24, 45, 3,4);
//    TestScheme::testmyMult(16, 24, 45, 3,1);
//Level=12
//    TestScheme::testmyMult(15, 12, 55, 3,12);
//    TestScheme::testmyMult(15, 12, 55, 3,6);
//    TestScheme::testmyMult(15, 12, 55, 3,4);
//    TestScheme::testmyMult(15, 12, 55, 3,3);
//    TestScheme::testmyMult(15, 12, 55, 3,2);
//    TestScheme::testmyMult(15, 12, 55, 3,1);
    }

//测试 KeySwitch 操作的正确性
    {
//    TestScheme::testmyMult(6, 2, 55, 1);
//    TestScheme::testmyMult(6, 2, 55, 1,2);
//    TestScheme::testmyMult(6, 4, 55, 1,2);
//    TestScheme::testmyMult(6, 4, 55, 1,4);

//    TestScheme::KeySwitchingExtTrue(15, 8, 55, 2, 3);
//    TestScheme::KeySwitchingExtFalse(15, 8, 55, 2, 3);
    }

//测试乘法、连续乘法的正确性
    {
//    TestScheme::testmyMult(6, 4, 55, 1,2);
//    TestScheme::testmyMult2(6, 4, 55, 1,2);
//    TestScheme::testmyMult3(6, 4, 55, 1,2);
//    TestScheme::testmyMult3(15, 12, 55, 1,4);
//    TestScheme::testmySquare(6, 4, 55, 1,2);
//    TestScheme::testmySquare(15, 12, 55, 2,4);
//    TestScheme::testmyMultAndEqual(15, 12, 55, 3,4);
    }

//测试同态计算切比雪夫多项式近似
    {
//    testLongDivisionChebyshev();
//    testEvalLinearWSumMutable(6, 4, 55, 1,1);
//    EvalLogisticExample();
//    EvalFunctionExample();
    }

//测试轮转操作的正确性 1. left, batch; 2. flat implementation; 3. use hybrid KeySwitch tech
    {
//    TestScheme::testmyFlatRotateBatch(6, 4, 55, 2, 5, 4, true);
//    TestScheme::testmyFlatRotateBatch(15, 4, 55, 1, 3, 4, true);
    }

//快速轮转相关测试
    {
//    testAutomorphismTransform();
//    testPrecomputeAutoMap();
//
//    testHRotate(6, 4, 55, 2, 3);//目前只支持 通过>0 的rotSlots值 表示左轮转
//
//    testFastRotate_demo(15, 6, 55, 2, 2, 4);
//    testFastRotate_demo(15, 6, 55, 2, -2, 4);
//    testEvalFastRotation(15, 8, 55, 2, 2, 3);
//    testEvalFastRotation(15, 8, 55, 2, -2, 3);
//
//    testEvalFastRotationExt(6, 4, 55, 1, 2, 3);
//    testEvalFastRotationExt(15, 8, 55, 2, 2, 3);
//    testEvalFastRotationExt(15, 8, 55, 2, -2, 3);
    }

//其它我用来理解工程提供的函数的功能的测试
//    TestScheme::testModUp(6, 2, 55, 1);
//    TestScheme::testModDown(6, 2, 55, 1);
//    TestScheme::testModUpAndDown(6, 2, 55, 1);//fixme: 加密->ModUp->ModDown->解密，failed
//    TestScheme::testModUpAndDown(15, 11, 55, 3);
//    TestScheme::testConstAdd(6, 4, 55, 1);
//    TestScheme::testConstMult(6, 4, 55, 1);
//    TestScheme::testConstMult(15, 11, 55, 1);//fixme: 如果cnstvec是double类型好像不能正确计算，见最上方debug区。
//    TestScheme::mytestExponentBatch();
//    TestScheme::testComplexConstMult(15, 11, 55, 1);

//-----------------------
//验证我写的功能函数的正确性
//-----------------------
//    TestScheme::testmySquare(15, 12, 55, 2, 4);
//    TestScheme::testmyMultAndEqual(15, 12, 55, 3, 4);
//    TestScheme::testmyMult(15, 12, 55, 3, 4);
//    TestScheme::testmyMult3(15, 12, 55, 1, 4);
//
//    testEvalLinearWSumMutable(6, 4, 55, 1, 1);
//    EvalLogisticExample();
//    EvalFunctionExample();
//    testDoubleAngleIterations();
//
//    testHRotate(6, 4, 55, 2, 3);//目前只支持 通过>0 的rotSlots值 表示左轮转
//    testFastRotate_demo(15, 6, 55, 2, 2, 4);
//    testFastRotate_demo(15, 6, 55, 2, -2, 4);
//    testEvalFastRotation(15, 8, 55, 2, 2, 3);
//    testEvalFastRotation(15, 8, 55, 2, -2, 3);
//    testEvalFastRotationExt(15, 8, 55, 2, 2, 3);
//    testEvalFastRotationExt(15, 8, 55, 2, -2, 3);
//
//    TestScheme::KeySwitchingExtTrue(15, 8, 55, 2, 3);
//    TestScheme::KeySwitchingExtFalse(15, 8, 55, 2, 3);
//
//    testModulusSwitch();
//    testPartialSum(6, 1, 55, 1, 1);
//    testPartialSum(6, 6, 55, 3, 4);
//    TestScheme::testConjugateBatch(15, 6, 55, 1);
//    testConjugate_demo(15, 6, 55, 2, 4);
//    TestScheme::testMultByIntegerInPlace();

//    TestScheme::testConstMult(13, 2, 30, 12);
//    TestScheme::testmyRescale(6, 2, 55, 3, 1);

//当前正在测试

    //原始工程提供的测试
//	TestScheme::testEncodeSingle(14, 1, 55);

//	TestScheme::testEncodeBatch(15, 6, 55, 3);

//  TestScheme::testBasic(15, 11, 55, 3);

//	TestScheme::testConjugateBatch(15, 6, 55, 1);

//	TestScheme::testRotateByPo2Batch(16, 26, 40, 1, 4, false);

//  TestScheme::testRotateBatch(15, 6, 55, 3, 4, true);

//	TestScheme::testimultBatch(16, 16, 55, 2);

//	TestScheme::testPowerOf2Batch(16, 15, 50, 2, 3);

//	TestScheme::testInverseBatch(14, 5, 55, 4, 3);

//	TestScheme::testExponentBatch(14, 5, 55, 7, 3);
//
//	TestScheme::testSigmoidBatch(16, 15, 55, 3, 3);

//	TestScheme::testSlotsSum(16, 15, 40, 3);

//	TestScheme::testMeanVariance(14, 3, 55, 13);

//	TestScheme::testHEML("data/uis.txt", 0, 5);

    return 0;
}