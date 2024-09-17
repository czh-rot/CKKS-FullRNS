/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "SecretKey.h"
#include "../data/testConstValue.h"
#include "../data/BSConstValue.h"

SecretKey::SecretKey(Context &context) {
    // 这里SK 默认是 (L+K) 大小！
    sx = new uint64_t[context.N * (context.L + context.K)]();
    context.sampleHWT(sx, context.L, context.K);
    context.NTTAndEqual(sx, context.L, context.K);
//    std::cout<<"cout sx"<<endl
//               <<"sx["<<context.N * (context.L + context.K)<<"]=";
//    for (int i = 0; i < context.N * (context.L + context.K); ++i) {
//        if ((i%10)==0){
//            std::cout<<std::endl;
//        }
//        std::cout<<sx[i]<< ",";
//    }
//    std::cout<<std::endl;
    DEBUG_SK = "NULL";
}

SecretKey::SecretKey(Context& context, string STRING) {
    // 这里SK 默认是 (L+K) 大小！

    //        sx[i]=SK_EVAL[i]; //debugModup时期的
    sx = new uint64_t[context.N * (context.L + context.K)]();
    context.sampleHWT(sx, context.L, context.K);
    context.NTTAndEqual(sx, context.L, context.K);

    if (STRING == "Rotation") {
        for (int i = 0; i < context.N * (context.L); i++) {
            sx[i] = Rotation_sk[i];
        }
    }
    else if (STRING == "Bootstrap") {
        for (int i = 0; i < (context.L + context.K); ++i) {
            uint64_t *sxj = sx + (i << context.logN);
            for (int j = 0; j < context.N; ++j) {
                sxj[j] = BS_sk[i][j];
            }
        }
    }
    else{
        cout<<"error! need to set sk for debugging"<<endl;
    }

}
