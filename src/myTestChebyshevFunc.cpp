//
// Created by EYx on 2023/4/30.
//
#include "Common.h"

#include "myTestChebyshevFunc.h"
#include "myChebyshevFunc.h"
#include "myUtils.h"

#include "TimeUtils.h"
#include "StringUtils.h"

//component for homomorphic function evaluation
void testLongDivisionChebyshev()
{
    //test1 以下测试数据和验证数据来自 EVAL LOGISTIC FUNCTION 中  EvalChebyshevSeriesPS 中的 auto divqr = LongDivisionChebyshev(f2, Tkm);
    double f2[] = {1,0.587681,-2.55872e-17,-0.121496,8.23994e-18,0.035318,-3.20924e-17,-0.0106942,6.00648e-17,0.00326191,1.51192e-16,-0.000998472,-1.3417e-16,0.000312841,-1.25036e-16,-0.000121295,0,0,0,0,0,1,};
    double Tkm[] = {0,0,0,0,0,0,0,0,0,0,0,0,1,};
    double gold_q[]={-2.6834e-16,0.000625682,-2.50071e-16,-0.00024259,0,0,0,0,0,2,};
    double gold_r[]={1,0.587681,-2.55872e-17,-1.1215,8.23994e-18,0.035318,-3.20924e-17,-0.0106942,6.00648e-17,0.00338321,2.76228e-16,-0.00131131,};
    double * quotient,*remainder;
    uint32_t quotient_len, remainder_len;
    LongDivisionChebyshev(quotient,quotient_len, remainder,remainder_len,
                          f2, 22,Tkm, 13);
    for (int i = 0; i < quotient_len; ++i) {
        if ((gold_q[i]-quotient[i])> 1e-5){
            cout<<"Key: "<<gold_q[i]<<"\t\t my result: "<<quotient[i]<<endl;
        }
    }

    for (int i = 0; i < remainder_len; ++i) {
        if ((gold_r[i]-remainder[i])> 1e-5){
            cout<<"Key: "<<gold_r[i]<<"\t\t my result: "<<remainder[i]<<endl;
        }
    }

//    //test2 以下测试数据和验证数据来自 EVAL LOGISTIC FUNCTION 中 EvalChebyshevSeriesPS 中的 auto divcs = LongDivisionChebyshev(r2, divqr->q);
//    double f[] = {1,0.587681,-2.55872e-17,-1.1215,8.23994e-18,0.035318,-3.20924e-17,-0.0106942,6.00648e-17,-0.996617,2.76228e-16,-0.00131131,};
//    double g[]={-2.6834e-16,0.000625682,-2.50071e-16,-0.00024259,0,0,0,0,0,2,};
//    double gold_r[] = {1,0.587993,-1.50429e-16,-1.12162,8.10948e-18,0.0353179,-3.20924e-17,-0.00938289,-2.16163e-16,};
//    double gold_q[] = {-0.996617,2.76228e-16,-0.00131131,};
//
//    double * quotient,*remainder;
//    uint32_t quotient_len, remainder_len;
//    LongDivisionChebyshev(quotient,quotient_len, remainder,remainder_len,
//                          f, 12, g, 10);
//
//
//    for (int i = 0; i < quotient_len; ++i) {
//        if ((gold_q[i] - quotient[i]) > 1e-5){
//            cout << "Key: " << gold_q[i] << "\t\t my result: " << quotient[i] << endl;
//        }
//    }
//
//    for (int i = 0; i < remainder_len; ++i) {
//        if ((gold_r[i]-remainder[i])> 1e-5){
//            cout<<"Key: "<<gold_r[i]<<"\t\t my result: "<<remainder[i]<<endl;
//        }
//    }

}

void testEvalLinearWSumMutable(long logN, long L, long logp, long logSlots, long dnum)
{
    cout << "!!! START TEST EvalLinearWSumMutable !!!" << endl;
    //-----------------------------------------
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);

    long slots =2;

    complex<double> mvec1[2]={0.25,0.5};
    complex<double> mvec2[2]={5,4};
    Ciphertext cipher1 = scheme.fullencrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.fullencrypt(mvec2, slots, L-2);

    Ciphertext wsum;
    uint32_t ciphertexts_num=2;
    Ciphertext *ciphertexts=new Ciphertext[ciphertexts_num];
    ciphertexts[0]=cipher1;
    ciphertexts[1]=cipher2;
    double  *constants=new double  [ciphertexts_num];
    constants[0]=2,constants[1]=4;

    EvalLinearWSumMutable(wsum,ciphertexts, ciphertexts_num, constants,scheme);

    complex<double>* dvec = scheme.decrypt(secretKey, wsum);

    complex<double> *golden_value = new complex<double> [ciphertexts_num];
    golden_value[0]=mvec1[0]*constants[0]+mvec2[0]*constants[1];
    golden_value[1]=mvec1[1]*constants[0]+mvec2[1]*constants[1];

    StringUtils::showcompare(golden_value, dvec, slots, "EvalLinearWSumMutable");

    if (dvec[0].real()-golden_value[0].real()>1e-10){
        cout<<"Check if dnum and myMult are both correctly set!"<<endl;
    }
}

void EvalLogisticExample()
{
    std::cout << "--------------------------------- EVAL LOGISTIC FUNCTION ---------------------------------"
              << std::endl;

    // Choosing a higher degree yields better precision, but a longer runtime.
    const uint32_t polyDegree = 16;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 6;

    long logN=10;
    long L=multDepth+1;
    long logp=55;   //fixme: check scaling modsize when test failed
    long dnum=7;

    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;

    SecretKey secretKey(context);
    DEBUG_SK="CURR_SK";
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    const long slots=8;//这个工程要求，slot必须是2的幂次！
    complex<double> input[slots]={-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    complex<double> goldenValue[slots]={0.0179885, 0.0474289, 0.119205, 0.268936, 0.5, 0.731064, 0.880795, 0.952571};


    Ciphertext cipher1 = scheme.fullencrypt(input, slots, L);

    double lowerBound = -5;
    double upperBound = 5;

//    auto result       = cc->EvalLogistic(ciphertext, lowerBound, upperBound, polyDegree);
//    Ciphertext result =EvalChebyshevFunction([](double x) -> double { return 1 / (1 + std::exp(-x)); }, cipher1, lowerBound, upperBound, polyDegree);
//    std::vector<double> coefficients = EvalChebyshevCoefficients(func, lowerBound, upperBound, polyDegree);
    double coefficients[polyDegree]={1, 0.587681, -2.55872e-17, -0.121496, 8.23994e-18, 0.035318, -3.20924e-17, -0.0106942, 6.00648e-17, 0.00326191, 1.51192e-16, -0.000998472, -1.3417e-16, 0.000312841, -1.25036e-16, -0.000121295};
    Ciphertext result;
    EvalChebyshevSeries(result,cipher1, coefficients, lowerBound, upperBound,polyDegree, scheme);


    complex<double>* dvec = scheme.decrypt(secretKey, result);
    StringUtils::showcompare(goldenValue, dvec, slots, "EvalLogisticExample");

}

void EvalFunctionExample()
{
    std::cout << "--------------------------------- EVAL SQUARE ROOT FUNCTION ---------------------------------"
              << std::endl;

    // Choosing a higher degree yields better precision, but a longer runtime.
    const uint32_t polyDegree = 50;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 7;

    long logN=10;
    long L=multDepth+1;
    long logp=55;   //fixme: check scaling modsize when test failed
    const long slots=8;//fixme: 这个工程要求，slot必须是2的幂次！
    long dnum=4;

    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;
    SecretKey secretKey(context);
    DEBUG_SK="CURR_SK";
    Scheme scheme(secretKey, context);

    //-----------------------------------------
    srand(time(NULL));
    //-----------------------------------------
    complex<double> input[slots]={1, 2, 3, 4, 5, 6, 7, 8,};
    complex<double> goldenValue_Square_Root[slots]={1, 1.414213, 1.732050, 2, 2.236067, 2.449489, 2.645751, 2.828427};
    complex<double> goldenValue_logx[slots]={0, 0.6931, 1.0986, 1.3862, 1.6094, 1.7917, 1.9459, 2.0794,};

    Ciphertext cipher1 = scheme.fullencrypt(input, slots, L);

    double lowerBound = 0;
    double upperBound = 10;

//    auto result_Square_Root       = cc->EvalLogistic(ciphertext, lowerBound, upperBound, polyDegree);
//    Ciphertext result_Square_Root =EvalChebyshevFunction([](double x) -> double { return 1 / (1 + std::exp(-x)); }, cipher1, lowerBound, upperBound, polyDegree);
//    std::vector<double> coefficients_Square_Root = EvalChebyshevCoefficients(func, lowerBound, upperBound, polyDegree);
    double coefficients_Square_Root[polyDegree]={4.0265, 1.34195, -0.268257, 0.114872, -0.0637436, 0.0405031, -0.0279885, 0.0204795, -0.0156203, 0.0122952,
                                                 -0.00991967, 0.00816345, -0.0068283, 0.00578941, -0.00496501, 0.00429971, -0.00375491, 0.00330303, -0.00292393, 0.00260266,
                                                 -0.00232789, 0.00209093, -0.00188503, 0.00170487, -0.00154619, 0.0014056, -0.00128032, 0.00116809, -0.00106704, 0.000975604,
                                                 -0.000892481, 0.000816571, -0.000746938, 0.000682785, -0.000623425, 0.000568266, -0.000516793, 0.000468553, -0.000423149, 0.000380231,
                                                 -0.000339485, 0.000300629, -0.000263409, 0.000227594, -0.000192971, 0.000159344, -0.000126527, 9.43477e-05, -6.26397e-05, 3.12425e-05,};
    Ciphertext result_Square_Root;
    EvalChebyshevSeries(result_Square_Root,cipher1, coefficients_Square_Root, lowerBound, upperBound, polyDegree, scheme);


    complex<double>* dvec_Square_Root = scheme.decrypt(secretKey, result_Square_Root);
    StringUtils::showcompare(goldenValue_Square_Root, dvec_Square_Root, slots, "EVAL SQUARE ROOT FUNCTION");

    double coefficients_logx[polyDegree]={
            1.86031,1.97227,-0.97226,0.638908,-0.472216,0.372184,-0.305477,0.257811,-0.222042,0.194202,
            -0.17191,0.15365,-0.138413,0.1255,-0.114409,0.104776,-0.0963248,0.0888455,-0.0821746,0.0761828,
            -0.0707668,0.0658427,-0.061342,0.0572079,-0.0533932,0.049858,-0.0465686,0.043496,-0.0406157,0.037906,
            -0.0353484,0.0329266,-0.0306261,0.0284343,-0.0263399,0.0243327,-0.0224036,0.0205444,-0.0187477,0.0170067,
            -0.015315,0.0136669,-0.012057,0.0104803,-0.008932,0.00740765,-0.00590295,0.00441379,-0.00293616,0.00146618,
    };
    Ciphertext result_logx;
    EvalChebyshevSeries(result_logx, cipher1,coefficients_logx,lowerBound,upperBound,polyDegree,scheme);
    complex<double>* dvec_logx = scheme.decrypt(secretKey, result_logx);
    StringUtils::showcompare(goldenValue_logx, dvec_logx, slots, "EVAL logx FUNCTION");
}

void testDoubleAngleIterations()
{
    std::cout << "--------------------------------- EVAL Double Angle Iterations  ---------------------------------"
              << std::endl;

    // Choosing a higher degree yields better precision, but a longer runtime.
    const uint32_t polyDegree = 16;

    // The multiplicative depth depends on the polynomial degree.
    // See the FUNCTION_EVALUATION.md file for a table mapping polynomial degrees to multiplicative depths.
    uint32_t multDepth = 14;

    long logN=10;
    long L=multDepth+1;
    long logp=55;   //fixme: check scaling modsize when test failed
    const long slots=8;//这个工程要求，slot必须是2的幂次！
    long dnum=1;

    TimeUtils timeutils;
    long K = L/dnum ; //Better, p10,Overview of Idea, Step0 //fixme: 是否需要上取整?, 判断参考OpenFHE
    Context context(logN, logp, L, K, dnum,512);
    cout<<"N: "<<context.N<<endl
        <<"L: "<<context.L<<endl
        <<"dnum: "<<context.dnum<<endl;

//    DEBUG_SK="Logistic";
//    SecretKey secretKey(context,DEBUG_SK);//fixme: 这是针对，logN=10,L=7,logp=55,dnum=1写死的私钥！
    SecretKey secretKey(context);
    DEBUG_SK="CURR_SK";
    Scheme scheme(secretKey, context);

    //------------
    //验证计算正确性
    //------------
    complex<double> input[slots]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
    complex<double> goldenValue[slots]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
    int32_t r = 6;
    for (int32_t j = 1; j < r + 1; j++) {
        for (int i = 0; i < slots; ++i) {
            goldenValue[i]*=goldenValue[i];
            goldenValue[i]*=2;
            double scalar = -1.0 / std::pow((2.0 * M_PI), std::pow(2.0, j - r));
            goldenValue[i]+=scalar;
        }
    }
    Ciphertext cipher1 = scheme.fullencrypt(input, slots, L);
    ApplyDoubleAngleIterations(cipher1,scheme);
    complex<double>* dvec = scheme.decrypt(secretKey, cipher1);
    StringUtils::showcompare(goldenValue, dvec, slots, "DoubleAngleIterations");

    //结合切比雪夫多项式验证功能。
//    complex<double> input[slots]={-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
//    complex<double> goldenValue[slots]={0.0179885, 0.0474289, 0.119205, 0.268936, 0.5, 0.731064, 0.880795, 0.952571};
//
//    Ciphertext cipher1 = scheme.fullencrypt(input, slots, L);
//
//    double lowerBound = -5;
//    double upperBound = 5;

//    auto result       = cc->EvalLogistic(ciphertext, lowerBound, upperBound, polyDegree);
//    Ciphertext result =EvalChebyshevFunction([](double x) -> double { return 1 / (1 + std::exp(-x)); }, cipher1, lowerBound, upperBound, polyDegree);
//    std::vector<double> coefficients = EvalChebyshevCoefficients(func, lowerBound, upperBound, polyDegree);
//    double coefficients[polyDegree]={1, 0.587681, -2.55872e-17, -0.121496, 8.23994e-18, 0.035318, -3.20924e-17, -0.0106942, 6.00648e-17, 0.00326191, 1.51192e-16, -0.000998472, -1.3417e-16, 0.000312841, -1.25036e-16, -0.000121295};
//    Ciphertext result = EvalChebyshevSeries(cipher1, coefficients, lowerBound, upperBound,polyDegree, scheme);
//    ApplyDoubleAngleIterations(result,scheme);
//
//    complex<double>* dvec = scheme.decrypt(secretKey, result);
//    StringUtils::showcompare(goldenValue, dvec, slots, "DoubleAngleIterations");

}