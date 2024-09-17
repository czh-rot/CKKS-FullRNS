# Readme

## TODO
- [x] ~~Bootstrap 同态编解码相关的部分代码还没上传~~ 已上传 ckksrns-fhe.h、ckksrns-utils.h 和 ckksrns-utils.cpp，设计基本同openfhe
- [x] ~~确认 现在 openFHE中，当 rescale=FIXEDMANULA时，bootstrap是否正确运行~~ 可以正确运行

## General

- 函数参数中的 `logp` 表示每个 RNS模数 `qi` 的大小。第一个RNS模数的位宽需大一些，在工程中通过宏定义 `#define Q0_BIT_SIZE` 设置。slot_bit_size=Q0_BIT_SIZE-logp-1 即为同态计算数据的位宽，例如 Q0_BIT_SIZE=61, logp=55, 则slot_bit_size=5，最多可正确计算绝对值小于32的数据进行同态计算。
- 基于 FullRNS-HEAAN 开源库，参考 Over100x 写的 BasisConv，参考 Better Bootstrap 写的 hybrid KS，参考 openfhe 写的 func-eval 和 fastRotation
- curr_limbs表示当前模数的个数。例 curr_limbs=1 表示只剩 1 个 q0 不能再做乘法了。
- 本工程中 limbs 是做一次乘法减一的，所以密文的 limbs 统一到最小值后才可进行同态计算，与 openfhe 中统一到最大值是相反的。
- scheme 中的 `Ciphertext Scheme::multByConst(Ciphertext& cipher, double cnst) ` 函数对于某些 cnst值无法正确计算，暂时采取先解密再加密的方法绕过这个函数的实现。
- 文件名中my开头的都是我写的
- 注释可能在cpp中或者头文件相应声明处
- 目前参考 openfhe中 ckks-bootstrap 示例程序进行实现。（就是 ckks-bootstrap中的main）

## KeySwitch

函数命名遵循 openfhe keyswitch-hybrid.cpp
实现是参考 over100x 和 better bootstrapping 写的
KeySwitchExt：

- 如果 某次密文乘法后不做重线性导致这里出现三个分量的输入就会出bug。openfhe是可以解决任意个分量的密文的，但是这里没有实现这种通用性

## Rotation

函数 EvalFastRotationPrecompute、EvalFastRotation、EvalFastRotationExt 命名和实现遵循 openfhe 

- 名字里含有 fastRotate 词 的函数使用方式见相应 test 和 头文件中函数声明处的注释
    - 例如：FastRotate_demo 需要先调用 FastRotate_KeyGen，见 testFastRotate_demo 函数

## Debug

openFHE工程最外层的那个CMakeLists.txt中搜索 -O3，把 -O3 删掉，加 -g（可选）

cmake选项：生成extra示例：-DBUILD_EXTRAS=ON 

环境变量：关闭所有多线程便于debug：OMP_NUM_THREADS=1 

Debug可能需要写死

- 调用构造函数 Context(string STRING,  ..., ...) ，固定 moduliQ, rootsQ, moduliP, rootsP, 
- sk, swk的ax和bx,；openfhe使用sk加密
- ct的ax和bx

## 其它

- 若初始化 `class Ciphertext`，使用 `Ciphertext(long N, long slots, long l);` 并调用 `myUtils.cpp` 中的 `SetZero` 函数进行置零初始化。
- 若openfhe函数返回值为 `class Ciphertext`， 应将返回值修改为 函数的参数。Ciphertext类 的内存分配应尽量在调用者处完成，否则应标注 fixme 说明。
- 若函数中需要使用scheme类或者context类中的属性，应尽量在函数的开始声明同名本地变量进行调用。例如：`alpha=scheme.context.alpha`。
    - 若调用 Ciphertext类或者Plaintext类，需见机行事。（还没想好= =
- 对于vector，目前皆修改为 (数组，数组名_len) 这样的组合。例如：`int a_len=3; int* a = new [a_len];` 

- 尽量不创建新的类，新的函数不考虑写在scheme类、context类里，这两个类里既有的可以直接调用。
- 文件名、函数名尽量和openfhe保持一致
## 参考资料

- openfhe代码

   - 核心代码文件

      ckksrns-fhe.cpp

      keyswitch-hybrid.cpp

- openfhe论坛 

   - https://openfhe.discourse.group/t/reference-of-evalcoeffstoslots-function/438

- openfhe这个组织发的论文（解释了ckks 四种 rescale方式，本工程目前采取 FIXEDMANUAL

- 优化bootstrapping部分组件的论文
   1. Better bootstrapping for approximate homomorphic encryption【优化同态近似取模】
   2. Improved bootstrapping for approximate homomorphic encryption【优化同态编解码】
   3. Improved Homomorphic Discrete Fourier Transforms and FHE Bootstrapping【优化同态编解码】
   4. Efficient Bootstrapping for Approximate Homomorphic Encryption with Non-Sparse Keys【优化rescale】
