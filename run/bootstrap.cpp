/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <time.h>       /* time_t, struct tm, time, localtime, asctime */
#include <iostream>

using namespace std;

#include "../src/myBootstrapExample.h"

int main() {
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("\n\n"
           "The current date/time is: %s\n\n", asctime(timeinfo));

//当前正在测试，
//    BootstrapExampleClean(1 << 6, 2, 1);//需要使用N=64, slots=2的data文件
//    BootstrapExampleClean(1 << 6, 8, 1);//需要使用N=64, slots=8的data文件
//    BootstrapExampleClean(1 << 6, 16, 1);//需要使用N=64, slots=16的data文件
    BootstrapExampleClean(1 << 6, 32, 1);//需要使用N=64, slots=32的data文件

    return 0;
}
