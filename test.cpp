#include <stdio.h>

// #include "fixed.h"

// using namespace numeric;
// typedef Fixed<16, 16> fixed;


// int main(){
//     fixed f(0.1);
//     float b = (f * 2).to_float();
//     float c = (1 / f).to_float();
//     float d = (f - 0.1).to_float();
//     printf("%f, %f, %f", b, c, d);
//     return 0;
// }

#include "flexfloat.hpp"
#include <iostream>
#include <iomanip>

typedef flexfloat<6, 17> floatc;
// typedef floatx<6, 12> floatc;

int main(){
    
    floatc f_1 = 500.1;
    floatc f_2 = -2000;
    floatc f_3 = 0.001;
    float b = float(f_1 * f_2);
    float c = float(f_3 / f_1);
    float d = float(f_1 - f_1);
    printf("%f, %f, %f\n", b, c, d);
    return 0;
}


// int main(){
    
//     FP_LONG f_1 = FromFloat(500.1);
//     FP_LONG f_2 = FromFloat(-2000);
//     FP_LONG f_3 = FromFloat(0.001);
//     FP_LONG f_4 = FromFloat(2147480648);
//     FP_LONG f_5 = FromFloat(-2147480648);
//     float b = ToFloat(Mul(f_1, f_2));           // 500.1 * (-2000) = -1000200.00        (-1000200.00)       correct
//     float c = ToFloat(DivPrecise(f_3, f_1));    // 0.001 / 500.1 = 0.000002             (0.000002)          correct
//     float d = ToFloat(Sub(f_1,f_1));            // 500.1 - 500.1 = 0.000000             (0.000000)          correct

//     float e = ToFloat(Add(f_4,f_1));            // 2147480648 + 500.1 = 2147481148.1    (2147481216.000000) 
//     float f = ToFloat(DivPrecise(f_4,f_1));     // 2147480648 / 500.1 = 4294102.5       (4294102.500000)    correct
//     float g = ToFloat(Mul(f_4,f_1));            // 2147480648 * 500.1 = overflow        (213289184.000000)

//     float h = ToFloat(Add(f_5,f_1));            // -2147480648 + 500.1 = -2147480149.9  (-2147480192.000000)
//     float i = ToFloat(DivPrecise(f_5,f_1));     // -2147480648 / 500.1 = -4294102.5     (-4294102.500000)   correct
//     float j = ToFloat(Mul(f_5,f_1));            // -2147480648 * 500.1 = overflow       (-213289184.000000)
//     printf("%f, %f, %f\n%f, %f, %f\n%f, %f, %f\n", b, c, d, e, f, g, h, i, j);
//     return 0;


// }
