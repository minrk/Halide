#include "Halide.h"
#include <stdio.h>

using namespace Halide;

void check(Func f, ImageParam in, int min, int extent) {
    Image<int> output(12345);
    output.set_min(-1234);

    Buffer buf;
    in.set(buf);
    f.infer_input_bounds(output);
    Image<int> im = in.get();

    if (im.dim(0).extent() != extent || im.dim(0).min() != min) {
        printf("Inferred size was [%d, %d] instead of [%d, %d]\n",
               im.dim(0).min(), im.dim(0).extent(), min, extent);
        exit(-1);
    }
}

int main(int argc, char **argv) {
    ImageParam input(Int(32), 1);
    Var x;
    Func f1 = lambda(x, input(abs(cast<int8_t>(x))));
    Func f2 = lambda(x, input(abs(cast<int16_t>(x))));
    Func f3 = lambda(x, input(cast<int32_t>(abs(cast<float>(x)))));

    // input should be required from 0 to 128 inclusive, because abs
    // of an int 8 can return 128. This is an extent of 129.
    check(f1, input, 0, 129);
    check(f2, input, 0, 32769);

    // cast from int to float is treated as lossless, so we get 12345 - 1234
    check(f3, input, 0, 11111);

    // test a reflect boundary condition between zero and 100
    Expr reflect_x = 100 - cast<int>(abs(100 - (x % 200)));
    Func f4 = lambda(x, input(reflect_x));
    check(f4, input, 0, 101);

    // Verify an undefined bound on one side of the range still results in
    // correct bounds from abs and not an undefined error in the logic or
    // failure to bound the negative branch to zero.
    Func f5;
    f5 = lambda(x, input(cast<int>(clamp(abs(1.0f / (x + .1f)), -50, 50))));
    check(f5, input, 0, 51);

    printf("Success!\n");
    return 0;
}
