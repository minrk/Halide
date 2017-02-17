// Halide tutorial lesson 17: Factoring a reduction using rfactor

// This lesson demonstrates basic usage of Halide as a JIT compiler for imaging.

// On linux, you can compile and run it like so:
// g++ lesson_18*.cpp -g -I ../include -L ../bin -lHalide -lpthread -ldl -o lesson_18 -std=c++11
// LD_LIBRARY_PATH=../bin ./lesson_18

// On os x:
// g++ lesson_18*.cpp -g -I ../include -L ../bin -lHalide -o lesson_18 -std=c++11
// DYLD_LIBRARY_PATH=../bin ./lesson_18

// If you have the entire Halide source tree, you can also build it by
// running:
//    make tutorial_lesson_18_parallel_reductions
// in a shell with the current directory at the top of the halide
// source tree.

// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"

// We'll also include stdio for printf.
#include <stdio.h>

using namespace Halide;

// Support code for loading pngs.
#include "halide_image_io.h"
using namespace Halide::Tools;

int main(int argc, char **argv) {
    // Declare some Vars to use below.
    Var x("x"), y("y"), u("u");

    // Load a grayscale image to use as an input.
    Buffer<uint8_t> input = load_image("images/gray.png");

    {
        // As mentioned previously in lesson 9, parallelizing variables that
        // are part of the reduction domain is tricky, since there may be data
        // dependencies across those variables.
        //
        // Consider the histogram example in lesson 9:
        //
        // Func histogram("histogram");
        // Var x("x");
        // histogram(x) = 0;
        // RDom r(0, input.width(), 0, input.height());
        // histogram(input(r.x, r.y)) += 1;
        //
        // Since there are data dependencies across r.x and r.y in the update
        // definition (i.e. the update refers to value computed in the previous
        // iteration), we can't parallelize r.x or r.y without introducing race
        // condition. Note, however, the histogram operation (i.e. the sum reduction)
        // is associative. A common trick to speed-up this type of reduction is to
        // split the reduction into chunks. Within a chunk, the reduction is not
        // serial, however, since chunks are independent of each other, we can
        // parallelize in chunk granularity.
        //
        // Going back to the histogram example, we split the reduction into chunks
        // by defining an intermediate function:
        //
        // Func intm("intm");
        // Var i("i"), y("y");
        // intm(i, y) = 0;
        // RDom rx(0, input.width());
        // intm(input(rx, y), y) += 1;
        //
        // This intermediate computes the sum reductions across the x dimension
        // (denoted by rx) for each y dimension independently. The histogram function
        // now sums over the partial results across the x dimension (denoted by i)
        // computed by the intermediate function:
        //
        // Func histogram("histogram");
        // histogram(i) = 0;
        // RDom ry(0, input.height());
        // histogram(i) += intm(i, ry);
        //
        // Since the intermediate no longer has data dependencies across the y dimension,
        // we can parallelize it over y:
        // intm.compute_root().update().parallel(y);
        //
        // This manual factorization of an associative reduction can be tedious and
        // bug-prone. Although it's fairly easy to do manually for histogram,
        // it can get complex pretty fast -- RDom may have predicates defined over
        // it or the value computed may be a multi-dimensional tuple. Halide
        // provides a way to do this type of factorization through the scheduling
        // directive 'rfactor'. rfactor splits an associative update definition
        // into an intermediate which computes the partial results over slices of a
        // reduction domain and replace the current update definition with a new
        // definition which merges those partial results.

        // Using rfactor, we don't need to change the algorithm at all:
        Func histogram("histogram");
        histogram(x) = 0;
        RDom r(0, input.width(), 0, input.height());
        histogram(input(r.x, r.y)) += 1;

        // The factoring of associative reduction is moved into the schedule,
        // via rfactor. To generate the same code as the manually-factored
        // version, we do the following:
        Func intm = histogram.update().rfactor({{r.y, y}});
        // rfactor takes as input a list of <RVar, Var> pairs, which contains all
        // reduction variables to be removed from the original function and lifted
        // to the intermediate function. The remaining reduction variables are
        // made 'pure' in the intermediate function, which make them race-condition
        // free and hence are parallelizable.

        // Since there is only one pair passed to rfactor, we could also write it
        // this way:
        // Func intm = histogram.update().rfactor(r.y, y);

        // Similar to the manual version, the intermediate is parallelizable
        // across the y dimension:
        intm.compute_root().update().parallel(y);

        // It is important to note that rfactor (or reduction factorization in
        // general) only works for associative reductions. Associtive reduction
        // has a nice property in which it does not change the result no matter
        // how we group the computation (i.e. splitting into chunks). If rfactor
        // can't prove the associativity of a reduction, it will throw an error.

        Buffer<int> halide_result = histogram.realize(256);

        // The equivalent C is:
        int c_intm[input.height()][256];
        for (int y = 0; y < input.height(); y++) {
            for (int x = 0; x < 256; x++) {
                c_intm[y][x] = 0;
            }
        }
        /* parallel */ for (int y = 0; y < input.height(); y++) {
            for (int r_x = 0; r_x < input.width(); r_x++) {
                c_intm[y][input(r_x, y)] += 1;
            }
        }

        int c_result[256];
        for (int x = 0; x < 256; x++) {
            c_result[x] = 0;
        }
        for (int x = 0; x < 256; x++) {
            for (int r_y = 0; r_y < input.height(); r_y++) {
                c_result[x] += c_intm[r_y][x];
            }
        }

        // Check the answers agree:
        for (int x = 0; x < 256; x++) {
            if (c_result[x] != halide_result(x)) {
                printf("halide_result(%d) = %d instead of %d\n",
                       x, halide_result(x), c_result[x]);
                return -1;
            }
        }
    }

    {
        // Since we can factor associative reductions through the scheduling
        // directive 'rfactor', it is really easy to explore various factorization
        // strategies. Given the same serial histogram code:
        Func histogram("histogram");
        histogram(x) = 0;
        RDom r(0, input.width(), 0, input.height());
        histogram(input(r.x, r.y)) += 1;

        // Instead of r.y, we call rfactor on r.x this time:
        Func intm = histogram.update().rfactor(r.x, u);
        // which produces an intermediate that is race-condition free across the
        // x dimension (denoted by 'u' in the intermediate). This allows us to
        // vectorize the reduction across the inner dimension:
        intm.compute_root().update().vectorize(u, 8);
        // Note that since vectorizing the inner dimension changes the order of
        // computations, this trick only works if the associative reduction is
        // also commutative. rfactor will ensure these properties hold and will
        // throw an error if it can't prove those properties.

        Buffer<int> halide_result = histogram.realize(256);

        // The equivalent C is:
        int c_intm[input.width()][256];
        for (int u = 0; u < input.width(); u++) {
            for (int x = 0; x < 256; x++) {
                c_intm[u][x] = 0;
            }
        }
        for (int r_y = 0; r_y < input.height(); r_y++) {
            for (int u = 0; u < input.width() / 8; u++) {
                /* vectorize */ for (int u_i = 0; u_i < 8; u_i++) {
                    c_intm[u*8 + u_i][input(u*8 + u_i, r_y)] += 1;
                }
            }
        }

        int c_result[256];
        for (int x = 0; x < 256; x++) {
            c_result[x] = 0;
        }
        for (int x = 0; x < 256; x++) {
            for (int r_x = 0; r_x < input.width(); r_x++) {
                c_result[x] += c_intm[r_x][x];
            }
        }

        // Check the answers agree:
        for (int x = 0; x < 256; x++) {
            if (c_result[x] != halide_result(x)) {
                printf("halide_result(%d) = %d instead of %d\n",
                       x, halide_result(x), c_result[x]);
                return -1;
            }
        }
    }


    printf("Success!\n");

    return 0;
}
