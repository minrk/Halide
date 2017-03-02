// Halide tutorial lesson 18: Factoring an associative reduction using rfactor

// This lesson demonstrates how to parallelize or vectorize an associative
// reduction using the scheduling directive 'rfactor'.

// On linux, you can compile and run it like so:
// g++ lesson_18*.cpp -g -I ../include -L ../bin -lHalide -lpthread -ldl -o lesson_18 -std=c++11
// LD_LIBRARY_PATH=../bin ./lesson_18

// On os x:
// g++ lesson_18*.cpp -g -I ../include -L ../bin -lHalide -o lesson_18 -std=c++11
// DYLD_LIBRARY_PATH=../bin ./lesson_18

// If you have the entire Halide source tree, you can also build it by
// running:
//    make tutorial_lesson_18_parallel_associative_reductions
// in a shell with the current directory at the top of the halide
// source tree.

// The only Halide header file you need is Halide.h. It includes all of Halide.
#include "Halide.h"

// We'll also include stdio for printf.
#include <stdio.h>

using namespace Halide;

int main(int argc, char **argv) {
    // Declare some Vars to use below.
    Var x("x"), y("y"), i("i"), u("u"), v("v");

    // Create an input with random values.
    Buffer<uint8_t> input(8, 8, "input");
    for (int y = 0; y < 8; ++y) {
        for (int x = 0; x < 8; ++x) {
            input(x, y) = (rand() % 256);
        }
    }

    {
        // As mentioned previously in lesson 9, parallelizing variables that
        // are part of the reduction domain is tricky, since there may be data
        // dependencies across those variables.

        // Consider the histogram example in lesson 9:
        Func histogram("hist_serial");
        histogram(x) = 0;
        histogram.vectorize(x, 8);
        RDom r(0, input.width(), 0, input.height());
        histogram(input(r.x, r.y) % 8) += 1;

        histogram.realize(8);
        // See figures/lesson_18_hist_serial.gif for a visualization of
        // what this did.

        // Since there are data dependencies across r.x and r.y in the update
        // definition (i.e. the update refers to value computed in the previous
        // iteration), we can't parallelize r.x or r.y without introducing a race
        // condition. 8ote, however, that the histogram operation (i.e. the sum
        // reduction) is associative. A common trick to speed-up this type of
        // reduction is to split the reduction into chunks. Within a chunk, the
        // reduction is serial, however, since chunks are independent of
        // each other, we can parallelize across chunks.
    }

    {
        // Going back to the histogram example, we split the reduction into chunks
        // by defining an intermediate function:
        Func intermediate("intm_par_manual");
        intermediate(i, y) = 0;
        RDom rx(0, input.width());
        intermediate(input(rx, y) % 8, y) += 1;

        // This intermediate computes the sum reductions across the x dimension
        // (denoted by rx) for each y dimension independently. The histogram function
        // now sums over the partial results across the x dimension (denoted by i)
        // computed by the intermediate function:
        Func histogram("merge_par_manual");
        histogram(i) = 0;
        RDom ry(0, input.height());
        histogram(i) += intermediate(i, ry);

        // Since the intermediate no longer has data dependencies across the y dimension,
        // we can parallelize it over y:
        intermediate.compute_root().update().parallel(y);

        // Vectorize the initializations to make things faster.
        intermediate.vectorize(i, 8);
        histogram.vectorize(i, 8);

        histogram.realize(8);
        // See figures/lesson_18_hist_manual_par.gif for a visualization of
        // what this did.
    }

    {
        // This manual factorization of an associative reduction can be tedious and
        // bug-prone. Although it's fairly easy to do manually for histogram,
        // it can get complex pretty fast -- RDom may have predicates defined over
        // it or the value computed may be a multi-dimensional tuple. Halide
        // provides a way to do this type of factorization through the scheduling
        // directive 'rfactor'. rfactor splits an associative update definition
        // into an intermediate which computes the partial results over slices of a
        // reduction domain and replaces the current update definition with a new
        // definition which merges those partial results.

        // Using rfactor, we don't need to change the algorithm at all:
        Func histogram("hist_rfactor_par");
        histogram(x) = 0;
        RDom r(0, input.width(), 0, input.height());
        histogram(input(r.x, r.y) % 8) += 1;

        // The factoring of associative reduction is moved into the schedule,
        // via rfactor. rfactor takes as input a list of <RVar, Var> pairs,
        // which contains list of reduction variables (RVars) to be made
        // "parallelizable". In the generated intermediate Func, all references
        // to this reduction variables are replaced with references to "pure"
        // variables (the Vars). Since, by construction, Vars are race-condition
        // free, the intermediate reduction is now parallelizable across those
        // dimensions. All reduction variables not in the list are removed from
        // the original function and "lifted" to the intermediate.

        // To generate the same code as the manually-factored version, we do the
        // following:
        Func intermediate = histogram.update().rfactor({{r.y, y}});
        // We pass {r.y, y} as the argument to rfactor to make the histogram
        // parallelizable across the y dimension, similar to the manually-factored
        // version.
        intermediate.compute_root().update().parallel(y);

        // Since there is only one pair passed to rfactor, we could also write
        // the previous rfactor this way:
        // Func intermediate = histogram.update().rfactor(r.y, y);

        // Vectorize the initializations to make things faster.
        intermediate.vectorize(x, 8);
        histogram.vectorize(x, 8);

        // It is important to note that rfactor (or reduction factorization in
        // general) only works for associative reductions. Associative reductions
        // have the nice property that their results are the same no matter how
        // the computation is grouped (i.e. split into chunks). If rfactor can't
        // prove the associativity of a reduction, it will throw an error.

        Buffer<int> halide_result = histogram.realize(8);

        // See figures/lesson_18_hist_rfactor_par.gif for a visualization of
        // what this did.

        // The equivalent C is:
        int c_intm[input.height()][8];
        for (int y = 0; y < input.height(); y++) {
            for (int x = 0; x < 8; x++) {
                c_intm[y][x] = 0;
            }
        }
        /* parallel */ for (int y = 0; y < input.height(); y++) {
            for (int r_x = 0; r_x < input.width(); r_x++) {
                c_intm[y][input(r_x, y) % 8] += 1;
            }
        }

        int c_result[8];
        for (int x = 0; x < 8; x++) {
            c_result[x] = 0;
        }
        for (int x = 0; x < 8; x++) {
            for (int r_y = 0; r_y < input.height(); r_y++) {
                c_result[x] += c_intm[r_y][x];
            }
        }

        // Check the answers agree:
        for (int x = 0; x < 8; x++) {
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
        Func histogram("hist_rfactor_vec");
        histogram(x) = 0;
        RDom r(0, input.width(), 0, input.height());
        histogram(input(r.x, r.y) % 8) += 1;

        // Instead of r.y, we call rfactor on r.x this time:
        Func intermediate = histogram.update().rfactor(r.x, u);
        // which produces an intermediate that is race-condition free across the
        // x dimension (denoted by 'u' in the intermediate). This allows us to
        // vectorize the reduction across the inner dimension:
        intermediate.compute_root().update().vectorize(u, 8);
        // Note that since vectorizing the inner dimension changes the order of
        // computations, this trick only works if the associative reduction is
        // also commutative. rfactor will ensure these properties hold and will
        // throw an error if it can't prove those properties.

        // Vectorize the initializations to make things faster.
        intermediate.vectorize(x, 8);
        histogram.vectorize(x, 8);

        Buffer<int> halide_result = histogram.realize(8);

        // See figures/lesson_18_hist_rfactor_vec.gif for a visualization of
        // what this did.

        // The equivalent C is:
        int c_intm[input.width()][8];
        for (int u = 0; u < input.width(); u++) {
            for (int x = 0; x < 8; x++) {
                c_intm[u][x] = 0;
            }
        }
        for (int r_y = 0; r_y < input.height(); r_y++) {
            for (int u = 0; u < input.width() / 8; u++) {
                /* vectorize */ for (int u_i = 0; u_i < 8; u_i++) {
                    c_intm[u*4 + u_i][input(u*8 + u_i, r_y) % 8] += 1;
                }
            }
        }

        int c_result[8];
        for (int x = 0; x < 8; x++) {
            c_result[x] = 0;
        }
        for (int x = 0; x < 8; x++) {
            for (int r_x = 0; r_x < input.width(); r_x++) {
                c_result[x] += c_intm[r_x][x];
            }
        }

        // Check the answers agree:
        for (int x = 0; x < 8; x++) {
            if (c_result[x] != halide_result(x)) {
                printf("halide_result(%d) = %d instead of %d\n",
                       x, halide_result(x), c_result[x]);
                return -1;
            }
        }
    }

    {
        // Assume the input is really large and won't fit in the cache so we
        // need to tile the computation.
        Func histogram("hist_rfactor_tile");
        histogram(x) = 0;
        RDom r(0, input.width(), 0, input.height());
        histogram(input(r.x, r.y) % 8) += 1;

        // Let's first split both r.x and r.y by a factor of four.
        RVar rx_outer("rx_outer"), rx_inner("rx_inner");
        RVar ry_outer("ry_outer"), ry_inner("ry_inner");
        histogram.update().split(r.x, rx_outer, rx_inner, 4);
        histogram.update().split(r.y, ry_outer, ry_inner, 4);

        // To make this run faster, let's parallelize the reduction across
        // tiles. But since there are data dependencies across tiles, we
        // need to use rfactor() to factor the reduction:
        Func intermediate = histogram.update().rfactor({{rx_outer, u}, {ry_outer, v}});
        // which produces an intermediate that is race-condition free across
        // the tiles (i.e. across rx_outer and ry_outer). This allows us to
        // parallelize the reduction across the tiles:
        intermediate.compute_root().update().parallel(u).parallel(v);

        // Vectorize the initializations to make things faster.
        intermediate.vectorize(x, 8);
        histogram.vectorize(x, 8);

        Buffer<int> halide_result = histogram.realize(8);

        // See figures/lesson_18_hist_rfactor_tile.gif for a visualization of
        // what this did.

        // The equivalent C is:
        int c_intm[input.height() / 2][input.width() / 2][8];
        for (int v = 0; v < input.height() / 2; v++) {
            for (int u = 0; u < input.width() / 2; u++) {
                for (int x = 0; x < 8; x++) {
                    c_intm[v][u][x] = 0;
                }
            }
        }
        /* parallel */ for (int v = 0; v < input.height() / 2; v++) {
                for (int ry_inner = 0; ry_inner < 2; ry_inner++) {
                /* parallel */ for (int u = 0; u < input.width() / 2; u++) {
                    for (int rx_inner = 0; rx_inner < 2; rx_inner++) {
                        c_intm[v][u][input(u*2 + rx_inner, v*2 + ry_inner) % 8] += 1;
                    }
                }
            }
        }

        int c_result[8];
        for (int x = 0; x < 8; x++) {
            c_result[x] = 0;
        }
        for (int x = 0; x < 8; x++) {
            for (int ry_outer = 0; ry_outer < input.height() / 2; ry_outer++) {
                for (int rx_outer = 0; rx_outer < input.width() / 2; rx_outer++) {
                    c_result[x] += c_intm[ry_outer][rx_outer][x];
                }
            }
        }

        // Check the answers agree:
        for (int x = 0; x < 8; x++) {
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
