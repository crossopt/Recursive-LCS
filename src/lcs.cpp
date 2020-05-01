#include "lcs.h"

#include <utility>
#include <vector>
#include <exception>

namespace LCS {
namespace kernel {

LCSKernel::LCSKernel(const std::string &a, const std::string &b): a(a), b(b),
    kernel(calculate_kernel(0, a.size(), 0, b.size())), kernel_sum(matrix::MongeMatrix(kernel)) {}

matrix::SubpermutationMatrix LCSKernel::calculate_kernel(unsigned a_l, unsigned a_r,
                                                        unsigned b_l, unsigned b_r) {
    unsigned sum_length = a_r - a_l + b_r - b_l;
    if (a_l >= a_r || b_l >= b_r) {
        return matrix::SubpermutationMatrix{1, 1, {1}};
    } else if (a_l + 1 == a_r && b_l + 1 == b_r) {
        if (a[a_l] == b[b_l]) {
            return matrix::SubpermutationMatrix{2, 2, {1, 2}};
        } else {
            return matrix::SubpermutationMatrix{2, 2, {2, 1}};
        }
    } else if (a_l + 1 < a_r) {  // split by row
        unsigned a_m  = (a_l + a_r) / 2;
        matrix::SubpermutationMatrix first_half = calculate_kernel(a_l, a_m, b_l, b_r);
        matrix::SubpermutationMatrix second_half = calculate_kernel(a_m, a_r, b_l, b_r);
        first_half.grow_front(sum_length);
        second_half.grow_back(sum_length);
        return first_half * second_half;
    } else { // first string has length 1 (split by column)
        unsigned b_m  = (b_l + b_r) / 2;
        matrix::SubpermutationMatrix first_half = calculate_kernel(a_l, a_r, b_l, b_m);
        matrix::SubpermutationMatrix second_half = calculate_kernel(a_l, a_r, b_m, b_r);
        first_half.grow_back(sum_length);
        second_half.grow_front(sum_length);
        return first_half * second_half;
    }
}

}  // namespace kernel
}  // namespace LCS
