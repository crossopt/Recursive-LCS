#include "lcs_kernel.h"

#include <utility>
#include <vector>
#include <exception>

namespace LCS {
namespace kernel {

LCSKernel::LCSKernel(const std::string &a, const std::string &b): a(a), b(b),
    kernel(calculate_kernel(0, a.size(), 0, b.size()).expand(a.size() + b.size(), a.size() + b.size())),
    kernel_sum(matrix::MongeMatrix(kernel)) {}

matrix::Permutation LCSKernel::calculate_kernel(unsigned a_l, unsigned a_r,
                                                    unsigned b_l, unsigned b_r) {
    unsigned sum_length = a_r - a_l + b_r - b_l;
    if (a_l >= a_r || b_l >= b_r) {
        return matrix::Permutation{{{1, 1}}, {{1, 1}}};
    } else if (a_l + 1 == a_r && b_l + 1 == b_r) {
        if (a[a_l] == b[b_l]) {
            return matrix::Permutation{{{2, 2}, {1, 1}}, {{1, 1}, {2, 2}}};
        } else {
            return matrix::Permutation{{{2, 1}, {1, 2}}, {{1, 2}, {2, 1}}};
        }
    } else if (a_l + 1 < a_r) {  // split by row
        unsigned a_m  = (a_l + a_r) / 2;
        matrix::Permutation first_half = calculate_kernel(a_l, a_m, b_l, b_r);
        matrix::Permutation second_half = calculate_kernel(a_m, a_r, b_l, b_r);
        first_half.grow_front(sum_length);
        second_half.grow_back(sum_length);
        return first_half * second_half;
    } else { // first string has length 1 (split by column)
        unsigned b_m  = (b_l + b_r) / 2;
        matrix::Permutation first_half = calculate_kernel(a_l, a_r, b_l, b_m);
        matrix::Permutation second_half = calculate_kernel(a_l, a_r, b_m, b_r);
        first_half.grow_back(sum_length);
        second_half.grow_front(sum_length);
        return first_half * second_half;
    }
}

unsigned LCSKernel::lcs_whole_a(unsigned b_l, unsigned b_r) {
    return b_r - b_l - kernel_sum(b_l + a.size(), b_r);
}

unsigned LCSKernel::lcs_whole_b(unsigned a_l, unsigned a_r) {
    return b.size() - kernel_sum(a.size() - a_l, a.size() + b.size() - a_r);
}

unsigned LCSKernel::lcs_suffix_a_prefix_b(unsigned a_l, unsigned b_r) {
    return b_r - kernel_sum(a.size() - a_l, b_r);
}

unsigned LCSKernel::lcs_prefix_a_suffix_b(unsigned a_r, unsigned b_l) {
    return b.size() - b_l - kernel_sum(b_l + a.size(), a.size() + b.size() - a_r);
}


}  // namespace kernel
}  // namespace LCS
