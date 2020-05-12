#include "lcs_kernel.h"

#include <iostream>
#include <algorithm>
#include <numeric>

namespace LCS {
namespace kernel {

unsigned dp_lcs(const std::string &a, const std::string &b) {
    std::vector <std::vector <unsigned>> lcs(a.size() + 1, std::vector <unsigned>(b.size() + 1, 0));
    for (unsigned i = 0; i < a.size(); ++i) {
        for (unsigned j = 0; j < b.size(); ++j) {
            lcs[i + 1][j + 1] = std::max(lcs[i][j + 1], lcs[i + 1][j]);
            if (a[i] == b[j]) {
                lcs[i + 1][j + 1] = std::max(lcs[i + 1][j + 1], lcs[i][j] + 1);
            }
        }
    }
    return lcs[a.size()][b.size()];
}

LCSKernel::LCSKernel(const std::string &a, const std::string &b, const matrix::MongeMatrix &kernel_sum): 
                                                                                    a(a),
                                                                                    b(b),
                                                                                    kernel_sum(kernel_sum) {}

RecursiveLCS::RecursiveLCS(const std::string &a, const std::string &b): LCSKernel(a, b,
    matrix::MongeMatrix(calculate_kernel(a, b, 0, a.size(), 0, b.size()).expand(a.size() + b.size(), a.size() + b.size()))) {}

IterativeLCS::IterativeLCS(const std::string &a, const std::string &b): LCSKernel(a, b, calculate_iterative_kernel(a, b)) {}

unsigned LCSKernel::lcs_whole_a(unsigned b_l, unsigned b_r) const {
    return b_r - b_l - kernel_sum(b_l + a.size(), b_r);
}

unsigned LCSKernel::lcs_whole_b(unsigned a_l, unsigned a_r) const {
    return b.size() - kernel_sum(a.size() - a_l, a.size() + b.size() - a_r);
}

unsigned LCSKernel::lcs_suffix_a_prefix_b(unsigned a_l, unsigned b_r) const {
    return b_r - kernel_sum(a.size() - a_l, b_r);
}

unsigned LCSKernel::lcs_prefix_a_suffix_b(unsigned a_r, unsigned b_l) const {
    return b.size() - b_l - kernel_sum(b_l + a.size(), a.size() + b.size() - a_r);
}

matrix::Permutation RecursiveLCS::calculate_recursion_base(const std::string &a,
                                                           const std::string &b,
                                                           unsigned a_l, unsigned a_r,
                                                           unsigned b_l, unsigned b_r) {
    std::vector <unsigned> last_row(b_r - b_l); //  The index of the braid strand at the end of row i.
    std::vector <unsigned> last_col(a_r - a_l); //  The index of the braid strand at the end of col i.
    std::iota(last_row.begin(), last_row.end(), a_r - a_l);
    for (unsigned i = a_l; i < a_r; ++i) {
        last_col[i] = a_r - i - 1;
        for (unsigned j = b_l; j < b_r; ++j) {
            // The braid strands should not cross if the string symbols match,
            // or if they have already crossed previously.
            // They have crossed previously if the natural ordering is ruined.
            if (a[i] == b[j] || last_col[i] > last_row[j]) {
                std::swap(last_col[i], last_row[j]);  // uncross the two strands
            }
        }
    }
    // Restore the kernel permutation from the braid.
    std::vector <unsigned> result(a_r - a_l + b_r - b_l);
    for (unsigned i = 0; i < a_r - a_l; ++i) {
        result[last_col[i]] = a_r - a_l + b_r - b_l - i;
    }
    for (unsigned i = 0; i < b_r - b_l; ++i) {
        result[last_row[i]] = i + 1;
    }
    return matrix::Permutation{result};
}

matrix::Permutation RecursiveLCS::calculate_kernel(const std::string &a,
                                                   const std::string &b,
                                                   unsigned a_l, unsigned a_r,
                                                   unsigned b_l, unsigned b_r) {
    unsigned sum_length = a_r - a_l + b_r - b_l;
    if (a_l >= a_r || b_l >= b_r) {
        return matrix::Permutation{{1}};
    } else if (a_l + 1 == a_r && b_l + 1 == b_r) {
        if (a[a_l] == b[b_l]) {
            return matrix::Permutation{{1, 2}};
        } else {
            return matrix::Permutation{{2, 1}};
        }
    } else if (a_l + 1 < a_r) {  // split by row
        unsigned a_m  = (a_l + a_r) / 2;
        matrix::Permutation first_half = calculate_kernel(a, b, a_l, a_m, b_l, b_r);
        matrix::Permutation second_half = calculate_kernel(a, b, a_m, a_r, b_l, b_r);
        first_half.grow_front(sum_length);
        second_half.grow_back(sum_length);
        return first_half * second_half;
    } else {  // first string has length 1 (split by column)
        unsigned b_m  = (b_l + b_r) / 2;
        matrix::Permutation first_half = calculate_kernel(a, b, a_l, a_r, b_l, b_m);
        matrix::Permutation second_half = calculate_kernel(a, b, a_l, a_r, b_m, b_r);
        first_half.grow_back(sum_length);
        second_half.grow_front(sum_length);
        return first_half * second_half;
    }
}

matrix::MongeMatrix IterativeLCS::calculate_iterative_kernel(const std::string &a, const std::string &b) {
    std::vector <unsigned> last_row(b.size()); //  The index of the braid strand at the end of row i.
    std::vector <unsigned> last_col(a.size()); //  The index of the braid strand at the end of col i.
    std::iota(last_row.begin(), last_row.end(), a.size());
    for (unsigned i = 0; i < a.size(); ++i) {
        last_col[i] = a.size() - i - 1;
        for (unsigned j = 0; j < b.size(); ++j) {
            // The braid strands should not cross if the string symbols match,
            // or if they have already crossed previously.
            // They have crossed previously if the natural ordering is ruined.
            if (a[i] == b[j] || last_col[i] > last_row[j]) {
                std::swap(last_col[i], last_row[j]);  // uncross the two strands
            }
        }
    }
    // Restore the kernel permutation from the braid.
    std::vector <unsigned> result(a.size() + b.size());
    for (unsigned i = 0; i < a.size(); ++i) {
        result[last_col[i]] = a.size() + b.size() - i;
    }
    for (unsigned i = 0; i < b.size(); ++i) {
        result[last_row[i]] = i + 1;
    }
    return matrix::MongeMatrix(matrix::PermutationMatrix(result.size(), result.size(), result));
}

}  // namespace kernel
}  // namespace LCS
