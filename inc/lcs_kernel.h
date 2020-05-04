#ifndef INC_LCS_KERNEL_H_
#define INC_LCS_KERNEL_H_

#include <iostream>
#include <string>
#include <algorithm>

#include "monge_matrix.h"

namespace LCS {
namespace kernel {

class LCSKernel {
public:
    // Initialize the LCS kernel for strings a and b.
    LCSKernel(const std::string &a, const std::string &b);
    // Count the lcs of the whole string a and the substring of b from b_l to b_r.
    unsigned lcs_whole_a(unsigned b_l, unsigned b_r);
    // Count the lcs of the substring of a from a_l to a_r and the whole string b.
    unsigned lcs_whole_b(unsigned a_l, unsigned a_r);
    // Count the lcs for the suffix of a from a_l and the prefix of b until b_r.
    unsigned lcs_suffix_a_prefix_b(unsigned a_l, unsigned b_r);
    // Count the lcs for the prefix of a until a_r and the suffix of b from b_l.
    unsigned lcs_prefix_a_suffix_b(unsigned a_r, unsigned b_l);
private:
    const std::string &a;
    const std::string &b;
    matrix::SubpermutationMatrix kernel;
    matrix::MongeMatrix kernel_sum;

    matrix::Permutation calculate_kernel(unsigned a_l, unsigned a_r,
                                         unsigned b_l, unsigned b_r);
};

}  // namespace kernel
}  // namespace LCS


#endif