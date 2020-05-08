#ifndef INC_LCS_KERNEL_H_
#define INC_LCS_KERNEL_H_

#include "monge_matrix.h"

namespace LCS {
namespace kernel {

// Counts the lcs of two strings using the O(|a||b|) dynamic programming algorithm.
unsigned dp_lcs(const std::string &a, const std::string &b);

// Class that calculates the LCS kernel to solve the semi-local LCS problem.
class LCSKernel {
public:
    // Initialize the LCS kernel for strings a and b.
    LCSKernel(const std::string &a, const std::string &b, const matrix::MongeMatrix &kernel_sum);
    // Count the lcs of the whole string a and the substring of b from b_l to b_r.
    unsigned lcs_whole_a(unsigned b_l, unsigned b_r) const;
    // Count the lcs of the substring of a from a_l to a_r and the whole string b.
    unsigned lcs_whole_b(unsigned a_l, unsigned a_r) const;
    // Count the lcs for the suffix of a from a_l and the prefix of b until b_r.
    unsigned lcs_suffix_a_prefix_b(unsigned a_l, unsigned b_r) const;
    // Count the lcs for the prefix of a until a_r and the suffix of b from b_l.
    unsigned lcs_prefix_a_suffix_b(unsigned a_r, unsigned b_l) const;
protected:
    const std::string a;
    const std::string b;
    const matrix::MongeMatrix kernel_sum;
};

// Class that calculates the LCS kernel for two strings using the basic recursive algorithm.
// This allows it to solve the semi-local LCS problem for two strings a and b in O(|a||b|) time.
class RecursiveLCS: public LCSKernel {
public:
    // Initialize the LCS kernel for strings a and b.
    RecursiveLCS(const std::string &a, const std::string &b);
private:
    // Count the LCS kernel for two substrings of a and b recursively.
    matrix::Permutation calculate_kernel(const std::string &a, const std::string &b, 
                                         unsigned a_l, unsigned a_r,
                                         unsigned b_l, unsigned b_r);
};

// Class that calculates the LCS kernel for two strings using iterative combing.
// This allows it to solve the semi-local LCS problem for two strings a and b in O(|a||b|) time.
class IterativeLCS: public LCSKernel {
public:
    // Initialize the LCS kernel for strings a and b.
    IterativeLCS(const std::string &a, const std::string &b);
private:
    // Count the LCS kernel for two substrings of a and b using iterative combing.
    matrix::MongeMatrix calculate_iterative_kernel(const std::string &a, const std::string &b);
};


}  // namespace kernel
}  // namespace LCS

#endif
