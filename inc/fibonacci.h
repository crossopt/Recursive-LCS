#ifndef INC_FIBONACCI_H_
#define INC_FIBONACCI_H_

#include "lcs_kernel.h"
#include "monge_matrix.h"

namespace LCS {
namespace fkernel {

// Class that calculates the LCS kernel to solve the semi-local LCS problem
// for a plain pattern and a fibonacci string text.
class FibonacciKernel {
public:
    // Initialize the LCS kernel for pattern p and text F_{fn}.
    FibonacciKernel(const std::string &p, int fn);
    const std::string p;
    const int fn;
    std::pair <matrix::Permutation, matrix::Permutation> kernel;
private:

    std::pair <matrix::Permutation, matrix::Permutation> calculate_kernel(const std::string &p, int fn);
    std::pair <matrix::Permutation, matrix::Permutation> calculate_char_kernel(const std::string &p, char c);
};


}  // namespace fkernel
}  // namespace LCS

#endif
