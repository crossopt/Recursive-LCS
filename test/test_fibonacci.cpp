#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "monge_matrix.h"
#include "lcs_kernel.h"
#include "fibonacci.h"

namespace LCS {
namespace fkernel {
namespace {


TEST(FibonacciTest, CalculateKernelLengthOneTest) {
    std::string s = "APATTERNA";
    auto other_kernel = kernel::IterativeLCS(s, "A");
    auto fkernel = FibonacciKernel(s, 0);
}

TEST(FibonacciTest, CalculateKernelLengthTwoTest) {
    std::string s = "APATBTERNAB";
    auto other_kernel = kernel::IterativeLCS(s, "AB");
    auto fkernel = FibonacciKernel(s, 1);
}

TEST(FibonacciTest, CalculateKernelLengthThreeTest) {
    std::string s = "APATBTERNAB";
    auto other_kernel = kernel::IterativeLCS(s, "ABA");
    auto fkernel = FibonacciKernel(s, 2);
}


}  // namespace
}  // namespace fkernel
}  // namespace LCS
