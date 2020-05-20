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
    for (auto i : fkernel.kernel.first.rows) std::cout << i.first << ' ' << i.second << '\n';
    for (auto i : fkernel.kernel.second.rows) std::cout << i.first << ' ' << i.second << '\n';
}

TEST(FibonacciTest, CalculateKernelLengthTwoTest) {
    std::string s = "APATBTERNAB";
    auto othera_kernel = kernel::IterativeLCS(s, "A");
    auto otherb_kernel = kernel::IterativeLCS(s, "B");
    auto other_kernel = kernel::IterativeLCS(s, "AB");
    auto fkernel = FibonacciKernel(s, 1);
    std::cout << "OK\n";
    for (auto i : fkernel.kernel.first.rows) std::cout << i.first << ' ' << i.second << '\n';
    std::cout << '\n';
    for (auto i : fkernel.kernel.second.rows) std::cout << i.first << ' ' << i.second << '\n';
}


}  // namespace
}  // namespace fkernel
}  // namespace LCS
