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

std::string fib_string(int n) {
	if (n == 0) {
		return "A";
	} else if (n == 1) {
		return "AB";
	} else {
		return fib_string(n - 1) + fib_string(n - 2);
	}
}

TEST(FibonacciTest, CalculateKernelLengthOneTest) {
    std::string s = "APATTERNA";
    auto fkernel = FibonacciKernel(s, 0);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(0)), fkernel.lcs);
}

TEST(FibonacciTest, CalculateKernelLengthTwoTest) {
    std::string s = "APATBTERNAB";
    auto fkernel = FibonacciKernel(s, 1);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(1)), fkernel.lcs);
}

TEST(FibonacciTest, CalculateKernelLengthThreeTest) {
    std::string s = "APATBTERNAB";
    auto fkernel = FibonacciKernel(s, 2);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(2)), fkernel.lcs);
}



TEST(FibonacciTest, CalculateKernelSameFibStringTest) {
    std::string s = "ABAABABAABAAB";
    auto fkernel = FibonacciKernel(s, 5);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(5)), fkernel.lcs);
}


TEST(FibonacciTest, CalculateKernelDifferentFibStringTest) {
    std::string s = "ABACABABDAABAAAB";
    auto fkernel = FibonacciKernel(s, 8);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(8)), fkernel.lcs);
}

TEST(FibonacciTest, CalculateKernelManyABsTest) {
    std::string s = "ABABABABABA";
    auto fkernel = FibonacciKernel(s, 5);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(5)), fkernel.lcs);
}

TEST(FibonacciTest, CalculateKernelNoMatchesTest) {
    std::string s = "QWERTY";
    auto fkernel = FibonacciKernel(s, 2);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(2)), fkernel.lcs);
}

TEST(FibonacciTest, CalculateKernelOneMatchATest) {
    std::string s = "QWAERTY";
    auto fkernel = FibonacciKernel(s, 7);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(7)), fkernel.lcs);
}


TEST(FibonacciTest, CalculateKernelOneMatchBTest) {
    std::string s = "QWERTBY";
    auto fkernel = FibonacciKernel(s, 7);
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(7)), fkernel.lcs);
}

}  // namespace
}  // namespace fkernel
}  // namespace LCS
