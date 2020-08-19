#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "monge_matrix.h"
#include "lcs_kernel.h"
#include "grammar_compressed.h"

namespace LCS {
namespace gc {
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

GrammarCompressed gc_fib_string(unsigned length) {
    if (length == 0) {
        return GrammarCompressed(length + 1, 'A');
    } else if (length == 1) {
        return GrammarCompressed(length + 1, gc_fib_string(0), GrammarCompressed(0, 'B'));
    } else {
        GrammarCompressed prev = gc_fib_string(length - 1);
        return GrammarCompressed(length + 1, prev, *prev.first_symbol);
    }
}

TEST(GrammarCompressedTest, FibonacciLengthOneTest) {
    std::string s = "APATTERNA";
    auto gc_string = gc_fib_string(0);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(0));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(0)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciLengthTwoTest) {
    std::string s = "APATBTERNAB";
    auto gc_string = gc_fib_string(1);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(1));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(1)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciLengthThreeTest) {
    std::string s = "APATBTERNAB";
    auto gc_string = gc_fib_string(2);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(2));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(2)), fkernel.lcs);
}



TEST(GrammarCompressedTest, FibonacciSameFibStringTest) {
    std::string s = "ABAABABAABAAB";
    auto gc_string = gc_fib_string(5);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(5));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(5)), fkernel.lcs);
}


TEST(GrammarCompressedTest, FibonacciDifferentFibStringTest) {
    std::string s = "ABACABABDAABAAAB";
    auto gc_string = gc_fib_string(8);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(8));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(8)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciManyABsTest) {
    std::string s = "ABABABABABA";
    auto gc_string = gc_fib_string(5);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(5));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(5)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciNoMatchesTest) {
    std::string s = "QWERTY";
    auto gc_string = gc_fib_string(2);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(2));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(2)), fkernel.lcs);
}

TEST(GrammarCompressedTest, FibonacciOneMatchATest) {
    std::string s = "QWAERTY";
    auto gc_string = gc_fib_string(7);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(7));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(7)), fkernel.lcs);
}


TEST(GrammarCompressedTest, FibonacciOneMatchBTest) {
    std::string s = "QWERTBY";
    auto gc_string = gc_fib_string(7);
    auto fkernel = GCKernel(s, gc_string);
    ASSERT_EQ(gc_string.decompress(), fib_string(7));
    ASSERT_EQ(kernel::dp_lcs(s, fib_string(7)), fkernel.lcs);
}

TEST(GrammarCompressedTest, LZ78ComputesCorrectlyTest) {
    std::string s = "ABACABA";
    auto gc_lz78_string = LZ78(s);
    // Test the general structure of the compressed string.
    ASSERT_EQ(gc_lz78_string.decompress(), s);
    auto left = gc_lz78_string.first_symbol, right = gc_lz78_string.second_symbol;
    ASSERT_EQ(left->decompress(), "ABACAB");
    ASSERT_EQ(right->decompress(), "A");
    auto lower_left = left->first_symbol, lower_right = left->second_symbol;
    ASSERT_EQ(lower_left->decompress(), "ABAC");
    ASSERT_EQ(lower_right->decompress(), "AB");
    auto final_left = lower_left->first_symbol, final_right = lower_left->second_symbol;
    ASSERT_EQ(final_left->decompress(), "AB");
    ASSERT_EQ(final_right->decompress(), "AC");
    // Test that the decompression works correctly.
    ASSERT_EQ(gc_lz78_string.number, 11);
    ASSERT_EQ(left->number, 9);
    ASSERT_EQ(right->number, 10);
    ASSERT_EQ(lower_left->number, 6);
    ASSERT_EQ(lower_right->number, 8);
    ASSERT_EQ(final_left->number, 3);
    ASSERT_EQ(final_right->number, 5);
}

TEST(GrammarCompressedTest, LZWComputesCorrectlyTest) {
    std::string s = "ABACABA";
    auto gc_lzw_string = LZW(s);
    // Test the general structure of the compressed string.
    ASSERT_EQ(gc_lzw_string.decompress(), s);
    auto left = gc_lzw_string.first_symbol, right = gc_lzw_string.second_symbol;
    ASSERT_EQ(left->decompress(), "ABAC");
    ASSERT_EQ(right->decompress(), "ABA");
    auto lower_left = left->first_symbol, lower_right = left->second_symbol;
    ASSERT_EQ(lower_left->decompress(), "AB");
    ASSERT_EQ(lower_right->decompress(), "AC");
    // Test that the decompression works correctly.
    ASSERT_EQ(gc_lzw_string.number, 26 + 5);
    ASSERT_EQ(left->number, 26 + 3);
    ASSERT_EQ(right->number, 26 + 4);
    ASSERT_EQ(lower_left->number, 26 + 1);
    ASSERT_EQ(lower_right->number, 26 + 2);
}

}  // namespace
}  // namespace gc
}  // namespace LCS
