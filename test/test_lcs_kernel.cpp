#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "monge_matrix.h"
#include "lcs_kernel.h"

namespace LCS {
namespace kernel {
namespace {

void test_whole_a(const LCSKernel &kernel, const std::string &a, const std::string &b) {
    for (unsigned i = 0; i <= b.size(); ++i) {
        for (unsigned j = i; j <= b.size(); ++j) {
            ASSERT_EQ(kernel.lcs_whole_a(i, j), dp_lcs(a, b.substr(i, j - i)));
        }
    }
}

void test_whole_b(const LCSKernel &kernel, const std::string &a, const std::string &b) {
    for (unsigned i = 0; i <= a.size(); ++i) {
        for (unsigned j = i; j <= a.size(); ++j) {
            ASSERT_EQ(kernel.lcs_whole_b(i, j), dp_lcs(a.substr(i, j - i), b));
        }
    }
}

void test_lcs_suffix_a_prefix_b(const LCSKernel &kernel, const std::string &a, const std::string &b) {
    for (unsigned i = 0; i <= a.size(); ++i) {
        for (unsigned j = 0; j <= b.size(); ++j) {
            ASSERT_EQ(kernel.lcs_suffix_a_prefix_b(i, j), dp_lcs(a.substr(i), b.substr(0, j)));
        }
    }
}

void test_lcs_prefix_a_suffix_b(const LCSKernel &kernel, const std::string &a, const std::string &b) {
    for (unsigned i = 0; i <= a.size(); ++i) {
        for (unsigned j = 0; j <= b.size(); ++j) {
            ASSERT_EQ(kernel.lcs_prefix_a_suffix_b(i, j), dp_lcs(a.substr(0, i), b.substr(j)));
        }
    }
}

TEST(KernelTest, CalculateRecursiveLCSWholeFirstStringTest) {
    test_whole_a(RecursiveLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_whole_a(RecursiveLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_whole_a(RecursiveLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateRecursiveLCSWholeSecondStringTest) {
    test_whole_b(RecursiveLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_whole_b(RecursiveLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_whole_b(RecursiveLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateRecursiveLCSPrefixSuffixTest) {
    test_lcs_prefix_a_suffix_b(RecursiveLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_lcs_prefix_a_suffix_b(RecursiveLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_lcs_prefix_a_suffix_b(RecursiveLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateRecursiveLCSSuffixPrefixTest) {
    test_lcs_suffix_a_prefix_b(RecursiveLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_lcs_suffix_a_prefix_b(RecursiveLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_lcs_suffix_a_prefix_b(RecursiveLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateIterativeLCSWholeFirstStringTest) {
    test_whole_a(IterativeLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_whole_a(IterativeLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_whole_a(IterativeLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateIterativeLCSWholeSecondStringTest) {
    test_whole_b(IterativeLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_whole_b(IterativeLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_whole_b(IterativeLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateIterativeLCSPrefixSuffixTest) {
    test_lcs_prefix_a_suffix_b(IterativeLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_lcs_prefix_a_suffix_b(IterativeLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_lcs_prefix_a_suffix_b(IterativeLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateIterativeLCSSuffixPrefixTest) {
    test_lcs_suffix_a_prefix_b(IterativeLCS("BAABCBCA", "BAABCABCABACA"), "BAABCBCA", "BAABCABCABACA");
    test_lcs_suffix_a_prefix_b(IterativeLCS("xvuy", "uyxv"), "xvuy", "uyxv");
    test_lcs_suffix_a_prefix_b(IterativeLCS("AAAAAAAAAAA", "AAAAAAAAA"), "AAAAAAAAAAA", "AAAAAAAAA");
}


}  // namespace
}  // namespace matrix
}  // namespace LCS
