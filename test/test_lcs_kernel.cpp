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

unsigned lcs(const std::string &a, const std::string &b) {
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

void test_whole_a(const std::string &a, const std::string &b) {
    LCSKernel kernel = LCSKernel(a, b);
    for (unsigned i = 0; i <= b.size(); ++i) {
        for (unsigned j = i; j <= b.size(); ++j) {
            ASSERT_EQ(kernel.lcs_whole_a(i, j), lcs(a, b.substr(i, j - i)));
        }
    }
}

void test_whole_b(const std::string &a, const std::string &b) {
    LCSKernel kernel = LCSKernel(a, b);
    for (unsigned i = 0; i <= a.size(); ++i) {
        for (unsigned j = i; j <= a.size(); ++j) {
            ASSERT_EQ(kernel.lcs_whole_b(i, j), lcs(a.substr(i, j - i), b));
        }
    }
}

void test_lcs_suffix_a_prefix_b(const std::string &a, const std::string &b) {
    LCSKernel kernel = LCSKernel(a, b);
    for (unsigned i = 0; i <= a.size(); ++i) {
        for (unsigned j = 0; j <= b.size(); ++j) {
            ASSERT_EQ(kernel.lcs_suffix_a_prefix_b(i, j), lcs(a.substr(i), b.substr(0, j)));
        }
    }
}

void test_lcs_prefix_a_suffix_b(const std::string &a, const std::string &b) {
    LCSKernel kernel = LCSKernel(a, b);
    for (unsigned i = 0; i <= a.size(); ++i) {
        for (unsigned j = 0; j <= b.size(); ++j) {
            ASSERT_EQ(kernel.lcs_prefix_a_suffix_b(i, j), lcs(a.substr(0, i), b.substr(j)));
        }
    }
}

TEST(KernelTest, CalculateLCSWholeFirstStringTest) {
    test_whole_a("BAABCBCA", "BAABCABCABACA");
    test_whole_a("xvuy", "uyxv");
    test_whole_a("AAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateLCSWholeSecondStringTest) {
    test_whole_b("BAABCBCA", "BAABCABCABACA");
    test_whole_b("xvuy", "uyxv");
    test_whole_b("AAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateLCSPrefixSuffixTest) {
    test_lcs_prefix_a_suffix_b("BAABCBCA", "BAABCABCABACA");
    test_lcs_prefix_a_suffix_b("xvuy", "uyxv");
    test_lcs_prefix_a_suffix_b("AAAAAAAAA", "AAAAAAAAA");
}

TEST(KernelTest, CalculateLCSSuffixPrefixTest) {
    test_lcs_suffix_a_prefix_b("BAABCBCA", "BAABCABCABACA");
    test_lcs_suffix_a_prefix_b("xvuy", "uyxv");
    test_lcs_suffix_a_prefix_b("AAAAAAAAA", "AAAAAAAAA");
}

}  // namespace
}  // namespace matrix
}  // namespace LCS
