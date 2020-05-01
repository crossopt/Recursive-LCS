#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "gtest/gtest.h"
#include "monge_matrix.h"
#include "lcs.h"

namespace LCS {
namespace kernel {
namespace {

// void print_matrix(const matrix::MatrixInterface &matrix) {
//     for (unsigned i = 0; i < matrix.get_rows(); ++i) {
//         for (unsigned j = 0; j < matrix.get_cols(); ++j) {
//             std::cout << matrix(i, j) << ' ';
//             } std::cout << '\n';
//         }
// }

// void print_lcs_matrix(const matrix::MatrixInterface &matrix, int m) {
//     for (unsigned i = 0; i < matrix.get_rows(); ++i) {
//         for (unsigned j = 0; j < matrix.get_cols(); ++j) {
//             int value = j;
//             value += m;
//             value -= i;
//             value -= matrix(i, j);
//             std::cout << value << ' ';
//             } std::cout << '\n';
//         }
// }


void test_matrices_match(const matrix::MatrixInterface &lcs_matrix, int m) {
    std::cout << lcs_matrix.get_rows() << ' ' << lcs_matrix.get_cols() << '\n';
    for (unsigned i = 0; i < lcs_matrix.get_rows(); ++i) {
        for (unsigned j = 0; j < lcs_matrix.get_cols(); ++j) {
            int value = j;
            value -= i;
            value += m;
            value -= lcs_matrix(i, j);
            if (-10 < value && value < 10) {
                if (value < 0)
                    std::cout << value << "  ";// ':' << second_matrix(i, j) << ' ';
                else 
                    std::cout << value << "   ";// ':' << second_matrix(i, j) << ' ';
            } else {
                if (value < 0)
                    std::cout << value << " ";// ':' << second_matrix(i, j) << ' ';
                else 
                    std::cout << value << "  ";// ':' << second_matrix(i, j) << ' ';
            }
            // ASSERT_EQ(lcs_matrix(i, j), second_matrix(i, j)); 
        }
        std::cout << '\n';
    }
}


TEST(KernelTest, CalculateKernelTest) {
    // LCSKernel kernel = LCSKernel("AC", "AB");
    // test_matrices_match(kernel.kernel_sum, 2);
    // // print_matrix(matrix::SubpermutationMatrix(3, 3, {3, 2, 1}));
    // print_lcs_matrix(matrix::MongeMatrix(matrix::SubpermutationMatrix(3, 3, {1, 3, 2})), 1); // A AB
    // // print_lcs_matrix(matrix::MongeMatrix(matrix::SubpermutationMatrix(3, 3, {3, 1, 2})), 1); // C AB
    LCSKernel kernel = LCSKernel("BAABCBCA", "BAABCABCABACA");
    test_matrices_match(kernel.kernel_sum, 8);
}




}  // namespace
}  // namespace matrix
}  // namespace LCS
