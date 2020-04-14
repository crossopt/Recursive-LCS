#include <string>
#include <algorithm>
#include <iostream>

#include "gtest/gtest.h"
#include "monge_matrix.h"

namespace LCS {
namespace matrix {
namespace {

void test_matrices_match(const MatrixInterface &first_matrix,
                         const MatrixInterface &second_matrix) {
    ASSERT_EQ(first_matrix.get_rows(), second_matrix.get_rows());
    ASSERT_EQ(first_matrix.get_cols(), second_matrix.get_cols());
    for (unsigned i = 0; i < first_matrix.get_rows(); ++i) {
        for (unsigned j = 0; j < first_matrix.get_cols(); ++j) {
            // std::cout << first_matrix(i, j) << ':' << second_matrix(i, j) << ' ';
            ASSERT_EQ(first_matrix(i, j), second_matrix(i, j));
        }
        // std::cout <<'\n';
    }
}

// void test_cross_difference(unsigned rows,
//                            unsigned cols,
//                            const std::vector <std::vector <unsigned>> &original_matrix,
//                            const std::vector <unsigned> &expected_permutation) {
//     MongeMatrix matrix = MongeMatrix{original_matrix};
//     SubpermutationMatrix cross_difference = SubpermutationMatrix(matrix);
//     SubpermutationMatrix expected =
//             SubpermutationMatrix{rows - 1, cols - 1, expected_permutation};
//     test_matrices_match(cross_difference, expected);
// }

// void test_dominance_sum(unsigned rows,
//                         unsigned cols,
//                         const std::vector <unsigned> &original_permutation,
//                         const std::vector <std::vector <unsigned>> &expected_matrix) {
//     SubpermutationMatrix matrix = SubpermutationMatrix{rows, cols, original_permutation};
//     MongeMatrix dominance_sum = MongeMatrix(matrix);
//     MongeMatrix expected = MongeMatrix{expected_matrix};
//     test_matrices_match(dominance_sum, expected);
// }

// void test_monge_multiplication(const MongeMatrix &first_matrix,
//                                const MongeMatrix &second_matrix,
//                                const MongeMatrix &expected_product) {
//     MongeMatrix product = first_matrix * second_matrix;
//     test_matrices_match(product, expected_product);

//     SubpermutationMatrix permutation_product =
//             SubpermutationMatrix(first_matrix) ^ SubpermutationMatrix(second_matrix);
//     test_matrices_match(permutation_product, SubpermutationMatrix(expected_product));
// }

// TEST(MongeMatrixTest, DominanceSumAndCrossDifference) {
//     std::vector <unsigned> permutation = {2, 1, 3};  // {{0, 1, 0},
//                                                      //  {1, 0, 0},
//                                                      //  {0, 0, 1}};
//     std::vector<std::vector <unsigned>> monge_matrix = {{0, 1, 2, 3},
//                                                         {0, 1, 1, 2},
//                                                         {0, 0, 0, 1},
//                                                         {0, 0, 0, 0}};
//     test_dominance_sum(3, 3, permutation, monge_matrix);
//     test_cross_difference(4, 4, monge_matrix, permutation);
// }

// TEST(MongeMatrixTest, DominanceSumAndCrossDifferenceForSubpermutation) {
//     std::vector <unsigned> permutation = {2, 1};  // {{0, 1, 0},
//                                                   //  {1, 0, 0},
//                                                   //  {0, 0, 0}};
//     std::vector<std::vector <unsigned>> monge_matrix = {{0, 1, 2, 2},
//                                                         {0, 1, 1, 1},
//                                                         {0, 0, 0, 0},
//                                                         {0, 0, 0, 0}};
//     test_dominance_sum(3, 3, permutation, monge_matrix);
//     test_cross_difference(4, 4, monge_matrix, permutation);
// }

// TEST(MongeMatrixTest, DominanceSumAndCrossDifferenceNonSquareMatrixWithExtraCols) {
//     std::vector <unsigned> permutation = {2, 3};  // {{0, 1, 0},
//                                                   //  {0, 0, 1}};
//     std::vector<std::vector <unsigned>> monge_matrix = {{0, 0, 1, 2},
//                                                         {0, 0, 0, 1},
//                                                         {0, 0, 0, 0}};
//     test_dominance_sum(2, 3, permutation, monge_matrix);
//     test_cross_difference(3, 4, monge_matrix, permutation);
// }

// TEST(MongeMatrixTest, DominanceSumAndCrossDifferenceNonSquareMatrixWithExtraRows) {
//     std::vector <unsigned> permutation = {2, 3, 0, 1};  // {{0, 1, 0},
//                                                         //  {0, 0, 1},
//                                                         //  {0, 0, 0},
//                                                         //  {1, 0, 0}};
//     std::vector<std::vector <unsigned>> monge_matrix = {{0, 1, 2, 3},
//                                                         {0, 1, 1, 2},
//                                                         {0, 1, 1, 1},
//                                                         {0, 1, 1, 1},
//                                                         {0, 0, 0, 0}};
//     test_dominance_sum(4, 3, permutation, monge_matrix);
//     test_cross_difference(5, 4, monge_matrix, permutation);
// }

// TEST(MongeMatrixTest, DominanceSumAndCrossDifferenceForTrivialPermutation) {
//     std::vector <unsigned> permutation = {1, 2, 3, 4};  // {{1, 0, 0, 0},
//                                                         //  {0, 1, 0, 0},
//                                                         //  {0, 0, 1, 0},
//                                                         //  {0, 0, 0, 1}};
//     std::vector<std::vector <unsigned>> monge_matrix = {{0, 1, 2, 3, 4},
//                                                         {0, 0, 1, 2, 3},
//                                                         {0, 0, 0, 1, 2},
//                                                         {0, 0, 0, 0, 1},
//                                                         {0, 0, 0, 0, 0}};
//     test_dominance_sum(4, 4, permutation, monge_matrix);
//     test_cross_difference(5, 5, monge_matrix, permutation);
// }

// TEST(MongeMatrixTest, DominanceSumAndCrossDifferenceForZeroPermutation) {
//     std::vector <unsigned> permutation = {4, 3, 2, 1};  // {{0, 0, 0, 1},
//                                                         //  {0, 0, 1, 0},
//                                                         //  {0, 1, 0, 0},
//                                                         //  {1, 0, 0, 0}};
//     std::vector<std::vector <unsigned>> monge_matrix = {{0, 1, 2, 3, 4},
//                                                         {0, 1, 2, 3, 3},
//                                                         {0, 1, 2, 2, 2},
//                                                         {0, 1, 1, 1, 1},
//                                                         {0, 0, 0, 0, 0}};
//     test_dominance_sum(4, 4, permutation, monge_matrix);
//     test_cross_difference(5, 5, monge_matrix, permutation);
// }

// TEST(MongeMatrixTest, MongeMatrixMultiplicationSquaredMatrix) {
//     std::vector<std::vector <unsigned>> matrix = {{0, 1, 2, 3},
//                                                   {0, 1, 1, 2},
//                                                   {0, 0, 0, 1},
//                                                   {0, 0, 0, 0}};

//     MongeMatrix monge_matrix = MongeMatrix{matrix};
//     test_monge_multiplication(monge_matrix, monge_matrix, monge_matrix);
// }

// TEST(MongeMatrixTest, MongeMatrixMultiplicationByZero) {
//     std::vector<std::vector <unsigned>> matrix = {{0, 1, 2, 3},
//                                                   {0, 1, 1, 2},
//                                                   {0, 0, 0, 1},
//                                                   {0, 0, 0, 0}};

//     MongeMatrix monge_matrix = MongeMatrix{matrix};
//     MongeMatrix zero_matrix = MongeMatrix(SubpermutationMatrix{3, 3, {3, 2, 1}});
//     test_monge_multiplication(monge_matrix, zero_matrix, zero_matrix);
//     test_monge_multiplication(zero_matrix, monge_matrix, zero_matrix);
// }

// TEST(MongeMatrixTest, MongeMatrixMultiplicationByOne) {
//     std::vector<std::vector <unsigned>> matrix = {{0, 1, 2, 3},
//                                                   {0, 1, 1, 2},
//                                                   {0, 0, 0, 1},
//                                                   {0, 0, 0, 0}};

//     MongeMatrix monge_matrix = MongeMatrix{matrix};
//     MongeMatrix one_matrix = MongeMatrix(SubpermutationMatrix{3, 3, {1, 2, 3}});
//     test_monge_multiplication(monge_matrix, one_matrix, monge_matrix);
//     test_monge_multiplication(one_matrix, monge_matrix, monge_matrix);
// }

// TEST(MongeMatrixTest, MongeMatrixMultiplicationNonSquareLargerResult) {
//     std::vector<std::vector <unsigned>> first_matrix = {{0, 1, 2},
//                                                         {0, 1, 1},
//                                                         {0, 0, 0},
//                                                         {0, 0, 0}};
//     std::vector<std::vector <unsigned>> second_matrix = {{0, 1, 1, 2},
//                                                          {0, 0, 0, 1},
//                                                          {0, 0, 0, 0}};
//     std::vector<std::vector <unsigned>> expected_product = {{0, 1, 1, 2},
//                                                             {0, 1, 1, 1},
//                                                             {0, 0, 0, 0},
//                                                             {0, 0, 0, 0}};
//     test_monge_multiplication(MongeMatrix{first_matrix},
//                               MongeMatrix{second_matrix},
//                               MongeMatrix{expected_product});
// }

// TEST(MongeMatrixTest, MongeMatrixMultiplicationNonSquareSmallerResult) {
//     std::vector<std::vector <unsigned>> first_matrix = {{0, 1, 2},
//                                                         {0, 1, 1},
//                                                         {0, 0, 0},
//                                                         {0, 0, 0}};
//     std::vector<std::vector <unsigned>> second_matrix = {{0, 1, 1, 2},
//                                                          {0, 0, 0, 1},
//                                                          {0, 0, 0, 0}};
//     std::vector<std::vector <unsigned>> expected_product = {{0, 1, 1},
//                                                             {0, 0, 0},
//                                                             {0, 0, 0}};
//     test_monge_multiplication(MongeMatrix{second_matrix},
//                               MongeMatrix{first_matrix},
//                               MongeMatrix{expected_product});
// }

// TEST(MongeMatrixTest, SubpermutationMatrixMultiplicationForAllPermutations) {
//     std::vector <unsigned> first_permutation = {1, 2, 3, 4, 5};
//     do {
//         std::vector <unsigned> second_permutation = {1, 2, 3, 4, 5};
//         do {
//             SubpermutationMatrix actual_product =
//                     SubpermutationMatrix{5, 5, first_permutation} *
//                     SubpermutationMatrix{5, 5, second_permutation};
//             SubpermutationMatrix expected_product = 
//                     SubpermutationMatrix{5, 5, first_permutation} ^
//                     SubpermutationMatrix{5, 5, second_permutation};
//             test_matrices_match(actual_product, expected_product);
//         } while (std::next_permutation(second_permutation.begin(), second_permutation.end()));
//     } while (std::next_permutation(first_permutation.begin(), first_permutation.end()));
// }

TEST(MongeMatrixTest, SubpermutationMatrixMultiplicationForAllPermutations) {
    std::vector <unsigned> first_permutation = {1, 2, 3, 4, 5, 6};
    do {
        std::vector <unsigned> second_permutation = {1, 2, 3, 4, 5, 6};
        do {
            unsigned len = 6;
            // std::cout << "permutations " << '\n';
            // for (auto i : first_permutation) {
            //     std::cout << i << ' ';
            // } std::cout << "           ";
            // for (auto i : second_permutation) {
            //     std::cout << i << ' ';
            // } std::cout << "           \n";
            SubpermutationMatrix actual_product =
                    SubpermutationMatrix{len, len, first_permutation} *
                    SubpermutationMatrix{len, len, second_permutation};
            SubpermutationMatrix expected_product = 
                    SubpermutationMatrix{len, len, first_permutation} ^
                    SubpermutationMatrix{len, len, second_permutation};
            // std::cout << "permutations " << '\n';
            // for (auto i : first_permutation) {
            //     std::cout << i << ' ';
            // } std::cout << "           ";
            // for (auto i : second_permutation) {
            //     std::cout << i << ' ';
            // } std::cout << "           \n";
            test_matrices_match(actual_product, expected_product);
        } while (std::next_permutation(second_permutation.begin(), second_permutation.end()));
    } while (std::next_permutation(first_permutation.begin(), first_permutation.end()));
}

// TEST(MongeMatrixTest, SubpermutationMatrixMultiplicationLargeMatrix) {
//     std::vector <unsigned> first_permutation = {1, 2, 3, 4, 5, 6, 7};
//     std::vector <unsigned> second_permutation = {1, 2, 3, 4, 5, 6, 7};
//     SubpermutationMatrix actual_product =
//                     SubpermutationMatrix{7, 7, first_permutation} *
//                     SubpermutationMatrix{7, 7, second_permutation};
//     SubpermutationMatrix expected_product = 
//                     SubpermutationMatrix{7, 7, first_permutation} ^
//                     SubpermutationMatrix{7, 7, second_permutation};
//     test_matrices_match(actual_product, expected_product);
// }

// TEST(MongeMatrixTest, SubpermutationMatrixMultiplicationLargeMatrix) {
//     std::vector <unsigned> first_permutation = {1, 2, 3, 4};
//     std::vector <unsigned> second_permutation = {1, 2, 3, 4};
//     SubpermutationMatrix actual_product =
//                     SubpermutationMatrix{4, 4, first_permutation} *
//                     SubpermutationMatrix{4, 4, second_permutation};
//     SubpermutationMatrix expected_product = 
//                     SubpermutationMatrix{4, 4, first_permutation} ^
//                     SubpermutationMatrix{4, 4, second_permutation};
//     test_matrices_match(actual_product, expected_product);
// }

// TEST(MongeMatrixTest, SubpermutationMatrixMultiplicationTestMatrix) {
//     std::vector <unsigned> first_permutation = {2, 1};
//     std::vector <unsigned> second_permutation = {2, 1};
//     SubpermutationMatrix actual_product =
//                     SubpermutationMatrix{2, 2, first_permutation} *
//                     SubpermutationMatrix{2, 2, second_permutation};
//     SubpermutationMatrix expected_product = 
//                     SubpermutationMatrix{2, 2, first_permutation} ^
//                     SubpermutationMatrix{2, 2, second_permutation};
//     test_matrices_match(actual_product, expected_product);
// }

}  // namespace
}  // namespace matrix
}  // namespace LCS
