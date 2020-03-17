#include "monge_matrix.h"
#include <iostream>

namespace LCS {
namespace matrix {

// Calculates the cross-difference of a Monge matrix.
// Throws MatrixException if the passed matrix was not
// a correct simple subunit-Monge matrix.
// (ie  if the resulting cross-difference is not a subpermutation).
// Time complexity is O(rows * cols). This function is currently being used
// in pair with the distribution-sum calculation that constructs and stores
// a rows * cols Monge matrix, so there is no need for optimization here.
SubpermutationMatrix::SubpermutationMatrix(
                    const MongeMatrix &distribution_matrix):
    MatrixInterface(distribution_matrix.get_rows() - 1,
                    distribution_matrix.get_cols() - 1) {
    matrix.assign(rows + 1, 0);
    std::vector <bool> was_used(cols + 1, false);
    for (unsigned i = 1; i < rows + 1; ++i) {
        for (unsigned j = 1; j < cols + 1; ++j) {
            unsigned cross_difference = distribution_matrix(i - 1, j)
                                      + distribution_matrix(i, j - 1)
                                      - distribution_matrix(i, j)
                                      - distribution_matrix(i - 1, j - 1);
            if (cross_difference != 0) {
                if (matrix[i] != 0 ||
                    cross_difference != 1 ||
                    was_used[j]) {
                    throw MatrixException("cross_difference calculation",
                            "cross_difference is not a subpermutation matrix");
                }
                matrix[i] = j;
                was_used[j] = true;
            }
        }
    }
}

// Calculates the distribution sum of a subpermutation.
// Time complexity is O(rows * cols) as the resulting rows * cols
// matrix is stored implicitly in the MongeMatrix class.
MongeMatrix::MongeMatrix(const SubpermutationMatrix &density_matrix):
        MatrixInterface(density_matrix.get_rows() + 1,
                        density_matrix.get_cols() + 1) {
    matrix.assign(rows, std::vector <unsigned> (cols, 0));

    for (unsigned i = rows - 1; i != 0; --i) {
        for (unsigned j = 1; j < cols; ++j) {
            matrix[i - 1][j] = density_matrix(i - 1, j - 1)
                             + matrix[i - 1][j - 1]
                             + matrix[i][j]
                             - matrix[i][j - 1];
        }
    }
}

unsigned SubpermutationMatrix::operator() (unsigned x, unsigned y) const {
    if (!is_element_correct(x, y)) {
        throw MatrixException("Subpermutation element query",
            "row " + std::to_string(x) + " column " + std::to_string(y));
    }
    return x + 1 < matrix.size() && matrix[x + 1] == y + 1;
}

unsigned MongeMatrix::operator() (unsigned x, unsigned y) const {
    if (!is_element_correct(x, y)) {
        throw MatrixException("MongeMatrix element query",
            "row " + std::to_string(x) + " column " + std::to_string(y));
    }
    return matrix[x][y];
}

SubpermutationMatrix SubpermutationMatrix::operator*(
                        const SubpermutationMatrix& m) const {
    MongeMatrix product = MongeMatrix(*this) * MongeMatrix(m);
    return SubpermutationMatrix(product);
}

// Multiplies two Monge matrices.
// No performance optimization is required as this is intended primarily for
// testing subpermutation sticky multiplication.
// All relevant operations with simple subunit-Monge matrices will be performed
// using their corresponding subpermutation matrices.
MongeMatrix MongeMatrix::operator*(const MongeMatrix& m) const {
    if (cols != m.get_rows()) {
        throw MatrixException("Tropical Multiplication",
                              "#arg1.columns != #arg2.rows");
    }
    std::vector<std::vector <unsigned>> result(rows,
                        std::vector<unsigned>(m.get_cols()));
    for (unsigned i = 0; i < rows; ++i) {
        for (unsigned k = 0; k < m.get_cols(); ++k) {
            result[i][k] = (*this)(i, 0) + m(0, k);
            for (unsigned j = 0; j < cols; ++j)
                result[i][k] = std::min(result[i][k], (*this)(i, j) + m(j, k));
        }
    }
    return MongeMatrix{result};
}

}  // namespace matrix
}  // namespace LCS
