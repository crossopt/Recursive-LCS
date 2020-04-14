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

// Multiplies two subpermutation matrices.
// No performance optimization is required as this is intended primarily for
// testing subpermutation sticky multiplication using the Steady Ant algorithm.
SubpermutationMatrix SubpermutationMatrix::operator^(
                        const SubpermutationMatrix& m) const {
    MongeMatrix product = MongeMatrix(*this) * MongeMatrix(m);
    return SubpermutationMatrix(product);
}

Permutation::Permutation(const SubpermutationMatrix &m) {
    std::vector <unsigned> count_for_cols(m.get_cols() + 1, 0);
    for (unsigned i = 0; i < m.get_rows(); ++i) {
        if (m.matrix[i + 1]) {
            row_view.push_back({i + 1, m.matrix[i + 1]});
        }
        count_for_cols[m.matrix[i + 1]] = i + 1;
    }
    for (unsigned i = 0; i < m.get_cols(); ++i) {
        if (count_for_cols[i + 1]) {
            col_view.push_back({i + 1, count_for_cols[i + 1]});
        }
    }
}

SubpermutationMatrix Permutation::expand(unsigned rows, unsigned cols) const {
    std::vector <unsigned> subpermutation(rows, 0);
    for (const auto &permutation_pair: row_view) {
        subpermutation[permutation_pair.first - 1] = permutation_pair.second;
    }
    return SubpermutationMatrix{rows, cols, subpermutation};
}

unsigned Permutation::get_row_split_index() const {
    return row_view[(row_view.size() - 1) / 2].first;
}

unsigned Permutation::get_col_split_index() const {
    return col_view[(col_view.size() - 1) / 2].first;
}

std::vector <std::pair<unsigned, unsigned>> Permutation::get_row_view() const {
    return row_view;
}

std::vector <std::pair<unsigned, unsigned>> Permutation::get_col_view() const {
    return col_view;
}

std::pair<Permutation, Permutation> Permutation::split_col(unsigned split_value) const {
    std::vector <std::pair <unsigned, unsigned> > row_first, row_second;
    std::vector <std::pair <unsigned, unsigned> > col_first, col_second;
    for (const auto &matched_pair: row_view) {
        if (matched_pair.second <= split_value) {
            row_first.push_back(matched_pair);
        }
        if (matched_pair.second > split_value) {
            row_second.push_back(matched_pair);
        }
    }
    for (const auto &matched_pair: col_view) {
        if (matched_pair.first <= split_value) {
            col_first.push_back(matched_pair);
        }
        if (matched_pair.first > split_value) {
            col_second.push_back(matched_pair);
        }
    }
    return {Permutation{row_first, col_first}, Permutation{row_second, col_second}};
}

std::pair<Permutation, Permutation> Permutation::split_row(unsigned split_value) const {
    std::vector <std::pair <unsigned, unsigned> > row_first, row_second;
    std::vector <std::pair <unsigned, unsigned> > col_first, col_second;
    for (const auto &matched_pair: row_view) {
        if (matched_pair.first <= split_value) {
            row_first.push_back(matched_pair);
        }
        if (matched_pair.first > split_value) {
            row_second.push_back(matched_pair);
        }
    }
    for (const auto &matched_pair: col_view) {
        if (matched_pair.second <= split_value) {
            col_first.push_back(matched_pair);
        }
        if (matched_pair.second > split_value) {
            col_second.push_back(matched_pair);
        }
    }
    return {Permutation{row_first, col_first}, Permutation{row_second, col_second}};
}

// How to join mappings: if joining for the product for two permutations, 
// where the first one's row x is actually i, and the other's row y is j.
// The resulting permutation square [x][y] should be [i][j].
Permutation steady_ant(const Permutation &p, const Permutation &q) {
    if (p.get_nonzero_amount() == 0 || q.get_nonzero_amount() == 0) {
        // If all elements are zeroes, the product is a zero as well.
        return Permutation({}, {});
    }
    // The recursion base: at most one non-zero in each permutation.
    if (p.get_nonzero_amount() == 1 && q.get_nonzero_amount() == 1) {
        // The only non-zero value in the product C = A * B is 
        // C[i][k] if A[i][j] and B[j][k] are both non-zero.
        std::pair<unsigned, unsigned> p_nonzero = p.get_row_view()[0];
        std::pair<unsigned, unsigned> q_nonzero = q.get_col_view()[0];
            return Permutation({{p_nonzero.first, q_nonzero.first}},
                               {{q_nonzero.first, p_nonzero.first}});
    }
    // The divide phrase.
    // Split the first matrix by cols and the second by rows on the same index
    //  and remove all zeroes from the permutation halves (ie permutation pairs
    // from one half into the other).
    auto p_split = p.split_col(p.get_col_split_index());
    auto q_split = q.split_row(q.get_row_split_index());
    // Recursively multiply two pairs of permutations.
    Permutation r_low = steady_ant(p_split.first, q_split.first);
    Permutation r_high = steady_ant(p_split.second, q_split.second);
    // Unmap the products by adding removed zero rows as zero rows in the result
    // from two matrices of size n / 2 to subpermutations of size n.
    // This is not required since the permutation is stored as an array of pairs.

    // The conquer phrase.
    // "Ant" scanline counting the amount of wrong elements in the sum-product.
    // Remove all incorrect elements: higher than the ant scan for the one matrix,
    // lower than the ant scan for the other matrix.
    SteadyAnt ant = SteadyAnt(r_low, r_high);
    ant.do_traversal();
    return Permutation(ant.good_elements_row, ant.good_elements_col);
}

SubpermutationMatrix SubpermutationMatrix::operator*(
                        const SubpermutationMatrix& m) const {
    if (cols != m.get_rows()) {
        throw MatrixException("Sticky Multiplication",
                              "#arg1.columns != #arg2.rows");
    }
    // Preprocessing before feeding everything to the recursive steady_ant function.
    // Convert the SubpermutationMatrix to the permutation format used by steady_ant.
    // It is not optimal for normal use as does not permit O(1) queries to elements
    // of the matrices, but it does not store zero rows and as such fits steady_ant.
    Permutation result = steady_ant(Permutation(*this), Permutation(m));
    return result.expand(rows, m.get_cols());
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
