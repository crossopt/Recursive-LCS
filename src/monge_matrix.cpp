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

std::vector <std::pair<unsigned, unsigned>> Permutation::get_row_view() const {
    return row_view;
}

std::pair<Permutation, Permutation> Permutation::split_row() const {
    unsigned split_value = row_view[(row_view.size() - 1) / 2].first;
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

std::vector <std::pair<unsigned, unsigned>> Permutation::get_col_view() const {
    return col_view;
}

std::pair<Permutation, Permutation> Permutation::split_col() const {
    unsigned split_value = col_view[(col_view.size() - 1) / 2].first;
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

class SteadyAnt {
private:
    std::vector <std::pair <unsigned, unsigned>> good_elements_row;
    std::vector <std::pair <unsigned, unsigned>> good_elements_col;

    std::vector <std::pair <unsigned, unsigned>> low_cols ;  // seen are lower-right
    std::vector <std::pair <unsigned, unsigned>> low_rows ;  // seen are lower-right
    std::vector <std::pair <unsigned, unsigned>> high_cols;  // seen are upper-left
    std::vector <std::pair <unsigned, unsigned>> high_rows;  // seen are upper-left
    unsigned low_row_position = 0;
    unsigned low_col_position = 0;
    unsigned high_row_position = 0;
    unsigned high_col_position = 0;
    unsigned ant_row, ant_col;
    unsigned max_row, max_col;


    // Checks whether the ant has passed the last row with permutation elements.
    bool have_rows_ended() const {
        return low_row_position == low_rows.size() &&
               high_row_position == high_rows.size();
    }

    // Checks whether the ant has passed the last col with permutation elements.
    bool have_cols_ended() const {
        return low_col_position == low_cols.size() &&
               high_col_position == high_cols.size();
    }

    // Try to move the ant up.
    // It might stop seeing a bad R_high value, or see a new bad R_low value.
    // If this happens, the balance is broken and the ant can not move up.
    bool can_move_up() const {
        unsigned new_low_position = low_row_position;
        unsigned new_high_position = high_row_position;
        while (new_high_position < high_rows.size()
                && high_rows[new_high_position].first == ant_row) {
            if (high_rows[new_high_position].second < ant_col) {
                return false;
            }
            new_high_position++;
        }
        while (new_low_position < low_rows.size()
                && low_rows[new_low_position].first == ant_row) {
            if (low_rows[new_low_position].second >= ant_col) {
                return false;
            }
            new_low_position++;
        }
        return !have_rows_ended();
    }

    // Try to move the ant right.
    // It might stop seeing a bad R_low value, or see a new bad R_high value.
    // If this happens, the balance is broken and the ant can not move right.
    bool can_move_right() const {
        unsigned new_low_position = low_col_position;
        unsigned new_high_position = high_col_position;
        while (new_high_position < high_cols.size() 
                && high_cols[new_high_position].first == ant_col) {
            if (high_cols[new_high_position].second <= ant_row) {
                return false;
            }
            new_high_position++;
        }
        while (new_low_position < low_cols.size() 
                && low_cols[new_low_position].first == ant_col) {
            if (low_cols[new_low_position].second > ant_row) {
                return false;
            }
            new_low_position++;
        }
        return !have_cols_ended();
    }

    // Moves the ant up to the next row with permutation elements.
    void move_up() {
        while (high_row_position < high_rows.size() 
                && high_rows[high_row_position].first == ant_row) {
            if (high_rows[high_row_position].second >= ant_col) {
                good_elements_row.push_back(high_rows[high_row_position]);
            }
            high_row_position++;
        }
        while (low_row_position < low_rows.size()
                && low_rows[low_row_position].first == ant_row) {
            if (low_rows[low_row_position].second < ant_col) {
                good_elements_row.push_back(low_rows[low_row_position]);
            }
            low_row_position++;
        }
        ant_row = get_next_row();
    }

    // Moves the ant right to the next col with permutation elements.
    void move_right() {
        while (high_col_position < high_cols.size()
                && high_cols[high_col_position].first == ant_col) {
            if (high_cols[high_col_position].second > ant_row) {
                good_elements_col.push_back(high_cols[high_col_position]);
            }
            high_col_position++;
        }
        while (low_col_position < low_cols.size()
                && low_cols[low_col_position].first == ant_col) {
            if (low_cols[low_col_position].second <= ant_row) {
                good_elements_col.push_back(low_cols[low_col_position]);
            }
            low_col_position++;
        }
        ant_col = get_next_col();
    }

    // Returns the next interesting row (with permutation elements) for the ant,
    // or a fixed row smaller than all existing ones if all such rows have been visited.
    unsigned get_next_row() const {
        return std::max(low_row_position == low_rows.size()
                            ? max_row : low_rows[low_row_position].first, 
                        high_row_position == high_rows.size()
                            ? max_row : high_rows[high_row_position].first);
    }

    // Returns the next interesting column (with permutation elements) for the ant,
    // or a fixed col smaller than all existing ones if all such cols have been visited.
    unsigned get_next_col() const {
        return std::min(low_col_position == low_cols.size()
                            ? max_col : low_cols[low_col_position].first, 
                        high_col_position == high_cols.size()
                            ? max_col : high_cols[high_col_position].first);
    }
public:
    // Initializes the steady ant traversal for a pair of r_low, r_high permutation matrices.
    SteadyAnt(const Permutation &r_low, const Permutation &r_high) {
        low_cols = r_low.get_col_view();
        low_rows = r_low.get_row_view();
        high_cols = r_high.get_col_view();
        high_rows = r_high.get_row_view();
        // The ant starts from high rows (lower-left corner).
        std::reverse(low_rows.begin(), low_rows.end());
        std::reverse(high_rows.begin(), high_rows.end());

        max_row = std::min(low_rows.size() ? low_rows.back().first : 1, 
                           high_rows.size() ? high_rows.back().first : 1) - 1;
        max_col = std::max(low_cols.size() ? low_cols.back().first : 1, 
                           high_cols.size() ? high_cols.back().first : 1) + 1;
    }

    // Does the main ant traversal and returns the fixed permutation product.
    Permutation restore_correct_product() {
        // The ant position.
        // This is the pair of indexes before which the ant is currently located.
        ant_row = get_next_row();
        ant_col = get_next_col();
        while (!have_rows_ended() || !have_cols_ended()) {
            if (can_move_up()) {
                move_up();
            } else if (can_move_right()) {
                move_right();
            } else {
                // A diagonal move adds an element to the resulting permutation.
                good_elements_row.push_back({ant_row, ant_col});
                good_elements_col.push_back({ant_col, ant_row});
                move_up();
                move_right();
            }
        }
        // Fixes the reversed row order.
        std::reverse(good_elements_row.begin(), good_elements_row.end());
        return Permutation{good_elements_row, good_elements_col};
    }
};

// Recursively multiplies two permutations using the steady ant algorithm.
Permutation multiply(const Permutation &p, const Permutation &q) {
    if (p.get_nonzero_amount() == 0 || q.get_nonzero_amount() == 0) {
        // If all elements are zeroes, the product is a zero as well.
        return Permutation({}, {});
    }
    // The recursion base: at most one non-zero in each permutation.
    if (p.get_nonzero_amount() == 1 && q.get_nonzero_amount() == 1) {
        // The only non-zero value in the product C = A * B is 
        // C[i][k] if A[i][_] and B[_][k] are both non-zero.
        std::pair<unsigned, unsigned> p_nonzero = p.get_row_view()[0];
        std::pair<unsigned, unsigned> q_nonzero = q.get_col_view()[0];
        return Permutation({{p_nonzero.first, q_nonzero.first}},
                           {{q_nonzero.first, p_nonzero.first}});
    }
    // The divide phrase.
    // Split the first matrix by cols and the second by rows on the same index
    //  and remove all zeroes from the permutation halves (ie permutation pairs
    // from one half into the other).
    auto p_split = p.split_col();
    auto q_split = q.split_row();
    // Recursively multiply two pairs of permutations.
    Permutation r_low = multiply(p_split.first, q_split.first);
    Permutation r_high = multiply(p_split.second, q_split.second);

    // The conquer phrase.
    // "Ant" scanline counting the amount of wrong elements in the sum-product.
    // Remove all incorrect elements: higher than the ant scan for the one matrix,
    // lower than the ant scan for the other matrix.
    SteadyAnt ant = SteadyAnt(r_low, r_high);
    return ant.restore_correct_product();
}

SubpermutationMatrix SubpermutationMatrix::operator*(
                        const SubpermutationMatrix& m) const {
    if (cols != m.get_rows()) {
        throw MatrixException("Sticky Multiplication",
                              "#arg1.columns != #arg2.rows");
    }
    // Preprocessing before feeding everything to the recursive multiply function.
    // Convert the SubpermutationMatrix to the permutation format used by multiply.
    // It is not optimal for normal use as does not permit O(1) queries to elements
    // of the matrices, but it does not store zero rows and as such fits multiply.
    Permutation result = multiply(Permutation(*this), Permutation(m));
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
