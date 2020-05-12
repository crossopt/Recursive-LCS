#include "monge_matrix.h"

namespace LCS {
namespace matrix {

// Calculates the cross-difference of a Monge matrix.
// Throws MatrixException if the passed matrix was not
// a correct simple subunit-Monge matrix.
// (ie  if the resulting cross-difference is not a subpermutation).
// Time complexity is O(rows * cols). This function is currently being used
// in pair with the distribution-sum calculation that constructs and stores
// a rows * cols Monge matrix, so there is no need for optimization here.
PermutationMatrix::PermutationMatrix(
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
MongeMatrix::MongeMatrix(const PermutationMatrix &density_matrix):
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

unsigned PermutationMatrix::operator() (unsigned x, unsigned y) const {
    if (!is_element_correct(x, y)) {
        throw MatrixException("Subpermutation element query",
            "row " + std::to_string(x) + " column " + std::to_string(y));
    }
    return x + 1 < matrix.size() && matrix[x + 1] == y + 1;
}

// Grows the permutation to the required size by adding trivial elements to its end.
// Takes O(new_size) time.
void Permutation::grow_back(unsigned new_cols) {
    unsigned max_row = rows[0].first;
    unsigned max_col = cols.back().first;
    if (new_cols <= max_col) {
        throw MatrixException("Growth query",
            "new column size " + std::to_string(new_cols));
    }
    rows.insert(rows.begin(), new_cols - max_col, {0, 0});
    for (unsigned add = 1; add <= new_cols - max_col; ++add) {
        cols.push_back({max_col + add, max_row + add});
        rows[new_cols - max_col - add] = {max_row + add, max_col + add};
    }
}

// Grows the permutation to the required size by adding trivial elements to its front.
// Takes O(new_size) time.
void Permutation::grow_front(unsigned new_rows) {
    unsigned max_row = rows[0].first;
    if (new_rows <= max_row) {
        throw MatrixException("Growth query",
            "new row size " + std::to_string(new_rows));
    }
    for (auto &i: rows) {
        i.first += (new_rows - max_row);
        i.second += (new_rows - max_row);
    }
    for (auto &i: cols) {
        i.first += (new_rows - max_row);
        i.second += (new_rows - max_row);
    }
    cols.insert(cols.begin(), new_rows - max_row, {0, 0});
    for (unsigned add = new_rows - max_row; add > 0; --add) {
        cols[add - 1] = {add, add};
        rows.push_back({add, add});
    }
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
PermutationMatrix PermutationMatrix::operator^(
                        const PermutationMatrix& m) const {
    MongeMatrix product = MongeMatrix(*this) * MongeMatrix(m);
    return PermutationMatrix(product);
}

Permutation::Permutation(const PermutationMatrix &m) {
    std::vector <unsigned> count_for_cols(m.get_cols() + 1, 0);
    for (unsigned i = m.get_rows(); i != 0; --i) {
        if (m.matrix[i]) {
            rows.push_back({i, m.matrix[i]});
        }
        count_for_cols[m.matrix[i]] = i;
    }
    for (unsigned i = 1; i <= m.get_cols(); ++i) {
        if (count_for_cols[i]) {
            cols.push_back({i, count_for_cols[i]});
        }
    }
}

Permutation::Permutation(const std::vector<unsigned> &permutation) {
    cols.resize(permutation.size());
    for (int row = permutation.size(); row > 0; --row) {
        rows.push_back({row, permutation[row - 1]});
        cols[permutation[row - 1] - 1] = {permutation[row - 1], row};
    }
}

PermutationMatrix Permutation::expand(unsigned row_amount, unsigned col_amount) const {
    std::vector <unsigned> subpermutation(row_amount, 0);
    for (const auto &permutation_pair: rows) {
        subpermutation[permutation_pair.first - 1] = permutation_pair.second;
    }
    return PermutationMatrix{row_amount, col_amount, subpermutation};
}

// Returns a pair of two permutations, where the first one has
// the (n + 1) / 2 smallest row element values.
std::pair<Permutation, Permutation> Permutation::split_row() const {
    unsigned split_value = rows[rows.size() / 2].first;
    std::vector <std::pair <unsigned, unsigned> > row_first, row_second;
    std::vector <std::pair <unsigned, unsigned> > col_first, col_second;
    for (const auto &matched_pair: rows) {
        if (matched_pair.first <= split_value) {
            row_first.push_back(matched_pair);
        }
        if (matched_pair.first > split_value) {
            row_second.push_back(matched_pair);
        }
    }
    for (const auto &matched_pair: cols) {
        if (matched_pair.second <= split_value) {
            col_first.push_back(matched_pair);
        }
        if (matched_pair.second > split_value) {
            col_second.push_back(matched_pair);
        }
    }
    return {Permutation{row_first, col_first}, Permutation{row_second, col_second}};
}

// Returns a pair of two permutations, where the first one has
// the (n + 1) / 2 smallest column element values.
std::pair<Permutation, Permutation> Permutation::split_col() const {
    unsigned split_value = cols[(cols.size() - 1) / 2].first;
    std::vector <std::pair <unsigned, unsigned> > row_first, row_second;
    std::vector <std::pair <unsigned, unsigned> > col_first, col_second;
    for (const auto &matched_pair: rows) {
        if (matched_pair.second <= split_value) {
            row_first.push_back(matched_pair);
        }
        if (matched_pair.second > split_value) {
            row_second.push_back(matched_pair);
        }
    }
    for (const auto &matched_pair: cols) {
        if (matched_pair.first <= split_value) {
            col_first.push_back(matched_pair);
        }
        if (matched_pair.first > split_value) {
            col_second.push_back(matched_pair);
        }
    }
    return {Permutation{row_first, col_first}, Permutation{row_second, col_second}};
}


// Try to move the ant up.
// It might stop seeing a bad R_high value, or see a new bad R_low value.
// If this happens, the balance is broken and the ant can not move up.
bool SteadyAnt::can_move_up() const {
    PermutationIterator new_high(high_it);
    PermutationIterator new_low(low_it);
    for (; !new_high.has_row_ended() && new_high.row() == ant_row; new_high.inc_row()) {
        if (new_high.matching_col() < ant_col) {
            return false;
        }
    }
    for (; !new_low.has_row_ended() && new_low.row() == ant_row; new_low.inc_row()) {
        if (new_low.matching_col() >= ant_col) {
            return false;
        }
    }
    return !high_it.has_row_ended() || !low_it.has_row_ended();
}

// Try to move the ant right.
// It might stop seeing a bad R_low value, or see a new bad R_high value.
// If this happens, the balance is broken and the ant can not move right.
bool SteadyAnt::can_move_right() const {
    PermutationIterator new_high(high_it);
    PermutationIterator new_low(low_it);
    for (; !new_high.has_col_ended() && new_high.col() == ant_col; new_high.inc_col()) {
        if (new_high.matching_row() <= ant_row) {
            return false;
        }
    }
    for (; !new_low.has_col_ended() && new_low.col() == ant_col; new_low.inc_col()) {
        if (new_low.matching_row() > ant_row) {
            return false;
        }
    }
    return !high_it.has_col_ended() || !low_it.has_col_ended();
}

// Moves the ant up. 
// The up move validity is not checked since it is confirmed in the main traversal.
void SteadyAnt::move_up() {
    for (; !high_it.has_row_ended() && high_it.row() == ant_row; high_it.inc_row()) {
        if (high_it.matching_col() >= ant_col) {
            good_elements_row.push_back(high_it.row_pair());
        }
    }
    for (; !low_it.has_row_ended() && low_it.row() == ant_row; low_it.inc_row()) {
        if (low_it.matching_col() < ant_col) {
            good_elements_row.push_back(low_it.row_pair());
        }
    }
    ant_row = get_next_row();
}

// Moves the ant right. 
// The right move validity is not checked since it is confirmed in the main traversal.
void SteadyAnt::move_right() {
    for (; !high_it.has_col_ended() && high_it.col() == ant_col; high_it.inc_col()) {
        if (high_it.matching_row() > ant_row) {
            good_elements_col.push_back(high_it.col_pair());
        }
    }
    for (; !low_it.has_col_ended() && low_it.col() == ant_col; low_it.inc_col()) {
        if (low_it.matching_row() <= ant_row) {
            good_elements_col.push_back(low_it.col_pair());
        }
    }
    ant_col = get_next_col();
}

unsigned SteadyAnt::get_next_row() const {
    return std::max(low_it.has_row_ended() ? min_row : low_it.row(), 
                    high_it.has_row_ended() ? min_row : high_it.row());
}

unsigned SteadyAnt::get_next_col() const {
    return std::min(low_it.has_col_ended() ? max_col : low_it.col(), 
                    high_it.has_col_ended() ? max_col : high_it.col());
}

SteadyAnt::SteadyAnt(const Permutation &r_low, const Permutation &r_high): 
                                                low_it(PermutationIterator(r_low)),
                                                high_it(PermutationIterator(r_high)) {

    min_row = std::min(r_low.rows.size() ? r_low.rows.back().first : 1, 
                       r_high.rows.size() ? r_high.rows.back().first : 1) - 1;
    max_col = std::max(r_low.cols.size() ? r_low.cols.back().first : 1, 
                       r_high.cols.size() ? r_high.cols.back().first : 1) + 1;
}

Permutation SteadyAnt::restore_correct_product() {
    // The ant position.
    // This is the pair of indexes before which the ant is currently located.
    ant_row = get_next_row();
    ant_col = get_next_col();
    while (ant_row != min_row || ant_col != max_col) {
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
    return Permutation{good_elements_row, good_elements_col};
}

Permutation multiply(const Permutation &p, const Permutation &q) {
    if (p.get_nonzero_amount() == 0 || q.get_nonzero_amount() == 0) {
        // If all elements are zeroes, the product is a zero as well.
        return Permutation({}, {});
    }
    // The recursion base: at most one non-zero in each permutation.
    if (p.get_nonzero_amount() == 1 && q.get_nonzero_amount() == 1) {
        // The only non-zero value in the product C = A * B is 
        // C[i][k] if A[i][_] and B[_][k] are both non-zero.
        std::pair<unsigned, unsigned> p_nonzero = p.rows[0];
        std::pair<unsigned, unsigned> q_nonzero = q.cols[0];
        return Permutation({{p_nonzero.first, q_nonzero.first}},
                           {{q_nonzero.first, p_nonzero.first}});
    }
    // The divide phrase.
    // Split the first matrix by cols and the second by rows on the same it
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

// Recursively multiplies two permutations using the steady ant algorithm.
Permutation Permutation::operator*(const Permutation &p) const {
    return multiply(*this, p);
}

PermutationMatrix PermutationMatrix::operator*(
                        const PermutationMatrix& m) const {
    if (cols != m.get_rows()) {
        throw MatrixException("Sticky Multiplication",
                              "#arg1.columns != #arg2.rows");
    }
    // Preprocessing before feeding everything to the recursive multiply function.
    // Convert the PermutationMatrix to the permutation format used by multiply.
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
