#ifndef INC_MONGE_MATRIX_H_
#define INC_MONGE_MATRIX_H_

#include <iostream>
#include <utility>
#include <vector>
#include <exception>
#include <string>
#include <algorithm>

namespace LCS {
namespace matrix {

// Exception thrown when errors in matrix computations occur.
class MatrixException: public std::exception {
public:
    MatrixException(const std::string &where, const std::string &what):
                    text(where + ": " + what + ".") {}
    virtual const char* what() const noexcept {
        return text.c_str();
    }
private:
    std::string text;
};

// Base class for implementation details shared by all matrices.
class MatrixInterface {
protected:
    const unsigned rows, cols;
    bool is_element_correct(unsigned x, unsigned y) const {
        return x < rows && y < cols;
    }
public:
    MatrixInterface(): rows(0), cols(0) {}
    MatrixInterface(unsigned rows, unsigned cols): rows(rows), cols(cols) {}

    unsigned get_rows() const {return rows; }
    unsigned get_cols() const {return cols; }

    // Queries the [x][y] element of a matrix.
    // It will be an unsigned integer as all simple subunit-Monge
    // and subpermutation matrix elements are non-negative.
    // Indexing for this operator starts from 0.
    virtual unsigned operator() (unsigned x, unsigned y) const = 0;
};

class MongeMatrix;

// Stores a subpermutation matrix in O(n) space, where n is the
// amount of non-zero elements in the matrix.
// A subpermutation matrix is a a zero-one matrix with
// no more than one non-zero element in every row and column.
class SubpermutationMatrix: public MatrixInterface {
private:
    // Stores a single element for each row: the corresponding non-zero column.
    std::vector <unsigned> matrix;

    friend class Permutation;

public:
    // Basic subpermutation matrix constructor.
    // Stores the passed permutation vector as a subpermutation matrix
    // without performing any explicit correctness checks.
    // The number of rows and columns is explicitly passed as a vector
    // containing the subpermutation details has no information about
    // the amount of columns in its corresponding subpermutation matrix.
    // As half-integer indexing is not practical for coding, normal 0-based
    // indexing is used. However, the permutation is stored in 1-based indexing
    // to allow a 0 to correspond to a non-existing element for the ith row.
    SubpermutationMatrix(unsigned rows,
                         unsigned cols,
                         std::vector <unsigned> permutation):
            MatrixInterface(rows, cols), matrix(std::move(permutation)) {
        matrix.insert(matrix.begin(), 0);
    }

    // Constructs a density matrix for a simple subunit-Monge matrix.
    // Alternatively, constructs a matrix from its distribution matrix.
    // The distribution matrix, aka dominance-sum matrix, is
    // a simple subunit-Monge matrix iff the resulting matrix is
    // a subpermutation matrix.
    // Time complexity is O(rows * cols).
    explicit SubpermutationMatrix(const MongeMatrix &distribution_matrix);

    // Slow sticky multiplication of two subpermutation matrices.
    // C = A * B if the cross-difference of the product
    // of the Monge matrices corresponding to A and B equals C.
    // The time complexity for the multiplication of
    // a n * m matrix by a m * k one is O(n * m * k).
    SubpermutationMatrix operator^(const SubpermutationMatrix &m) const;

    // Sticky multiplication of two subpermutation matrices.
    // O(n log n) using the Steady Ant algorithm, where n is the
    // maximum amount of nonzeroes in the subpermutation matrices.
    SubpermutationMatrix operator*(const SubpermutationMatrix &m) const;

    unsigned operator() (unsigned x, unsigned y) const override;
};


class Permutation {
private:
    std::vector <std::pair <unsigned, unsigned>> row_view;
    std::vector <std::pair <unsigned, unsigned>> col_view;

public:
    explicit Permutation(const std::vector <std::pair <unsigned, unsigned> > &permutation_vector,
                         const std::vector <std::pair <unsigned, unsigned> > &rev_permutation_vector):
                        row_view(permutation_vector),
                        col_view(rev_permutation_vector) {}
    explicit Permutation(const SubpermutationMatrix &m);

    unsigned get_nonzero_amount() const {return row_view.size(); }
    unsigned get_row_split_index() const;
    unsigned get_col_split_index() const;
    std::pair<Permutation, Permutation> split_col(unsigned split_value) const;
    std::pair<Permutation, Permutation> split_row(unsigned split_value) const;
    SubpermutationMatrix expand(unsigned rows, unsigned cols) const;

    std::vector <std::pair<unsigned, unsigned> > get_row_view() const;
    std::vector <std::pair<unsigned, unsigned> > get_col_view() const;
};

// Explicitly stores a simple subunit-Monge matrix.
// A Monge matrix is a matrix the cross-difference matrix
// of which is non-negative.
// A matrix is (sub)unit-Monge if its cross-difference is a
// (sub)permutation matrix.
// A Monge matrix is simple if it is equal to the distribution-sum
// matrix of its cross-difference matrix. This is equivalent to it
// having all zeros in its left column and bottom row.
class MongeMatrix: public MatrixInterface {
private:
    // Stores the matrix of values explicitly.
    std::vector <std::vector<unsigned>> matrix;

public:
    // Basic Monge matrix constructor.
    // Stores the passed matrix as a Monge matrix
    // without performing any explicit correctness checks.
    explicit MongeMatrix(std::vector <std::vector<unsigned>> matrix):
            MatrixInterface(matrix.size(), matrix[0].size()),
            matrix(matrix) {}

    // Explicitly constructs a matrix from its density matrix.
    // The density matrix, aka cross-difference matrix, is
    // a subpermutation matrix iff the resulting matrix is
    // a simple subunit-Monge matrix.
    // Time complexity is O(rows * cols).
    explicit MongeMatrix(const SubpermutationMatrix &density_matrix);

    // Multiplies two Monge matrices in the tropical semiring.
    // Addition is given by the min operator, multiplication by +.
    // As such, the formula for matrix multiplication becomes
    // C[i][k] = min_j(A[i][j] + B[j][k]).
    // The time complexity for the multiplication of
    // a n * m matrix by a m * k one is O(n * m * k).
    MongeMatrix operator*(const MongeMatrix &m) const;
    unsigned operator() (unsigned x, unsigned y) const override;
};


class SteadyAnt {
private:
    std::vector <std::pair <unsigned, unsigned>> low_cols ;  // seen are lower-right
    std::vector <std::pair <unsigned, unsigned>> low_rows ;  // seen are lower-right
    std::vector <std::pair <unsigned, unsigned>> high_cols; // seen are upper-left
    std::vector <std::pair <unsigned, unsigned>> high_rows; // seen are upper-left
    unsigned low_row_position = 0;
    unsigned low_col_position = 0;
    unsigned high_row_position = 0;
    unsigned high_col_position = 0;
    unsigned ant_row, ant_col;

public:
    std::vector <std::pair <unsigned, unsigned>> good_elements_row;
    std::vector <std::pair <unsigned, unsigned>> good_elements_col;

    SteadyAnt(const Permutation &r_low, const Permutation &r_high) {
        low_cols = r_low.get_col_view();  // seen are lower-right
        low_rows = r_low.get_row_view();  // seen are lower-right
        high_cols = r_high.get_col_view(); // seen are upper-left
        high_rows = r_high.get_row_view(); // seen are upper-left

        // The ant starts from high rows (lower-left corner).
        std::reverse(low_rows.begin(), low_rows.end());
        std::reverse(high_rows.begin(), high_rows.end());
    }

    bool have_rows_ended() const {
        return low_row_position == low_rows.size() &&
               high_row_position == high_rows.size();
    }

    bool have_cols_ended() const {
        return low_col_position == low_cols.size() &&
               high_col_position == high_cols.size();
    }

    // Try to move the ant up.
    // It might stop seeing a bad R_high value, or see a new bad R_low value.
    bool can_move_up() const {
        unsigned new_low_position = low_row_position;
        unsigned new_high_position = high_row_position;
        while (new_high_position < high_rows.size() && high_rows[new_high_position].first == ant_row) {
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
    bool can_move_right() const {
        unsigned new_low_position = low_col_position;
        unsigned new_high_position = high_col_position;
        while (new_high_position < high_cols.size() && high_cols[new_high_position].first == ant_col) {
            if (high_cols[new_high_position].second <= ant_row) {
                return false;
            }
            new_high_position++;
        }
        while (new_low_position < low_cols.size() && low_cols[new_low_position].first == ant_col) {
            if (low_cols[new_low_position].second > ant_row) {
                return false;
            }
            new_low_position++;
        }
        return !have_cols_ended();
    }

    void move_up() {
        while (high_row_position < high_rows.size()
            && high_rows[high_row_position].first == ant_row) {
            // check badness before moving on
            if (high_rows[high_row_position].second >= ant_col) {
                good_elements_row.push_back(high_rows[high_row_position]);
            }
            high_row_position++;
        }
        while (low_row_position < low_rows.size()
            && low_rows[low_row_position].first == ant_row) {
            // check badness before moving on
            if (low_rows[low_row_position].second < ant_col) {
                good_elements_row.push_back(low_rows[low_row_position]);
            }
            low_row_position++;
        }
        ant_row = get_next_row();
    }

    void move_right() {
        while (high_col_position < high_cols.size()
            && high_cols[high_col_position].first == ant_col) {
            if (high_cols[high_col_position].second > ant_row) {
                good_elements_col.push_back(high_cols[high_col_position]);
            }
            high_col_position++;
            // check badness before moving on
        }
        while (low_col_position < low_cols.size()
            && low_cols[low_col_position].first == ant_col) {
            // check badness before moving on
            if (low_cols[low_col_position].second <= ant_row) {
                good_elements_col.push_back(low_cols[low_col_position]);
            }
            low_col_position++;
        }
        ant_col = get_next_col();
    }

    unsigned get_next_row() const {
        if (have_rows_ended()) {
            return std::min(low_rows.back().first, high_rows.back().first) - 1;
        } else if (low_row_position == low_rows.size()) {
            return high_rows[high_row_position].first;
        } else if (high_row_position == high_rows.size()) {
            return low_rows[low_row_position].first;
        } else {
            return std::max(low_rows[low_row_position].first,
                            high_rows[high_row_position].first);
        }
    }

    unsigned get_next_col() const {
        if (have_cols_ended()) {
            return std::max(low_cols.back().first, high_cols.back().first) + 1;
        } else if (low_col_position == low_cols.size()) {
            return high_cols[high_col_position].first;
        } else if (high_col_position == high_cols.size()) {
            return low_cols[low_col_position].first;
        } else {
            return std::min(low_cols[low_col_position].first,
                            high_cols[high_col_position].first);
        }
    }

    void do_traversal() {
        if (!low_cols.size()) {
            good_elements_row = high_rows;
            good_elements_col = high_cols;
            return;
        }
        if (!high_cols.size()) {
            good_elements_row = low_rows;
            good_elements_col = low_cols;
            return;
        }

        // The ant is *before* these indexes.
        ant_row = get_next_row();
        ant_col = get_next_col();
        while (!have_rows_ended() || !have_cols_ended()) {
            if (can_move_up()) {
                move_up();
            } else if (can_move_right()) {
                move_right();
            } else {
                good_elements_row.push_back({ant_row, ant_col});
                good_elements_col.push_back({ant_col, ant_row});
                move_up();
                move_right();
            }
        }
        std::reverse(good_elements_row.begin(), good_elements_row.end());
    }

};

}  // namespace matrix
}  // namespace LCS

#endif  // INC_MONGE_MATRIX_H_
