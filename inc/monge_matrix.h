#ifndef INC_MONGE_MATRIX_H_
#define INC_MONGE_MATRIX_H_

#include <vector>
#include <exception>
#include <string>

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
    unsigned rows, cols;
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
class PermutationMatrix: public MatrixInterface {
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
    PermutationMatrix(unsigned rows,
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
    explicit PermutationMatrix(const MongeMatrix &distribution_matrix);

    // Slow sticky multiplication of two subpermutation matrices.
    // C = A * B if the cross-difference of the product
    // of the Monge matrices corresponding to A and B equals C.
    // The time complexity for the multiplication of
    // a n * m matrix by a m * k one is O(n * m * k).
    PermutationMatrix operator^(const PermutationMatrix &m) const;

    // Sticky multiplication of two subpermutation matrices.
    // O(n log n) using the Steady Ant algorithm, where n is the
    // maximum amount of nonzeroes in the subpermutation matrices.
    PermutationMatrix operator*(const PermutationMatrix &m) const;

    unsigned operator() (unsigned x, unsigned y) const override;
};

// Utility class for storing a permutation as a list of pairs. 
class Permutation {
public:
    explicit Permutation(const std::vector <std::pair <unsigned, unsigned> > &permutation_vector,
                         const std::vector <std::pair <unsigned, unsigned> > &rev_permutation_vector):
                        rows(permutation_vector),
                        cols(rev_permutation_vector) {}
    explicit Permutation(const PermutationMatrix &m);
    explicit Permutation() {}

    // Returns the amount of non-zero elements in the permutation.
    unsigned get_nonzero_amount() const {return rows.size(); }
    // Splits the permutation into two roughly equal halves by their column values.
    std::pair<Permutation, Permutation> split_col() const;
    // Splits the permutation into two roughly equal halves by their row values.
    std::pair<Permutation, Permutation> split_row() const;
    // Expands the permutation into a PermutationMatrix
    // with the necessary amount of rows and columns.
    PermutationMatrix expand(unsigned row_amount, unsigned col_amount) const;

    // The list of permutation pairs, sorted by row index in descending order.
    std::vector <std::pair <unsigned, unsigned>> rows;
    // The list of permutation pairs, sorted by col index in ascending order.
    std::vector <std::pair <unsigned, unsigned>> cols;

    // Sticky multiplication of two permutation.
    // O(n log n) using the Steady Ant algorithm, where n is the
    // amount of elements in the permutation vectors.
    Permutation operator*(const Permutation &p) const;

    // Adds an Id matrix to the beginning of the permutation
    // so the new amount of rows in it is new_rows.
    // All existing elements' indexes are incremented.
    void grow_front(unsigned new_rows);

    // Adds an Id matrix to the beginning of the permutation
    // so the new amount of cols in it is new_cols.
    // All existing elements' indexes are left unchanged.
    void grow_back(unsigned new_cols);
};

// Utility class for iterating over permutations in the Steady Ant algorithm.
class PermutationIterator {
private:
    const Permutation &permutation;
    unsigned row_it;
    unsigned col_it;
public:
    explicit PermutationIterator(const Permutation &permutation): permutation(permutation),
                                                                  row_it(0),
                                                                  col_it(0) {}

    // Checks whether the iterator has passed the last row with permutation elements.
    bool has_row_ended() const {
        return row_it == permutation.rows.size();
    }

    // Checks whether the iterator has passed the last col with permutation elements.
    bool has_col_ended() const {
        return col_it == permutation.cols.size();
    }

    // Increments the row iterator.
    void inc_row() {
        row_it++;
    }
    // Increments the col iterator.
    void inc_col() {
        col_it++;
    }

    // Returns the element the row iterator is pointing to.
    std::pair <unsigned, unsigned> row_pair() const {
        return permutation.rows[row_it];
    }
    // Returns the element the col iterator is pointing to.
    std::pair <unsigned, unsigned> col_pair() const {
        return permutation.cols[col_it];
    }

    // Returns the row of the element the row iterator is pointing to.
    unsigned row() const {
        return row_pair().first;
    }
    // Returns the col of the element the col iterator is pointing to.
    unsigned col() const {
        return col_pair().first;
    }

    // Returns the matching row of the element the col iterator is pointing to.
    unsigned matching_row() const {
        return col_pair().second;
    }
    // Returns the matching col of the element the row iterator is pointing to.
    unsigned matching_col() const {
        return row_pair().second;
    }
};

// Utility class encapsulating the main logic of the ant traversal for two permutations.
class SteadyAnt {
private:
    std::vector <std::pair <unsigned, unsigned>> good_elements_row;
    std::vector <std::pair <unsigned, unsigned>> good_elements_col;

    PermutationIterator low_it;  // ant-visible elements are lower-right
    PermutationIterator high_it;  // ant-visible elements are upper-left

    unsigned ant_row, ant_col;
    unsigned min_row, max_col;

    // Returns true if the ant can be moved up to the next row with permutation elements.
    bool can_move_up() const;
    // Returns true if the ant can be moved right to the next col with permutation elements.
    bool can_move_right() const;

    // Moves the ant up to the next row with permutation elements.
    void move_up();
    // Moves the ant right to the next col with permutation elements.
    void move_right();

    // Returns the next interesting row (with permutation elements) for the ant,
    // or a fixed row smaller than all existing ones if all such rows have been visited.
    unsigned get_next_row() const;
    // Returns the next interesting column (with permutation elements) for the ant,
    // or a fixed col smaller than all existing ones if all such cols have been visited.
    unsigned get_next_col() const;
public:
    // Initializes the steady ant traversal for a pair of r_low, r_high permutation matrices.
    SteadyAnt(const Permutation &r_low, const Permutation &r_high);

    // Does the main ant traversal and returns the fixed permutation product.
    Permutation restore_correct_product();
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
    explicit MongeMatrix(const PermutationMatrix &density_matrix);

    // Multiplies two Monge matrices in the tropical semiring.
    // Addition is given by the min operator, multiplication by +.
    // As such, the formula for matrix multiplication becomes
    // C[i][k] = min_j(A[i][j] + B[j][k]).
    // The time complexity for the multiplication of
    // a n * m matrix by a m * k one is O(n * m * k).
    MongeMatrix operator*(const MongeMatrix &m) const;
    unsigned operator() (unsigned x, unsigned y) const override;
};

}  // namespace matrix
}  // namespace LCS

#endif  // INC_MONGE_MATRIX_H_
