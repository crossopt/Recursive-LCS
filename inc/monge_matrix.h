#ifndef INC_MONGE_MATRIX_H_
#define INC_MONGE_MATRIX_H_

#include <utility>
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

    // Sticky multiplication of two subpermutation matrices.
    // C = A * B if the cross-difference of the product
    // of the Monge matrices corresponding to A and B equals C.
    // Currently NOT the Steady-Ant algorithm. TODO fix this!
    SubpermutationMatrix operator*(const SubpermutationMatrix &m) const;
    unsigned operator() (unsigned x, unsigned y) const override;
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

}  // namespace matrix
}  // namespace LCS

#endif  // INC_MONGE_MATRIX_H_