#ifndef PROJETO_AM_INCLUDE_MATH_LIB_H_
#define PROJETO_AM_INCLUDE_MATH_LIB_H_

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace math {

class MathException : public std::exception {
private:
  const std::string kErrorMsg;
public:
  explicit MathException(const std::string &error) : kErrorMsg(error) {}

  const char* what() const;
};

class BadIndexException : public MathException {
public:
  explicit BadIndexException(const std::string &error) : MathException(error) {}
};

class BadDimensionException : public MathException {
public:
  explicit BadDimensionException(const std::string &error) : MathException(error) {}
};

class MatrixDimensionMismatchException : public MathException {
public:
  explicit MatrixDimensionMismatchException(const std::string &error) : MathException(error) {}
};

class Matrix {
private:
  const uint32_t kDefaultSize = 10;
  uint32_t rows_;
  uint32_t cols_;
  double** data_;

  void CopyFrom(const Matrix &M);
  Matrix Add(const Matrix &A, const Matrix &B);
  Matrix Multiply(const Matrix &A, const Matrix &B);
  void SwapRows(const uint32_t i, const uint32_t j);
protected:
  // Finds the row echelon form (an upper-triangular matrix U) of this matrix.
  void ApplyForwardElimination(Matrix *U);
  // Finds the reduced echelon form of this matrix (it must receive the row echelon form matrix U -
  // see ApplyForwardElimination).
  void ApplyBackSubstitution(Matrix *R);
public:
  explicit Matrix(const std::vector<std::vector<double>> &matrix);
  // Copy constructor.
  Matrix(const Matrix &M);
  // Creates a null matrix.
  Matrix(const uint32_t &m, const uint32_t &n);
  // Creates a null square matrix of order n = kDefaultSize.
  Matrix() : Matrix(this->kDefaultSize, this->kDefaultSize) {}
  ~Matrix();

  bool operator==(const Matrix &M) const;
  Matrix& operator=(const Matrix &M);
  Matrix& operator+=(const Matrix &M);
  Matrix& operator-=(const Matrix &M);
  Matrix& operator*=(const Matrix &M);
  Matrix operator+(const Matrix &M);
  Matrix operator*(const Matrix &M);
  Matrix operator-(const Matrix &M);
  Matrix operator-() const;

  inline double& operator()(const uint32_t &i, const uint32_t &j) const {
    return this->data_[i][j];
  }

  double& At(const uint32_t &i, const uint32_t &j);
  // Returns the reduced row echelon form of this matrix by applying a Gaussian Elimination.
  Matrix ApplyGaussianElimination();
  std::string ToString() const;

  inline uint32_t rows() const { return this->rows_; }
  inline uint32_t cols() const { return this->cols_; }
  inline double** data() const { return this->data_; }
};

Matrix operator*(const double &k, const Matrix &M);
Matrix operator*(const Matrix &M, const double &k);

class ColumnMatrix : public Matrix {
private:
  const int kDefaultDimension = 10;

  std::vector<std::vector<double>> ToColumnMatrix(const std::vector<double> &vector);
public:
  // Given a standard vector (from std), constructs a n x 1 matrix, where n is 
  // the argument vector's size.
  explicit ColumnMatrix(const std::vector<double> &vector) : Matrix(this->ToColumnMatrix(vector)) {}
  // Copy constructor.
  ColumnMatrix(const ColumnMatrix &v) : Matrix(v) {}
  // Creates a null vector.
  ColumnMatrix(const uint32_t &dimension) : Matrix(dimension, 1) {}
  // Creates a null vector with default dimension n = kDefaultDimension.
  ColumnMatrix() : ColumnMatrix(this->kDefaultDimension) {}
  ~ColumnMatrix() {}
};

class SquareMatrix : public Matrix {
private:
public:
  explicit SquareMatrix(const std::vector<std::vector<double>> &matrix);
  // Copy constructor.
  SquareMatrix(const SquareMatrix &M) : Matrix(M) {}
  // Creates a null square matrix with a given order.
  SquareMatrix(const uint32_t &order) : Matrix(order, order) {}
  // Creates a null square matrix of order n = kDefaultSize (see Matrix class definition).
  SquareMatrix() : Matrix() {}
  ~SquareMatrix() {}

  // Returns this square matrix determinant.
  double GetDeterminant();
};

}  // namespace math_lib

#endif  // PROJETO_AM_INCLUDE_MATH_LIB_H_