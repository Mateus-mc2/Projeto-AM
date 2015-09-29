#ifndef PROJETO_AM_INCLUDE_MATH_LIB_H_
#define PROJETO_AM_INCLUDE_MATH_LIB_H_

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
  const int kDefaultSize = 10;
  uint32_t rows_;
  uint32_t cols_;
  double** data_;

  void CopyFrom(const Matrix &M);
  Matrix add(const Matrix &A, const Matrix &B);
  Matrix multiply(const Matrix &A, const Matrix &B);
public:
  explicit Matrix(const std::vector<std::vector<double>> &matrix);
  // Copy constructor.
  Matrix(const Matrix &M);
  // Creates a null matrix.
  Matrix(const uint32_t &m, const uint32_t &n);
  // Creates a square null matrix of order n = kDefaultSize.
  Matrix();
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
  std::string ToString() const;

  inline uint32_t rows() const { return this->rows_; }
  inline uint32_t cols() const { return this->cols_; }
  inline double** data() const { return this->data_; }
};

Matrix operator*(const double &k, const Matrix &M);
Matrix operator*(const Matrix &M, const double &k);

}  // namespace math_lib

#endif  // PROJETO_AM_INCLUDE_MATH_LIB_H_