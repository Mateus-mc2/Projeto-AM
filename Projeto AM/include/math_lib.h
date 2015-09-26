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
	int rows_;
	int cols_;
  
	double** data_;
public:
  explicit Matrix(const std::vector<std::vector<double>> &matrix);
  // Creates a null matrix.
  Matrix(const int &m, const int &n);
  ~Matrix();

  bool operator==(const Matrix& M) const;
  Matrix& operator=(const Matrix& M);
  Matrix& operator+=(const Matrix& M);
  Matrix& operator*=(const Matrix& M);

  inline double& operator()(const uint32_t &i, const uint32_t &j) const {
    return this->data_[i][j];
  }

  inline Matrix& operator+(const Matrix& B) {
    *this += B;
    return *this;
  }

  inline friend Matrix& operator*(const double& k, Matrix& M) {
    for (int i = 0; i < M.rows_; ++i) {
      for (int j = 0; j < M.cols_; ++j) {
        M(i, j) *= k;
      }
    }

    return M;
  }

  inline Matrix& operator*(const Matrix& M) {
    *this *= M;
    return *this;
  }

  double& At(const uint32_t &i, const uint32_t &j);
  std::string ToString();

  inline int rows() const { return this->rows_; }
  inline int cols() const { return this->cols_; }
  inline double** data() const { return this->data_; }
};

}  // namespace math_lib

#endif  // PROJETO_AM_INCLUDE_MATH_LIB_H_