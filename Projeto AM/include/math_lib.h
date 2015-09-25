#ifndef PROJETO_AM_INCLUDE_MATH_LIB_H_
#define PROJETO_AM_INCLUDE_MATH_LIB_H_

#include <string>
#include <vector>

namespace math {

class Matrix {
public:
  Matrix(const std::vector<std::vector<double>> &matrix, const int &m, const int &n);
  // Creates a null matrix.
  Matrix(const int &m, const int &n);
  ~Matrix();

  bool operator==(const Matrix& to_compare) const;
  Matrix& operator=(const Matrix& other) const;
  Matrix& operator=(Matrix&& other) const;
  Matrix& operator+=(const Matrix& matrix) const;
  Matrix& operator*=(const Matrix& matrix) const;

  friend Matrix operator+(const Matrix& A, const Matrix& B) {
    A += B;
    return A;
  }

  friend Matrix operator*(const int& k, Matrix& M) {
    for (int i = 0; i < M.rows(); ++i) {
      for (int j = 0; j < M.cols(); ++j) {
        M.At(i, j) *= k;
      }
    }

    return M;
  }

  friend Matrix operator*(Matrix& A, const Matrix& B) {
    A *= B;
    return A;
  }

  double& At(int i, int j);
  std::string ToString();

  inline int rows() const { return this->kRows; }
  inline int cols() const { return this->kCols; }
  inline double** matrix() const { return this->matrix_; }
private:
  const int kRows;
  const int kCols;

  void CopyFrom(const Matrix& matrix);

  double** matrix_;
};



}  // namespace math_lib

#endif  // PROJETO_AM_INCLUDE_MATH_LIB_H_