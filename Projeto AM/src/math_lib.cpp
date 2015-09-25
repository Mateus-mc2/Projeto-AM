#include "math_lib.h"

namespace math {

Matrix::Matrix(const std::vector<std::vector<double>> &matrix, const int &m, const int &n)
  : kRows(m), 
    kCols(n) {
  if (matrix.size() != m) {
    // throw MatrixDimensionMismatchException("Input matrix dimensions must be consistent.");
  }

  for (int i = 0; i < m; ++i) {
    if (matrix[i].size() != n) {
      // throw MatrixDimensionMismatchException("Input matrix dimensions must be consistent.");
    }
  }

  // Now we initialize our matrix properly.
  this->matrix_ = new double*[m];

  for (int i = 0; i < m; ++i) {
    this->matrix_[i] = const_cast<double*>(matrix[i].data());
  }
}

Matrix::Matrix(const int &m, const int &n) : kRows(m), kCols(n) {
  this->matrix_ = new double*[m];
  
  for (int i = 0; i < m; ++i) {
    this->matrix_[i] = new double[n];

    for (int j = 0; j < n; ++j) {
      this->matrix_[i][j] = 0;
    }
  }
}

Matrix::~Matrix() {
  for (int i = 0; i < this->kRows; ++i) {
    delete[] this->matrix_[i];
  }

  delete[] this->matrix_;
}

}  // namespace math_lib