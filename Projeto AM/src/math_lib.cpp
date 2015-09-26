#include "math_lib.h"

namespace math {

const char* MathException::what() const {
  return this->kErrorMsg.c_str();
}

Matrix::Matrix(const std::vector<std::vector<double>> &matrix){
  if (matrix.empty()) {
    throw BadDimensionException("Empty matrix.");
  }

  this->rows_ = matrix.size();
  this->cols_ = matrix[0].size();

  for (int i = 0; i < this->rows_; ++i) {
    if (matrix[i].size() != this->cols_) {
      throw MatrixDimensionMismatchException("Input matrix dimensions must be consistent.");
    }
  }

  // Now we initialize our matrix properly.
  this->data_ = new double*[this->rows_];

  for (int i = 0; i < this->rows_; ++i) {
    this->data_[i] = new double[this->cols_];

    for (int j = 0; j < this->cols_; ++j) {
      this->data_[i][j] = matrix[i][j];
    }
    
  }
}

Matrix::Matrix(const int &m, const int &n) : rows_(m), cols_(n) {
  if (m <= 0 || n <= 0){
    throw BadDimensionException("Matrix dimensions must be positive.");
  }

  this->data_ = new double*[m];
  
  for (int i = 0; i < m; ++i) {
    this->data_[i] = new double[n];

    for (int j = 0; j < n; ++j) {
      this->data_[i][j] = 0;
    }
  }
}

Matrix::~Matrix() {
  for (int i = 0; i < this->rows_; ++i) {
    delete[] this->data_[i];
  }

  delete[] this->data_;
}

bool Matrix::operator==(const Matrix &M) const {
  if (this != &M) {
    if (this->rows_ != M.rows() || this->cols_ != M.cols()) {
      return false;
    }

    for (int i = 0; i < this->rows_; ++i) {
      for (int j = 0; j < this->cols_; ++j) {
        if (this->data_[i][j] != M(i, j)) {
          return false;
        }
      }
    }

    return true;
  }
  
  return true;
}

Matrix& Matrix::operator=(const Matrix &M) {
  if (this != &M) {
    if (this->rows_ != M.rows() || this->rows_ != M.cols()) {
      for (int i = 0; i < this->rows_; ++i) {
        delete[] this->data_[i];
      }

      delete[] this->data_;
      
      this->rows_ = M.rows();
      this->cols_ = M.cols();      
      this->data_ = new double*[this->rows_];

      for (int i = 0; i < this->rows_; ++i) {
        this->data_[i] = new double[this->cols_];
      }
    }

    for (int i = 0; i < this->rows_; ++i) {
      for (int j = 0; j < this->cols_; ++j) {
        this->data_[i][j] = M(i, j);
      }
    }
  }

  return *this;
}

Matrix& Matrix::operator+=(const Matrix &M) {
  if (this->rows_ != M.rows() || this->cols_ != M.cols()) {
    throw MatrixDimensionMismatchException("Matrices dimensions do not agree in addition operator.");
  }


}

double& Matrix::At(const uint32_t &i, const uint32_t &j) {
  if (i >= this->rows_ || j >= this->cols_) {
    throw BadIndexException("Index out of bounds.");
  }

  return this->data_[i][j];
}

std::string Matrix::ToString() {
  std::string result = "[";

  for (int i = 0; i < this->rows_; ++i) {
    for (int j = 0; j < this->cols_; ++j) {
      result += std::to_string(this->data_[i][j]) + ", ";
    }

    result.erase(result.size() - 2);
    result += ";\n ";
  }

  result.erase(result.size() - 3);
  result += "]";

  return result;
}

}  // namespace math_lib