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

Matrix::Matrix(const Matrix &M) : rows_(M.rows()), cols_(M.cols()) {
  this->data_ = new double*[this->rows_];

  for (int i = 0; i < this->rows_; ++i) {
    this->data_[i] = new double[this->cols_];

    for (int j = 0; j < this->cols_; ++j) {
      this->data_[i][j] = M(i, j);
    }
  }
}

Matrix::Matrix(const uint32_t &m, const uint32_t &n) : rows_(m), cols_(n) {
  if (m == 0 || n == 0){
    throw BadDimensionException("Matrix dimensions must be positive.");
  }

  this->data_ = new double*[this->rows_];

  for (int i = 0; i < this->rows_; ++i) {
    this->data_[i] = new double[this->cols_];

    for (int j = 0; j < this->cols_; ++j) {
      this->data_[i][j] = 0;
    }
  }
}

Matrix::Matrix() : rows_(this->kDefaultSize), cols_(this->kDefaultSize) {
  this->data_ = new double*[this->rows_];

  for (int i = 0; i < this->rows_; ++i) {
    this->data_[i] = new double[this->cols_];

    for (int j = 0; j < this->cols_; ++j) {
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

void Matrix::CopyFrom(const Matrix &M) {
  if (this->rows_ != M.rows() || this->rows_ != M.cols()) {
    // Must resize the matrix, causing reallocation.
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

Matrix Matrix::add(const Matrix &A, const Matrix &B) {
  if (A.rows() != B.rows() || A.cols() != B.cols()) {
    throw MatrixDimensionMismatchException(
      "Matrices dimensions do not agree in addition operator."
    );
  }

  Matrix result(A.rows(), A.cols());

  for (int i = 0; i < A.rows(); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      result(i, j) = A(i, j) + B(i, j);
    }
  }

  return result;
}

Matrix Matrix::multiply(const Matrix &A, const Matrix &B) {
  if (A.cols() != B.rows()) {
    throw MatrixDimensionMismatchException("Matrices dimensions do not agree in multiplication.");
  }

  Matrix C(A.rows(), B.cols());

  for (int i = 0; i < C.rows(); ++i) {
    for (int j = 0; j < C.cols(); ++j) {
      for (int k = 0; k < A.cols(); ++k) {
        C(i, j) += A(i, k)*B(k, j);
      }
    }
  }

  return C;
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
    this->CopyFrom(M);
  }

  return *this;
}

Matrix& Matrix::operator+=(const Matrix &M) {
  if (this->rows_ != M.rows() || this->cols_ != M.cols()) {
    throw MatrixDimensionMismatchException("Matrices dimensions do not agree in addition operator.");
  }

  for (int i = 0; i < this->rows_; ++i) {
    for (int j = 0; j < this->cols_; ++j) {
      this->data_[i][j] += M(i, j);
    }
  }

  return *this;
}

Matrix& Matrix::operator-=(const Matrix &M) {
  return *this += -M;
}

Matrix& Matrix::operator*=(const Matrix &M) {
  if (this->cols_ != M.rows()) {
    throw MatrixDimensionMismatchException("Matrices dimensions do not agree in multiplication.");
  }

  Matrix result(this->rows_, M.cols());

  for (int i = 0; i < result.rows(); ++i) {
    for (int j = 0; j < result.cols(); ++j) {
      for (int k = 0; k < this->cols_; ++k) {
        result(i, j) = this->data_[i][k] + M(k, j);
      }
    }
  }

  this->CopyFrom(result);
  return *this;
}

Matrix Matrix::operator+(const Matrix &M) {
  Matrix copy(*this);
  return this->add(copy, M);
}

Matrix Matrix::operator*(const Matrix &M) {
  Matrix copy(*this);
  return this->multiply(copy, M);
}

Matrix Matrix::operator-(const Matrix &M) {
  Matrix copy(*this);
  return this->add(copy, -M);
}

Matrix Matrix::operator-() const {
  Matrix copy(*this);

  for (int i = 0; i < copy.rows(); ++i) {
    for (int j = 0; j < copy.cols(); ++j) {
      copy(i, j) = -copy(i, j);
    }
  }

  return copy;
}

double& Matrix::At(const uint32_t &i, const uint32_t &j) {
  if (i >= this->rows_ || j >= this->cols_) {
    throw BadIndexException("Index out of bounds.");
  }

  return this->data_[i][j];
}

std::string Matrix::ToString() const {
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

Matrix operator*(const double &k, const Matrix &M) {
  Matrix result(M.rows(), M.cols());

  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) {
      result(i, j) = k*M(i, j);
    }
  }
  
  return result;  
}

Matrix operator*(const Matrix &M, const double &k) {
  return k*M;
}

}  // namespace math_lib