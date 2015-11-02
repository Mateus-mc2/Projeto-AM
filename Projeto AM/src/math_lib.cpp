#include <iostream>
#include <cassert>

#include "math_lib.h"

namespace math {

const char* MathException::what() const {
  return this->kErrorMsg.c_str();
}

Matrix::Matrix(const std::vector<std::vector<double>> &matrix) : precision_(this->GetPrecision()) {
  if (matrix.empty() || matrix[0].empty()) {
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
      if (fabs(matrix[i][j]) < kEps) {
        this->data_[i][j] = 0.0;
      }
      else {
        this->data_[i][j] = matrix[i][j];
      }
    }
    
  }
}

Matrix::Matrix(const Matrix &M) : rows_(M.rows()), cols_(M.cols()), precision_(this->GetPrecision()) {
  this->data_ = new double*[this->rows_];

  for (int i = 0; i < this->rows_; ++i) {
    this->data_[i] = new double[this->cols_];

    for (int j = 0; j < this->cols_; ++j) {
      this->data_[i][j] = M(i, j);
    }
  }
}

Matrix::Matrix(const uint32_t &m, const uint32_t &n) : rows_(m), cols_(n), precision_(this->GetPrecision()) {
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

Matrix::~Matrix() {
  for (int i = 0; i < this->rows_; ++i) {
    delete[] this->data_[i];
  }

  delete[] this->data_;
}

int Matrix::GetPrecision() {
  double eps = kEps;
  int precision = 0;

  while (eps != 1.0) {
    ++precision;
    eps *= 10;
  }

  return precision;
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

Matrix Matrix::Add(const Matrix &A, const Matrix &B) {
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

Matrix Matrix::Multiply(const Matrix &A, const Matrix &B) {
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

void Matrix::SwapRows(const uint32_t i, const uint32_t j) {
  double* aux = this->data_[i];
  this->data_[i] = this->data_[j];
  this->data_[j] = aux;
}

void Matrix::ApplyForwardElimination(Matrix *U, int *num_permutations) {
  int last_pivot_col = 0;
  int min = (U->rows() <= U->cols()) ? U->rows() : U->cols();
  bool precision_flag = false;

  for (int i = 0; i < min; ++i) {
    bool found_next_pivot = false;

    for (int j = last_pivot_col; j < U->cols() && !found_next_pivot; ++j) {
      if (!(*U)(i, j)) {
        for (int k = i + 1; k < U->rows() && !found_next_pivot; ++k) {
          if ((*U)(k, j)) {
            U->SwapRows(i, k);
            last_pivot_col = j;
            found_next_pivot = true;
            ++(*num_permutations);
          }
        }
      } else {
        last_pivot_col = j;
        found_next_pivot = true;
      }
    }

    if (found_next_pivot) {
      // It also says that (*U)(i, last_pivot_col) != 0.
      assert((*U)(i, last_pivot_col));

      for (int k = i + 1; k < U->rows(); ++k) {
        double scalar_factor = (*U)(k, last_pivot_col) / (*U)(i, last_pivot_col);

        for (int j = last_pivot_col; j < U->cols(); ++j) {
          int abs = ((*U)(k, j) - scalar_factor*(*U)(i, j) >= 0.0) ?
                    (*U)(k, j) - scalar_factor*(*U)(i, j) :
                    -(*U)(k, j) + scalar_factor*(*U)(i, j);

          if (abs > 0.0 && abs < kEps) {
            (*U)(k, j) = 0.0;
            precision_flag = true;
          } else {
            (*U)(k, j) -= scalar_factor*(*U)(i, j);
          }
        }
      }
    } else {
      break;
    }
  }

  if (precision_flag) {
    std::cout << "  *** WARNING *** " << std::endl << "Aqui vai o aviso pros desavisados: escalonamento impreciso, hohoho. :D" << std::endl;
  }
}

void Matrix::ApplyBackSubstitution(Matrix *R) {
  int last_pivot_col = R->cols();

  for (int i = R->rows() - 1; i >= 0; --i) {
    // Find next pivot.
    bool found = false;

    for (int j = 0; j < last_pivot_col; ++j) {
      if ((*R)(i, j)) {
        last_pivot_col = j;
        found = true;
      }
    }

    if (found) {
      // Set all elements from the pivot's column to zero (in fact, if i is the pivot's line, 
      // it is enough to set all elements (*R)(k, j), k < i, to zero, since R is already at row
      // echelon form).
      for (int k = i - 1; k >= 0; --k) {
        // To perform it by subtracting rows, just do
        // (Line k) <- (Line k) - ((*R)(k, last_pivot_col)/(*R)(i, last_pivot_col)*(Line i),
        // thus (*R)(k, last_pivot_col) <- (*R)(k, last_pivot_col) -
        // - ((*R)(k, last_pivot_col)/(*R)(i, last_pivot_col)) - (*R)(i, last_pivot_col) =
        // = (*R)(k, last_pivot_col) - (*R)(k, last_pivot_col) = 0. The result must be the same.
        (*R)(k, last_pivot_col) = 0.0;
      }
    }
  }
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
  return this->Add(copy, M);
}

Matrix Matrix::operator*(const Matrix &M) {
  Matrix copy(*this);
  return this->Multiply(copy, M);
}

Matrix Matrix::operator-(const Matrix &M) {
  Matrix copy(*this);
  return this->Add(copy, -M);
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

Matrix Matrix::ApplyGaussianElimination() {
  Matrix result(*this);
  int num_permutations = 0;  // Not needed here.

  this->ApplyForwardElimination(&result, &num_permutations);
  this->ApplyBackSubstitution(&result);

  return result;
}

std::string Matrix::ToString() const {
  std::ostringstream out;
  out.precision(this->precision_);
  std::string result = "[";

  for (int i = 0; i < this->rows_; ++i) {
    for (int j = 0; j < this->cols_; ++j) {
      out << this->data_[i][j];
      result += out.str() + ", ";
      out.str("");
      out.clear();
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

SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> &matrix) : Matrix(matrix) {
  if (matrix.size() != matrix[0].size()) {
    throw MatrixDimensionMismatchException("Input matrix must be a square matrix.");
  }
}

double SquareMatrix::GetDeterminant() {
  Matrix U(*this);
  int num_permutations = 0;
  this->ApplyForwardElimination(&U, &num_permutations);
  double diag_product = 1.0;

  for (int i = 0; i < U.rows(); ++i) {
    diag_product *= U(i, i);
  }

  return pow(-1, num_permutations)*diag_product;
}

}  // namespace math_lib