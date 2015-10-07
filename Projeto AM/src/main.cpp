#include <cstdlib>
#include <iostream>

#include "math_lib.h"

int main() {
  std::cout.precision(8);
  std::cout << "Testando o projeto inicialmente." << std::endl;
  math::Matrix A(3, 3), B(2, 1);
  
  for (int i = 0; i < 3; ++i) {
    A(i, i) = 1.5;
  }
  
  A = 3*A;
  std::cout << A.ToString() << std::endl;
  
  A = 2*B;
  std::cout << (A == B) << std::endl;

  std::vector<std::vector<double>> input1(3);
  input1[0] = std::vector<double>(2, 1);
  input1[1] = std::vector<double>({1.4, 0.000000001});
  input1[2] = std::vector<double>({2.5, -2});

  math::Matrix C(input1);
  std::cout << C.ToString() << std::endl;

  B = 2*C;

  try {
    A = B + 0.5*C;
    std::cout << A.ToString() << std::endl;
    std::cout << B.ToString() << std::endl;
  } catch (math::MatrixDimensionMismatchException& e) {
    std::cout << e.what() << std::endl;
  }
  
  try {
    input1[1] = std::vector<double>(0, 10);
    math::Matrix D(input1);
    std::cout << D.ToString() << std::endl;
  } catch (math::BadDimensionException& e) {
    std::cout << e.what() << std::endl;
  } catch (math::MatrixDimensionMismatchException& e) {
    std::cout << e.what() << std::endl;
  }

  std::vector<std::vector<double>> input2(2);
  input2[0] = std::vector<double>({1, 0, -1, 0, 1, -1});
  input2[1] = std::vector<double>({2, 1, 1, -1, 0, 1});

  std::vector<std::vector<double>> input3(6);
  input3[0] = std::vector<double>(1, 1);
  input3[1] = std::vector<double>(1, 0);
  input3[2] = std::vector<double>(1, -2);
  input3[3] = std::vector<double>(1, 1);
  input3[4] = std::vector<double>(1, 3);
  input3[5] = std::vector<double>(1, 2);

  try {
    A = math::Matrix(input2);
    math::Matrix x(input3);
    math::Matrix b;

    std::cout << (A.ApplyGaussianElimination()).ToString() << std::endl;

    b = A*x;
    std::cout << b.ToString() << std::endl;  // b = (4, 1)^T

    b = x*A;  // Exception.
    std::cout << b.ToString() << std::endl;
  } catch (math::MatrixDimensionMismatchException& e) {
    std::cout << e.what() << std::endl;
  }

  std::vector<std::vector<double>> input4(3);
  input4[0] = std::vector<double>({0.0, -1.0, 3.0});
  input4[1] = std::vector<double>({2.0, 1.0, 1.0});
  input4[2] = std::vector<double>({1.0, 0.0, 2.0});

  try {
    math::SquareMatrix M(input4);

    std::cout << M.ToString() << std::endl;
    std::cout << "  Determinante: " << M.GetDeterminant() << std::endl;  // Should be 0.
    std::cout << (M.ApplyGaussianElimination()).ToString() << std::endl;
  } catch (math::MatrixDimensionMismatchException& e) {
    std::cout << e.what() << std::endl;
  } catch (math::BadDimensionException& e) {
    std::cout << e.what() << std::endl;
  }
	
  system("pause");
	return 0;
}