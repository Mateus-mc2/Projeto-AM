// Projeto da disciplina de Aprendizagem de Máquina (IF699), ministrada pelos
// professores Dr. Francisco de Assis Tenório Carvalho (fatc@cin.ufpe.br) e 
// Dr. Cleber Zanchettin (cz@cin.ufpe.br), no período de 2015.2.
// Equipe:
//	- Ana Caroline Ferreira de França (acff@cin.ufpe.br)
//	- Mateus de Freitas Leite (mfl3@cin.ufpe.br)

#include <cstdlib>
#include <iostream>

#include "math_lib.h"

int main() {
	std::cout << "Testando o projeto inicialmente." << std::endl;
  math::Matrix A(3, 3), B(2, 1);
  
  for (int i = 0; i < 3; ++i) {
    A(i, i) = 1.5;
  }
  
  A = 3*A;
  std::cout << A.ToString() << std::endl;
  
  A = 2*B;
  std::cout << (A == B) << std::endl;

  std::vector<std::vector<double>> input(3);
  input[0] = std::vector<double>(2, 1);
  input[1] = std::vector<double>({1.4, 0});
  input[2] = std::vector<double>({2.5, -2});

  math::Matrix C(input);
  std::cout << C.ToString() << std::endl;

  try {
    input[1] = std::vector<double>(0, 10);
    math::Matrix D(input);
    std::cout << D.ToString() << std::endl;
  } catch (math::BadDimensionException& e) {
    std::cout << e.what() << std::endl;
  }
  catch (math::MatrixDimensionMismatchException& e) {
    std::cout << e.what() << std::endl;
  }  
	
  system("pause");
	return 0;
}