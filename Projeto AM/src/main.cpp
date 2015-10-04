// Projeto da disciplina de Aprendizagem de Máquina (IF699), ministrada pelos
// professores Dr. Francisco de Assis Tenório Carvalho (fatc@cin.ufpe.br) e 
// Dr. Cleber Zanchettin (cz@cin.ufpe.br), no período de 2015.2.
// Equipe:
//	- Ana Caroline Ferreira de França (acff@cin.ufpe.br)
//	- Mateus de Freitas Leite (mfl3@cin.ufpe.br)

#include <cstdlib>
#include <iostream>
#include "ReadingFile.h"

using namespace std;

int main() {
	matrixDissimilarity(matrixExamples(), 10);
	cout << "Testando o projeto inicialmente." << endl;
	system("pause");

	return 0;
}