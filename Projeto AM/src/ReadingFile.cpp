#include "ReadingFile.h"

//x,x,x,x,o,o,x,o,o,positive to [1,1,1,1,0,0,1,0,0]
std::vector<double> lineToVectorA(const std::string &line) {
	int count = 0;
	std::vector<double> returnVector(9);

	for (int i = 0; i < 18; i = i + 2) {
		if (line[i] == 'x'){
			returnVector[count] = 1;
		} else if (line[i] == 'o') {
			returnVector[count] = 0;
		}	else {
			returnVector[count] = -1;
		}

		++count;
	}

	return returnVector;
}

//return matrix [num. lines][9]
math::Matrix  matrixExamplesA(const int &totalExamples) {
	std::string strInput;
	int line = 0;
	std::vector<std::vector<double>> matrixPrinc(totalExamples);

	std::ifstream inf("TestFile.txt");

	if (!inf) {
		std::cout << "//could not be opened for reading!\\" << std::endl;
		exit(1);
	}

	while (getline(inf, strInput)) {
		matrixPrinc[line] = lineToVectorA(strInput);
		++line;
	}

	math::Matrix C(matrixPrinc);
	return C;
}

math::Matrix  matrixDissimilarityA(const math::Matrix &m1, const int &totalLines){
	int auxLine = 0;
	math::Matrix m2 = m1;//auxiliar para comparacao
	math::Matrix result(totalLines, totalLines);

	while (auxLine < totalLines) {
		for (int i = 0; i < totalLines; i++) {
			for (int j = 0; j < 9; j++)	{
				if (m2(auxLine, j) != m1(i, j)) {
					result(auxLine, i) += 1;
				}
			}
			std::cout << " " << result(auxLine, i);
		}

		std::cout << std::endl;
		++auxLine;
	}

	return result;
}