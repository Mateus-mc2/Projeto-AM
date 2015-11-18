#include "ReadingFile.h"

//x,x,x,x,o,o,x,o,o,positive to [1,1,1,1,0,0,1,0,0]
std::vector<double> lineToVector(const std::string &line) {
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
math::Matrix  matrixExamples() {
	std::string strInput;
	int line = 0;
	std::vector<std::vector<double>> matrixPrinc;

	std::ifstream inf("C:/Users/Suporte/Documents/Visual Studio 2013/Projects/Projeto AM/bin/TestFile.txt");

	if (!inf) {
		std::cout << "//could not be opened for reading!\\" << std::endl;
    system("pause");
		exit(1);
	}

	while (getline(inf, strInput)) {
		matrixPrinc.push_back(lineToVector(strInput));
	}

	math::Matrix C(matrixPrinc);
	return C;
}

math::SquareMatrix matrixDissimilarity(const math::Matrix &m1){
	math::Matrix m2 = m1;  // Auxiliar para comparacao
	math::SquareMatrix result(m1.rows());

  for (int i = 0; i < result.rows(); ++i) {
    for (int l = 0; l < result.cols(); ++l) {
      double sum = 0.0;
      int p = m1.cols();

      for (int j = 0; j < p; ++j)	{
				if (m1(i, j) != m2(l, j)) {
          sum += 1.0;
				}
			}

      result(i, l) = sum;
		}
	}

  // std::cout << result.ToString() << std::endl;

	return result;
}