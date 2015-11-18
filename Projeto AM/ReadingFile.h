#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "math_lib.h"

std::vector<double> lineToVectorA(std::string line);
math::Matrix  matrixExamplesA(int totalExamples);
math::Matrix  matrixDissimilarityA(math::Matrix m1, int totalLines);
