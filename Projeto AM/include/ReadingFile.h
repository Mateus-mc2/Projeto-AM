#ifndef READING_FILE_H_
#define READING_FILE_H_

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "math_lib.h"

std::vector<double> lineToVector(const std::string &line);
math::Matrix matrixExamples(const int &total_lines);
math::Matrix matrixDissimilarity(const math::Matrix &m1, const int &total_lines);

#endif  // READING_FILE_H_
