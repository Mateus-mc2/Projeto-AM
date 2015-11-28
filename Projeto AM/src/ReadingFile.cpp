#include "ReadingFile.h"

namespace reader {

FileReader::FileReader() {
  this->prior_hard_partition_.resize(2);
  math::Matrix data;
  this->ReadFile(&data);
  this->delta_ = this->GetDissimilarityMatrix(data);
}

FileReader::~FileReader() {
  for (int i = 0; i < this->prior_hard_partition_.size(); ++i) {
    this->prior_hard_partition_[i].clear();
  }
}

std::vector<double> FileReader::LineToVector(const std::string &line, int *index) {
  const int kNumTokens = 10;
  std::vector<double> return_vector(kNumTokens - 1);

  const std::string kTokenX("x");
  const std::string kTokenO("o");
  const std::string kTokenB("b");
  const std::string kTokenPositive("positive");
  const std::string kTokenNegative("negative");

  std::istringstream iss(line);

  for (int i = 0; i < kNumTokens; ++i) {
    std::string token;
    iss >> token;

    if (!token.compare(kTokenX)) {
      return_vector[i] = 1;
    } else if (!token.compare(kTokenO)) {
      return_vector[i] = 0;
    } else if (!token.compare(kTokenB)) {
      return_vector[i] = -1;
    } else if (!token.compare(kTokenPositive)) {
      *index = 0;
    } else if (!token.compare(kTokenNegative)) {
      *index = 1;
    }
  }

	return return_vector;
}

//return matrix [num. lines][9]
void FileReader::ReadFile(math::Matrix *data) {
	std::string str_input;
	std::vector<std::vector<double>> princ_matrix;

	std::ifstream inf("../../TestFile.txt");

	if (!inf) {
		std::cout << "//could not be opened for reading!\\" << std::endl;
		exit(1);
	}

  for (int i = 0; getline(inf, str_input); ++i) {
    std::replace(str_input.begin(), str_input.end(), ',', ' ');
    int index;
    std::vector<double> curr_line = this->LineToVector(str_input, &index);

    this->prior_hard_partition_[index].insert(i);
		princ_matrix.push_back(curr_line);
	}

	math::Matrix C(princ_matrix);
  *data = C;
}

math::SquareMatrix FileReader::GetDissimilarityMatrix(const math::Matrix &m1) {
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

}  // namespace reader