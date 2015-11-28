#ifndef READING_FILE_H_
#define READING_FILE_H_

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "math_lib.h"

namespace reader {

class FileReader {
  private:
    std::vector<std::unordered_set<int>> prior_hard_partition_;
    math::SquareMatrix delta_;

    std::vector<double> LineToVector(const std::string &line, int *index);
    void ReadFile(math::Matrix *delta);
    math::SquareMatrix GetDissimilarityMatrix(const math::Matrix &m1);
  public:
    FileReader();
    ~FileReader();

    // Accessors.
    inline math::SquareMatrix delta() const { return this->delta_; }
    inline std::vector<std::unordered_set<int>> prior_hard_partition() const {
      return this->prior_hard_partition_; 
    }
};


}  // namespace reader


#endif  // READING_FILE_H_
