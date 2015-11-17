#ifndef FUZZY_ALGORITHM_H_
#define FUZZY_ALGORITHM_H_

#include <chrono>
#include <random>

#include "math_lib.h"
#include "ReadingFile.h"

namespace project {

class FuzzyClustering {
  // Parameters for the algorithm.
  const int K = 2;  // Number of partitions.
  const int m = 2;  // Parameter which controls the fuzziness of membership for each object.
  const int T = 100;  // Maximum number of iterations.
  const int q = 2;  // Cardinality of prototypes sets.
  const double eps = 1.0e-10;

  math::SquareMatrix delta_;  // Dissimilarity matrix.
  
  std::vector<std::vector<int>> GenerateRandomPrototypes();
  double GetMembershipDegree(const std::vector<std::vector<int>> &prototypes,
    const int &i, const int &k);
  double GetAdequacyCriterion(const math::Matrix &fuzzy_partition,
    const std::vector<std::vector<int>> &prototypes);
  // Execute procedure to return the prototype which minimizes its respective clustering criterion.
  std::vector<int> UpdatePrototype(const math::Matrix &fuzzy_partition);

  public:
    // TODO(Mateus): criar construtor com os métodos para a leitura da base de dados.
    explicit FuzzyClustering(math::SquareMatrix delta) : delta_(delta) {}
    ~FuzzyClustering() {}

    math::Matrix ExecuteClusteringAlgorithm();
};

}   // namespace project

#endif  // FUZZY_ALGORITHM_H_