#ifndef FUZZY_ALGORITHM_H_
#define FUZZY_ALGORITHM_H_

#include <chrono>
#include <random>

#include "math_lib.h"
#include "ReadingFile.h"

namespace project {

class FuzzyClustering {
  private:
    // Parameters for the algorithm.
    const int K = 2;  // Number of partitions.
    const int m = 2;  // Parameter which controls the fuzziness of membership for each object.
    const int T = 150;  // Maximum number of iterations.
    const int q = 2;  // Cardinality of prototypes sets.
    const double eps = 1.0e-10;

    math::SquareMatrix delta_;
    math::Matrix fuzzy_partition_;
    double adequacy_criterion_;

    int BinomialCoefficient(const int &n, const int &k) const;
    std::vector<std::vector<int>> GenerateRandomPrototypes();
    double GetMembershipDegree(const std::vector<std::vector<int>> &prototypes,
      const int &i, const int &k);
    double GetAdequacyCriterion(const math::Matrix &fuzzy_partition,
      const std::vector<std::vector<int>> &prototypes);
    // Execute procedure to return the prototype which minimizes its respective clustering criterion.
    std::vector<int> UpdatePrototype(const math::Matrix &fuzzy_partition, const int &index);
  
  public:
    explicit FuzzyClustering(const math::SquareMatrix &delta) : delta_(delta) {}
    ~FuzzyClustering() {}

    // Returns the matrix U which corresponds to the fuzzy partition of this data set. 
    math::Matrix ExecuteClusteringAlgorithm();
    std::vector<std::unordered_set<int>> GetHardPartition(const math::Matrix &fuzzy_partition);
    std::vector<int> GetMedoids(const std::vector<std::unordered_set<int>> &hard_partition);
    double GetCorrectedRandIndex(
      const std::vector<std::unordered_set<int>> &prior_hard_partition,
      const std::vector<std::unordered_set<int>> &posterior_hard_partition
    );

    // Accessors.
    inline math::Matrix fuzzy_partition() const { return this->fuzzy_partition_; }
    inline math::SquareMatrix dataset_reader() const { return this->delta_; }
    inline double adequacy_criterion() const { return this->adequacy_criterion_; }
};

}   // namespace project

#endif  // FUZZY_ALGORITHM_H_