#include "fuzzy_algorithm.h"

namespace project {

std::vector<std::vector<int>> FuzzyClustering::GenerateRandomPrototypes() {
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  
  std::vector<std::vector<int>> prototypes(this->K);
  
  /*for (int i = 0; i < this->K; ++i) {
    for (int j = 0; j < this->q; ++j) {
      int partition_index = generator() % this->K;
      int element_index = generator() % this->partitions_[partition_index].size();
      prototypes[i].push_back(this->partitions_[partition_index][element_index]);
    }
  }*/

  for (int i = 0; i < this->K; ++i) {
    for (int j = 0; j < this->q; ++j) {
      int element = generator() % this->delta_.rows();
      prototypes[i].push_back(element);
    }
  }

  return prototypes;
}

double FuzzyClustering::GetMembershipDegree(const std::vector<std::vector<int>> &prototypes,
                                            const int &i,
                                            const int &k) {
  double dist_k = 0.0;
  double membership_degree = 0.0;

  for (int j = 0; j < this->q; ++j) {
    dist_k += this->delta_(i, prototypes[k][j]);
  }

  for (int h = 0; h < this->K; ++h) {
    double dist_h = 0.0;

    for (int j = 0; j < this->q; ++j) {
      dist_h += this->delta_(i, prototypes[h][j]);
    }

    membership_degree += pow(dist_k / dist_h, 1 / static_cast<double>(this->m - 1));
  }

  membership_degree = 1 / membership_degree;

  return membership_degree;
}

double FuzzyClustering::GetAdequacyCriterion(const math::Matrix &fuzzy_partition,
                                             const std::vector<std::vector<int>> &prototypes) {
  double criterion = 0.0;

  for (int k = 0; k < fuzzy_partition.cols(); ++k) {
    for (int i = 0; i < fuzzy_partition.rows(); ++i) {
      double dist_k = 0.0;

      for (int j = 0; j < this->q; ++j) {
        dist_k += this->delta_(i, prototypes[k][j]);
      }

      criterion += pow(fuzzy_partition(i, k), this->m)*dist_k;
    }
  }

  return criterion;
}

}  // namespace project