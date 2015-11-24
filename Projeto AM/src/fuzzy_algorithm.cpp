#include "fuzzy_algorithm.h"
#include <Windows.h>

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
    std::unordered_set<int> prototype_set;
    int j = 0;

    while(j < this->q) {
      int element = generator() % this->delta_.rows();

      if (prototype_set.find(element) == prototype_set.end()) {
        prototypes[i].push_back(element);
        prototype_set.insert(element);
        ++j;
      }        
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

std::vector<int> FuzzyClustering::UpdatePrototype(const math::Matrix &fuzzy_partition,
                                                  const int &index) {
  std::vector<int> new_prototype;
  std::unordered_set<int> prototype_set;  // To keep track of elements inserted at constant time.

  for (int j = 0; j < this->q; ++j) {
    int arg_min = 0;
    double min = 0.0;

    for (int h = 0; h < this->delta_.rows(); ++h) {
      if (prototype_set.find(h) == prototype_set.end()) {
        double sum = 0.0;

        for (int i = 0; i < this->delta_.rows(); ++i) {
          sum += pow(fuzzy_partition(i, index), this->m)*this->delta_(i, h);
        }
        
        if (h == 0 || min > sum) {
          min = sum;
          arg_min = h;
        }
      }      
    }

    new_prototype.push_back(arg_min);
    prototype_set.insert(arg_min);
  }

  return new_prototype;
}

math::Matrix FuzzyClustering::ExecuteClusteringAlgorithm() {
  // Initialization.
  std::vector<std::vector<int>> G = this->GenerateRandomPrototypes();
  math::Matrix U(this->delta_.rows(), this->K);

  for (int i = 0; i < U.rows(); ++i) {
    for (int k = 0; k < U.cols(); ++k) {
      U(i, k) = this->GetMembershipDegree(G, i, k);
    }
  }

  double J = this->GetAdequacyCriterion(U, G);

  // Optimization.

  for (int t = 0; t < this->T; ++t) {
    double prev_J = J;

    // Step 1: computation of the best prototypes.

    for (int k = 0; k < G.size(); ++k) {
      G[k] = this->UpdatePrototype(U, k);
    }

    // Step 2: definition of the best fuzzy partition.

    for (int i = 0; i < U.rows(); ++i) {
      for (int k = 0; k < U.cols(); ++k) {
        U(i, k) = this->GetMembershipDegree(G, i, k);
      }
    }
        
    // Analysing topping criterion.

    J = this->GetAdequacyCriterion(U, G);

    if (fabs(J - prev_J) <= this->eps) {
      break;
    }
  }

  return U;
}

std::vector<std::unordered_set<int>> FuzzyClustering::GetHardPartition(
    const math::Matrix &fuzzy_partition) {
  std::vector<std::unordered_set<int>> hard_partition(fuzzy_partition.cols());

  for (int i = 0; i < fuzzy_partition.rows(); ++i) {
    int cluster_index = 0;
    double arg_max = 0.0;

    for (int j = 0; j < fuzzy_partition.cols(); ++j) {
      if (fuzzy_partition(i, j) > arg_max) {
        arg_max = fuzzy_partition(i, j);
        cluster_index = j;
      }
    }

    hard_partition[cluster_index].insert(i);
  }

  return hard_partition;
}

std::vector<int> FuzzyClustering::GetMedoids(
    const std::vector<std::unordered_set<int>> &hard_partition) {
  std::vector<int> medoids(hard_partition.size());

  for (int i = 0; i < hard_partition.size(); ++i) {
    int object = 0;
    double arg_min = 0.0;

    for (std::unordered_set<int>::iterator it = hard_partition[i].begin();
         it != hard_partition[i].end();
         ++it) {
      double sum = 0.0;

      for (int j = 0; j < this->delta_.cols(); ++j) {
        sum += this->delta_((*it), j);
      }

      if (it == hard_partition[i].begin() || sum < arg_min) {
        object = (*it);
        arg_min = sum;
      }
    }

    medoids[i] = object;
  }

  return medoids;
}

}  // namespace project