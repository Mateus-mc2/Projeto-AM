#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>

#include "fuzzy_algorithm.h"

int main() {
  project::FuzzyClustering clustering_obj;
  std::vector<std::vector<std::unordered_set<int>>> hard_partitions(2);

  for (int k = 0; k < hard_partitions.size(); ++k) {
    math::Matrix fuzzy_partition = clustering_obj.ExecuteClusteringAlgorithm();
    hard_partitions[k] = clustering_obj.GetHardPartition(fuzzy_partition);

    std::string header = "---------------- Membership matrix for the minimal U obtained ----------------";
    std::cout.width(header.size());
    std::cout << header << std::endl;

    std::cout << "Sample\t";

    for (int j = 0; j < fuzzy_partition.cols(); ++j) {
      std::cout << "\tColumn " << j + 1 << "\t";
    }

    std::cout << "\tHard Partition" << std::endl;

    for (int i = 0; i < fuzzy_partition.rows(); ++i) {
      std::cout.width(6);
      std::cout << std::right << i + 1 << "\t\t";
      std::ostringstream str_stream;

      for (int j = 0; j < fuzzy_partition.cols(); ++j) {
        str_stream << "Column " << j + 1 << "\t\t";
        std::string str = str_stream.str();

        str_stream.str(std::string());
        str_stream.clear();

        std::cout.width(str.size());
        std::cout << std::left << fuzzy_partition(i, j);
        std::cout << "\t\t";
      }

      for (int j = 0; j < hard_partitions[k].size(); ++j) {
        if (hard_partitions[k][j].find(i) != hard_partitions[k][j].end()) {
          str_stream << "Hard Partition";
          std::string str = str_stream.str();

          str_stream.str(std::string());
          str_stream.clear();

          std::cout << std::right << j + 1 << std::endl;
          break;
        }
      }
    }

    std::cout << "--- Medoids ---" << std::endl;
    std::vector<int> medoids = clustering_obj.GetMedoids(hard_partitions[k]);

    for (int i = 0; i < medoids.size(); ++i) {
      std::cout << "Partition " << i + 1 << ": " << medoids[i] << std::endl;
    }
  }

  double rand_index = clustering_obj.GetCorrectedRandIndex(hard_partitions[0], hard_partitions[1]);
  std::cout << "Corrected Rand Index = " << rand_index << std::endl;

  system("pause");

	return 0;
}