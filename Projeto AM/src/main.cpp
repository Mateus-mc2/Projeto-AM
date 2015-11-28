#include <cstdlib>
#include <iomanip>

#include "fuzzy_algorithm.h"

int main() {
  const int N = 100;
  reader::FileReader dataset_reader;
  project::FuzzyClustering clustering_obj(dataset_reader.delta());
  double arg_min = 0.0;
  math::Matrix fuzzy_partition;

  for (int i = 0; i < N; ++i) {
    math::Matrix curr_partition;
    std::cout << "Getting partition number "<< i + 1 << "..." << std::endl;
    curr_partition = clustering_obj.ExecuteClusteringAlgorithm();

    if (i == 0 || arg_min > clustering_obj.adequacy_criterion()) {
      fuzzy_partition = curr_partition;
      arg_min = clustering_obj.adequacy_criterion();
    }
  }

  std::vector<std::unordered_set<int>> hard_partition = clustering_obj.GetHardPartition(fuzzy_partition);
  
  std::string header = "---------------- Membership matrix for the best partition obtained obtained ----------------";
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

    for (int j = 0; j < hard_partition.size(); ++j) {
      if (hard_partition[j].find(i) != hard_partition[j].end()) {
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
  std::vector<int> prior_medoids = clustering_obj.GetMedoids(dataset_reader.prior_hard_partition());
  std::cout << "Prior partition's medoids:";

  for (int i = 0; i < prior_medoids.size(); ++i) {
    std::cout << " " << prior_medoids[i];
  }
  
  std::cout << std::endl;

  std::vector<int> posterior_medoids = clustering_obj.GetMedoids(hard_partition);
  std::cout << "Posterior partition's medoids:";

  for (int i = 0; i < posterior_medoids.size(); ++i) {
    std::cout << " " << posterior_medoids[i];
  }

  std::cout << std::endl;

  double rand_index = clustering_obj.GetCorrectedRandIndex(dataset_reader.prior_hard_partition(),
                                                           hard_partition);
  std::cout << "---------------\nCorrected Rand Index = " << rand_index << std::endl;

  system("pause");

	return 0;
}