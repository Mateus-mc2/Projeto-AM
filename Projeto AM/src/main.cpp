#include <cstdlib>
#include <iostream>
#include "fuzzy_algorithm.h"

int main() {
  project::FuzzyClustering clustering_obj;
  math::Matrix fuzzy_partition = clustering_obj.ExecuteClusteringAlgorithm();
  std::cout << fuzzy_partition.ToString() << std::endl;

  system("pause");

	return 0;
}