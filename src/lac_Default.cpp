/**
 * @file lac_Default.xpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief The Default functions
 * @version 0.1
 * @date 2025-06-18
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_Error.hpp"
#include <mpi.h>
#include <cstring>

int _lac_GatherFullArray(const double *a, const int n, double **full_a, int *full_n){
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int *each_n = new int[size];
  int *disps = new int[size];
  std::memset(each_n, 0, sizeof(int) * size);
  MPI_Gather(&n, 1, MPI_INT, each_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  disps[0] = 0;
  for (int ii = 1; ii < size; ++ii)
    disps[ii] = disps[ii - 1] + each_n[ii - 1];
  *full_n = disps[size - 1] + each_n[size - 1];
  *full_a = (rank == 0) ? new double[*full_n] : nullptr;
  MPI_Gatherv(a, n, MPI_DOUBLE, *full_a, each_n, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  delete[] each_n;
  delete[] disps;

  return lac_OK;

} // _lac_GatherFullArray

