/**
 * @file lac_Norms.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Functions for computing matrix/vector norms
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_MPI.hpp"
#include "lac_Error.hpp"
#include "lac_Norms.hpp"
#include "mpi.h"
#include <cmath>
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

int lac_L2Norm(const double *a, const int n, double *norm){
  *norm = 0;
  for (int ii = 0; ii < n; ++ii)
    (*norm) += a[ii]*a[ii];
  (*norm) = sqrt(*norm);

  return lac_OK;

} // lac_L2Norm

int lac_L2NormAllReduce(const double *a, const int n, double *norm){

  double local_norm = 0;
  for (int ii = 0; ii < n; ++ii)
    local_norm += a[ii]*a[ii];

  *norm = 0;
  lac_Allreduce(&local_norm, norm, 1, MPI_SUM);
  
  (*norm) = sqrt(*norm);

  return lac_OK;

} // lac_L2NormAllReduce
 
int lac_DeterministicL2NormAllReduce(const double *a, const int n, double *norm){
  int full_n;
  double *full_a;
  int ierr = _lac_GatherFullArray(a, n, &full_a, &full_n);
  if (ierr != lac_OK) return ierr;

  ierr = lac_L2Norm(full_a, full_n, norm);
  if (ierr != lac_OK) return ierr;

  MPI_Bcast(norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  delete[] full_a;
  return lac_OK;

} // lac_L2NormAllReduce

int lac_L1Norm(const double *a, const int n, double *norm){
  *norm = 0;
  for (int ii = 0; ii < n; ++ii)
    (*norm) += fabs(a[ii]);

  return lac_OK;

} // lac_L1Norm

int lac_L1NormAllReduce(const double *a, const int n, double *norm){
  double local_norm;
  lac_L1Norm(a, n, &local_norm);

  lac_Allreduce(&local_norm, norm, 1, MPI_SUM);

  return lac_OK;

} // lac_L1NormAllReduce
  
int lac_DeterministicL1NormAllReduce(const double *a, const int n, double *norm){
  int full_n;
  double *full_a;
  int ierr = _lac_GatherFullArray(a, n, &full_a, &full_n);
  if (ierr != lac_OK) return ierr;

  ierr = lac_L1Norm(full_a, full_n, norm);
  if (ierr != lac_OK) return ierr;

  MPI_Bcast(norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  delete[] full_a;
  return lac_OK;

} // lac_L2NormAllReduce

