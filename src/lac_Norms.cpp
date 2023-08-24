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

#include "lac_Error.hpp"
#include "lac_Norms.hpp"
#include "mpi.h"
#include <cmath>

int lac_L2Norm(const double *a, const int n, double *norm){
  *norm = 0;
  for (int ii = 0; ii < n; ++ii)
    (*norm) += a[ii]*a[ii];
  (*norm) = sqrt(*norm);

  return lac_OK;

} // lac_L2Norm

int lac_L2NormAllReduce(const double *a, const int n, double *norm){
  *norm = 0;
  double local_norm = 0;
  for (int ii = 0; ii < n; ++ii)
    local_norm += a[ii]*a[ii];

  MPI_Allreduce(&local_norm, norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  (*norm) = sqrt(*norm);

  return lac_OK;

} // lac_L2NormAllReduce
 
int lac_L1Norm(const double *a, const int n, double *norm){
  *norm = 0;
  for (int ii = 0; ii < n; ++ii)
    (*norm) += fabs(a[ii]);

  return lac_OK;

} // lac_L1Norm

