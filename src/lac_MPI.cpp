/**
 * @file lac_MPI.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief MPI wrappers for LAC
 * @version 0.1
 * @date 2023-08-25
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_MPI.hpp"
#include "lac_Error.hpp"

int lac_Allreduce(const double *send, double *recv, const int count, MPI_Op op){
  
  MPI_Allreduce(send, recv, count, MPI_DOUBLE, op, MPI_COMM_WORLD);

  return lac_OK;

} // lac_Allreduce
