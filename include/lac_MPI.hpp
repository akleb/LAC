/**
 * @file lac_MPI.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief MPI wrappers for LAC
 * @version 0.1
 * @date 2023-08-25
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_MPI_HPP_
#define _LAC_MPI_HPP_

#include <mpi.h>

/**
 * @brief does a global reduce on all procs of a double across MPI_COMM_WORLD
 *
 * @param send the send buffer
 * @param recv the recv buffer
 * @param count the number of items we are reducing
 * @param op the operation to perform for the reduction
 */
int lac_Allreduce(const double *send, double *recv, const int count, MPI_Op op);

#endif
