/**
 * @file lac_Norms.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Computing various matrix/vector norms
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_NORMS_HPP_
#define _LAC_NORMS_HPP_

/**
 * @brief computes the L2 norm of a
 *
 * @param a the vector
 * @param n the size of the vector
 * @param norm the computed L2 norm
 */
int lac_L2Norm(const double *a, const int n, double *norm);

/**
 * @brief computes the L2 norm of a vector across all procs in MPI_COMM_WORLD.
 * This causes a sync across all procs.
 *
 * @param a the vector to compute the L2 norm across
 * @param n the number of elements in the array on this processor
 * @param norm where the final reduced sum is stored
 */
int lac_L2NormAllReduce(const double *a, const int n, double *norm);

/**
 * @brief computes the L2 norm of a vector across all procs in MPI_COMM_WORLD.
 * This causes a sync across all procs. It always sums the data in the same
 * order, so it will match a single proc call of lac_L1Norm. Should only be used
 * for debugging
 *
 * @param a the vector to compute the L2 norm across
 * @param n the number of elements in the array on this processor
 * @param norm where the final reduced sum is stored
 */
int lac_DeterministicL2NormAllReduce(const double *a, const int n, double *norm);

/**
 * @brief computes the L2 norm of a
 *
 * @param a the vector
 * @param n the size of the vector
 * @param norm the computed L1 norm
 */
int lac_L1Norm(const double *a, const int n, double *norm);

/**
 * @brief computes the L1 norm of a vector across all procs in MPI_COMM_WORLD.
 * This causes a sync across all procs.
 *
 * @param a the vector to compute the L1 norm across
 * @param n the number of elements in the array on this processor
 * @param norm where the final reduced sum is stored
 */
int lac_L1NormAllReduce(const double *a, const int n, double *norm);

/**
 * @brief computes the L1 norm of a vector across all procs in MPI_COMM_WORLD.
 * This causes a sync across all procs. It always sums the data in the same
 * order, so it will match a single proc call of lac_L1Norm. Should only be used
 * for debugging
 *
 * @param a the vector to compute the L1 norm across
 * @param n the number of elements in the array on this processor
 * @param norm where the final reduced sum is stored
 */
int lac_DeterministicL1NormAllReduce(const double *a, const int n, double *norm);

#endif

