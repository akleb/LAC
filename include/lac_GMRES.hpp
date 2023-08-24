/**
 * @file lac_GMRES.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Generalized Minimum Residual Linear Solver
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_GMRES_HPP_
#define _LAC_GMRES_HPP_

#include "lac_Default.hpp"
#include "lac_Error.hpp"
#include "lac_Structs.hpp"
#include "lac_MatrixMath.hpp"
#include "lac_Norms.hpp"

#include <cstring>
#include <cmath>

/**
 * @brief Applies previous Givens rotations to H_col and comptues the new
 * rotation to remove the last entry in the column
 *
 * @param col the column index we are acting on inside of H
 * @param H_col the col column in the upper Hessenberg without Givens rotations
 * applied
 * @param e the RHS for the least-squares that needs the new Givens rotation
 * applied to it as well
 * @param F stores all the previous Givens rotations that need to be applied to
 * H_col before computing the new rotation, the new rotation is stored in
 * F[2*col] <= c and F[2*col + 1] <= s
 */
int lac_GivensRotation(const int col, double *H_col, double *e, double *F);

/**
 * @brief computes the solution to Ax = b using a matrix free GMRES approach
 *
 * @param linSys a structs that follows the interface provided by 
 * lac_MatrixFreeLinearSystem
 *      int MatVecProd(const double *v, double *Av)
 *      int RightPrecondition(const double *x, double *Mx)
 * @param b the right hand side we are solving for
 * @param n the number of unkowns
 * @param nRst the number of Krylov vectors to build before reseting
 * @param tol the tolerance to converge the residual to relative to initial
 * linear residual
 * @param precondition true to right precondition, false otherwise
 * @param x the solution vector
 * @param size the size of the MPI_COMM_WORLD
 * @param rank the rank of the processor for a parallel solve
 * @param verbose true for extra print out, false for less
 * @param iters set to the number of interior iterations to solve the system
 */
int lac_GMRES(lac_MatrixFreeLinearSystem *linSys, const double *b, const int n, const int nRst, 
               const double tol, const bool precondition, double *x,
               const int size = 0, const int rank = 1, bool verbose = false, 
               int *iters = nullptr);

#endif

