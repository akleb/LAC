/**
 * @file lac_Direct.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Direct Solver using Gaussian Elimination with partial pivoting to do
 * direct solves
 * @version 0.1
 * @date 2023-07-26
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_DIRECT_HPP_
#define _LAC_DIRECT_HPP_

#include "lac_Default.hpp"
#include "lac_Error.hpp"

/**
 * @brief Computes the PLU factorization given a square matrix A of size nxn
 *
 * @param A the square matrix to factor, size nxn
 * @param n the size of the matrix
 * @param P the permutation matrix stored as a 1D array with indexes into the
 * order of rows
 * @param LU the factorization, lower/upper triangular matrices are stored
 * together in the same nxn matrix because lower triangular will have all 1's on
 * the diagonals, these are omitted from matrix LU and instead it holds the
 * diagonal of the upper matrix
 */
int lac_PLUFactorization(const double *A, const int n, double *P, double *LU);

/**
 * @brief given a PLU factorization and RHS, b, this solves for x (PLU x = b)
 *
 * @param LU the LU factorization in the format described for 
 * lac_PLUFactorization, size nxn 2D matrix
 * @param P the permutation matrix in the format described for 
 * lac_PLUFactorization, size n array
 * @param n defines the size of the system
 * @param b the RHS vector to solve for x, size n array
 * @param x the solution vector we solved for, size n array
 */
int lac_PLUForwardBackwardSub(const double *LU, const int *P, const int n, 
                              const double *b, double *x);

/**
 * @brief computes the inverse of A (Gaussian elimination to get PLU
 * facotrization and then solving RHS's for the identity matrix)
 *
 * @param A an nxn matrix
 * @param n the size of the system
 * @param A_inv set to the inverse of the matrix
 */
int lac_InvertMatrix(const double *A, const int n, double *A_inv);

/**
 * @brief Computes the solution to Ax = b using a direct sole (Gaussian
 * elimination with partial pivoting for numerical stability)
 *
 * @param A nxn matrix
 * @param b the RHS, size n array
 * @param n the size of the system
 * @param x the solution array, size n
 */
int lac_DirectSolve(const double *A, const double *b, const int n, double *x);

#endif

