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

#include "lac_Default.hpp"
#include "lac_Error.hpp"
#include "lac_Direct.hpp"
#include <cstring>
#include <cmath>

int lac_PLUFactorization(const double *A, const int n, int *P, double *LU){
  for (unsigned long ii = 0; ii < n; ++ii)
    P[ii] = ii;

  // LU is going to be where we store the upper/lower factorization
  // track this so we can delete it
  std::memcpy(LU, A, sizeof(double)*n*n);

  // Loop over the number of equations
  for (unsigned long ii = 0; ii < n; ++ii) {
    // find the max value in current column to be the next row
    double max = 0;
    unsigned long max_i = n;
    for (unsigned long row = ii; row < n; ++row) {
      if (fabs(LU[row*n + ii]) > max) {
        max = fabs(LU[row * n + ii]);
        max_i = row;
      } // if
    }   // for max value
    if (max_i == n) {
      ERR("Inverting singular matrix.\n");
      for (unsigned long row = 0; row < n; ++row) {
        CONT("|");
        for (unsigned long col = 0; col < n; ++col) {
          printf(" %10.3e", A[row * n + col]);
        } // for
        printf("|\n");
      } // for
      return lac_SINGULAR_MATRIX;
    } // if

    // mark the selected permutation in the permutation vector
    double temp = P[ii];
    P[ii] = P[max_i];
    P[max_i] = temp;

    // permute the rows
    for (int jj = 0; jj < n; ++jj){
      temp = LU[ii*n + jj];
      LU[ii*n + jj] = LU[max_i*n + jj];
      LU[max_i*n + jj] = temp;
    } // for

    // Now the matrix is set up to cancel the first column
    for (unsigned long row = ii + 1; row < n; ++row)
      LU[row * n + ii] = LU[row*n + ii] / LU[ii*n + ii];

    // Now we update the block matrix at smaller size with the gaussian
    // elimination
    for (unsigned long row = ii + 1; row < n; ++row) {
      for (unsigned long col = ii + 1; col < n; ++col) {
        LU[row*n + col] = LU[row*n + col] - LU[row*n + ii] * LU[ii*n + col];
      } // for columns
    }   // for rows
  }     // for factorization

  return lac_OK;

} // lac_PLUFactorization

int lac_PLUForwardBackwardSub(const double *LU, const int *P, const int n, 
                              const double *b, double *x){

  // Solving LUx=b
  // Ly = perm*b
  for (unsigned long row = 0; row < n; ++row) {
    x[row] = b[P[row]];
    for (unsigned long col = 0; col < row; ++col)
      x[row] -= x[col] * LU[row*n + col];
  } // for

  // Ux = y
  for (unsigned long row = n - 1; row >= 0; --row) {
    for (unsigned long col = row + 1; col < n; ++col) {
      x[row] -= x[col] * LU[row*n + col];
    } // for columns
    x[row] /= LU[row*n + row];

    if (row == 0)
      break;
  } // for rows

  return lac_OK;

} // lac_PLUForwardBackwardSub

int lac_InvertMatrix(const double *A, const int n, double *A_inv){

  double *LU = new double[n*n];
  int *P = new int[n];

  // Get the PLU factorization
  int ierr = lac_PLUFactorization(A, n, P, LU);
  if (ierr != lac_OK) {delete[] LU; delete[] P; return ierr;}

  // Now we solve for each column in an identity matrix 
  double *b = new double[n];
  std::memset(b, 0, sizeof(double)*n);
  double *x = new double[n];
  for (int ii = 0; ii < n; ++ii){
    if (ii > 0) b[ii - 1] = 0;
    b[ii] = 1;

    // solve the system
    ierr = lac_PLUForwardBackwardSub(LU, P, n, b, x);
    if (ierr != lac_OK) {
      delete[] LU;
      delete[] P;
      delete[] b;
      delete[] x;

      return ierr;

    } // if

    for (int row = 0; row < n; ++row)
      A_inv[row * n + ii] = x[row];

  } // for

  delete[] LU;
  delete[] P;
  delete[] x;
  delete[] b;

  return lac_OK;

} // lac_InvertMatrix
  
int lac_DirectSolve(const double *A, const double *b, const int n, 
                              double *x){
  double *LU = new double[n*n];
  int *P = new int[n];

  // PLU factor the matrix
  int ierr = lac_PLUFactorization(A, n, P, LU);
  if (ierr != lac_OK) {delete[] LU; delete[] P; return ierr; }

  // Solve with the factorization
  ierr = lac_PLUForwardBackwardSub(LU, P, n, b, x);

  delete[] LU;
  delete[] P;

  if (ierr != lac_OK) return ierr;

  return lac_OK;

} // lac_DirectSolve
  
