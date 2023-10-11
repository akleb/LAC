/**
 * @file lac_MatrixMath.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Matrix Math operations in LAC
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_MPI.hpp"
#include "lac_MatrixMath.hpp"
#include "lac_Error.hpp"
#include <cstring>

int lac_MatVecMultCol(const double *A, const double *x, const int n, 
                        const int m, double *b){

  std::memset(b, 0, sizeof(double)*n);
  for (int col = 0; col < m; ++col){
    for (int row = 0; row < n; ++row){
      b[row] += A[n*col + row] * x[col];
    } // for
  } // for

  return lac_OK;

} // lac_MatVecMultCol

int lac_MatVecMultRow(const double *A, const double *x, const int n, 
                        const int m, double *b){
  std::memset(b, 0, sizeof(double)*n);
  for (int row = 0; row < n; ++row){
    for (int col = 0; col < m; ++col){
      b[row] += A[m*row + col] * x[col];
    } // for
  } // for

  return lac_OK;

} // lac_MatVecMultRow

int lac_DotProduct(const double *a, const double *b, const int n, double *dot){
  *dot = 0;
  for (int ii = 0; ii < n; ++ii)
    (*dot) += a[ii]*b[ii];

  return lac_OK;

} // lac_DotProduct

int lac_DotProductAllReduce(const double *a, const double *b, const int n, double *dot){
  double local_dot;
  int ierr = lac_DotProduct(a, b, n, &local_dot);
  if (ierr != lac_OK) return ierr;
  ierr = lac_Allreduce(&local_dot, dot, 1, MPI_SUM);
  if (ierr != lac_OK) return ierr;

  return lac_OK;

} // lac_DotProductAllReduce

