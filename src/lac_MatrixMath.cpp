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

int lac_MatVecMultBlockCRS(const lac_BlockCRSMatrix *A, const int block_n,
                           const int block_m, const double *x, double *b){
  const int n_block = A->n_block;
  const int n = A->n;
  const int *n_col = A->n_col;
  const int *col_index = A->col_index;
  const double *data = A->data;

  std::memset(b, 0, sizeof(double) * n * block_n);
  // loop over the entries in b
  for (int row = 0; row < n; ++row){
    double *b_data = b + block_n * row;
    // loop over the entries in associated row of A
    for (int jj = n_col[row]; jj < n_col[row + 1]; ++jj){
      const double *block_data = data + n_block * jj;
      const double *x_data = x + block_m * col_index[jj];
      // Add the product of these two entries to the associated entries in b
      for (int block_row = 0; block_row < block_n; ++block_row){
        for (int block_col = 0; block_col < block_m; ++block_col){
          b_data[block_row] += block_data[block_m * block_row + block_col] * x_data[block_col];
        } // for
      } // for
    } // for
  } // for

  return lac_OK;

} // lac_MatVecMultBlockCRS


int lac_MatRowMatRowMult(const double *A, const double *B, const int n, 
                                          const int m, const int k, double *C){
  std::memset(C, 0, sizeof(double) * n * m);
  for (int row = 0; row < n; ++row){
    for (int col = 0; col < m; ++col){
      for (int dot = 0; dot < k; ++dot)
        C[m * row + col] += A[k * row + dot] * B[m * dot + col];
    } // for
  } // for

  return lac_OK;

} // lac_MatRowMatRowMult
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

