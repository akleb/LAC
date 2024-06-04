/**
 * @file lac_BlockCRSMatrix.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Block CRS Matrix format
 * @version 0.1
 * @date 2023-12-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_BlockCRSMatrix.hpp"
#include "lac_Direct.hpp"
#include "lac_Default.hpp"
#include "lac_Error.hpp"
#include "lac_MatrixMath.hpp"

#include <cstring>

void lac_BlockCRSInit(lac_BlockCRSMatrix **pp_Mat, const int n_block, const int n,
                      const int m, const int *col_index, const int *n_col){
  
  *pp_Mat = new lac_BlockCRSMatrix;
  (*pp_Mat)->n_nzero = n_col[n];
  (*pp_Mat)->n_block = n_block;
  (*pp_Mat)->n = n;
  (*pp_Mat)->m = m;

  (*pp_Mat)->data = new double[n_block * n_col[n]];
  std::memset((*pp_Mat)->data, 0, sizeof(double) * n_block * n_col[n]);

  (*pp_Mat)->col_index = new int[n_col[n]];
  std::memcpy((*pp_Mat)->col_index, col_index, sizeof(int) * n_col[n]);

  (*pp_Mat)->n_col = new int[n + 1];
  std::memcpy((*pp_Mat)->n_col, n_col, sizeof(int) * (n + 1));

  return;

} // lac_BlockCRSInit

void lac_BlockCRSFree(lac_BlockCRSMatrix **pp_Mat){

  delete[] (*pp_Mat)->data;
  delete[] (*pp_Mat)->n_col;
  delete[] (*pp_Mat)->col_index;
  delete *pp_Mat;
  *pp_Mat = nullptr;

  return;

} // lac_BlockCRSFree

int lac_BlockCRSAddData(lac_BlockCRSMatrix *p_Mat, const int row, const int col, const double *data){
  double *mat_data;
  int err = lac_BlockCRSGetData(p_Mat, row, col, &mat_data);
  if (err != lac_OK) return err;

  for (int ii = 0; ii < p_Mat->n_block; ++ii)
    mat_data[ii] += data[ii];

  return lac_OK;

} // lac_BlockCRSAddData

int lac_BlockCRSSetData(lac_BlockCRSMatrix *p_Mat, const int row, const int col, const double *data){
  double *mat_data;
  int err = lac_BlockCRSGetData(p_Mat, row, col, &mat_data);
  if (err != lac_OK) return err;

  std::memcpy(mat_data, data, sizeof(double)*p_Mat->n_block);

  return lac_OK;

} // lac_BlockCRSSetData

int lac_BlockCRSZeroData(lac_BlockCRSMatrix *p_Mat){

  std::memset(p_Mat->data, 0, sizeof(double) * p_Mat->n_block * p_Mat->n_nzero);
  return lac_OK;

} // lac_BlockCRSZeroData

int lac_BlockCRSGetData(lac_BlockCRSMatrix *p_Mat, const int row, const int col, double **data){
  // Bounds check
  if (row >= p_Mat->n) return lac_INDEX_OUT_OF_BOUNDS;
  if (col >= p_Mat->m) return lac_INDEX_OUT_OF_BOUNDS;

  // need to do a binary search on the column data in specified row
  int *lower_bound = p_Mat->col_index + p_Mat->n_col[row];
  int *upper_bound = p_Mat->col_index + p_Mat->n_col[row + 1];
  int *mid = lower_bound + (upper_bound - lower_bound) / 2;
  while (lower_bound < upper_bound){
    if (col > *mid){
      if (lower_bound == mid) break;
      lower_bound = mid;
    } // if
    else if (col < *mid)
      upper_bound = mid;
    else{
      *data = p_Mat->data + p_Mat->n_block * (mid - p_Mat->col_index);
      return lac_OK;
    } // else
    mid = lower_bound + (upper_bound - lower_bound) / 2;
  } // while

  return lac_INDEX_ZERO_ENTRY;

} // lac_BlockCRSGetData

int lac_BlockCRSTranspose(lac_BlockCRSMatrix *p_Mat){
  const int n_nzero = p_Mat->n_nzero;
  const int n_block = p_Mat->n_block;
  const double *data = p_Mat->data;
  const int *col_index = p_Mat->col_index;
  const int *n_col = p_Mat->n_col;
  const int m = p_Mat->m;
  const int n = p_Mat->n;

  // new data arrays
  double *T_data = new double[n_nzero * n_block];
  int *T_n_col = new int[m + 1];
  std::memset(T_n_col, 0, sizeof(int) * (m + 1));
  int *T_col_index = new int[n_nzero];

  // First let's set the n_col for the new matrix
  for (int ii = 0; ii < n_nzero; ++ii)
    T_n_col[col_index[ii] + 1]++;
  for (int ii = 0; ii < m; ++ii)
    T_n_col[ii + 1] += T_n_col[ii];

  // Now we initialize pointers to all the location in the new col_index data to
  // set which column each piece is in
  int *cur_index = new int[m];
  for (int ii = 0; ii < m; ++ii)
    cur_index[ii] = T_n_col[ii];

  // Iterate over each old row to determine which columns in the new matrix have
  // entries in the associated columns, copy data over too
  for (int ii = 0; ii < n; ++ii){
    for (int jj = n_col[ii]; jj < n_col[ii + 1]; ++jj){
      const int row = col_index[jj];
      T_col_index[cur_index[row]] = ii;
      std::memcpy(T_data + n_block * cur_index[row], data + n_block * jj, sizeof(double) * n_block);
      cur_index[row]++;
    } // for
  } // for
    
  delete[] col_index;
  delete[] cur_index;
  delete[] n_col;
  delete[] data;

  p_Mat->m = n;
  p_Mat->n = m;
  p_Mat->data = T_data;
  p_Mat->n_col = T_n_col;
  p_Mat->col_index = T_col_index;

  return lac_OK;

} // lac_BlockCRSTranspose

int lac_BlockCRSILU0(const lac_BlockCRSMatrix *p_A, const int block_n, lac_BlockCRSMatrix *p_LU){
  const int n_block = p_A->n_block;
  const int n = p_A->n;
  const int *n_col = p_A->n_col;
  const int *col_index = p_A->col_index;

  std::memcpy(p_LU->data, p_A->data, sizeof(double) * p_A->n_block * p_A->n_nzero);
  double *a_ik= nullptr;
  double *a_ij = nullptr;
  double *a_kj = nullptr;
  double *a_kk_i = new double[block_n * block_n];
  double *working = new double[block_n * block_n];

  // For each block row
  for (int ii = 1; ii < n; ++ii){
    for (int kk_i = n_col[ii], kk = col_index[kk_i]; kk < ii; kk = col_index[++kk_i]){

      // skip if (i, k) is zero
      if (lac_BlockCRSGetData(p_LU, ii, kk, &a_ik) == lac_INDEX_ZERO_ENTRY)
        return lac_INDEX_ZERO_ENTRY;

      if (lac_BlockCRSGetData(p_LU, kk, kk, &a_ij) == lac_INDEX_ZERO_ENTRY)
        return lac_INDEX_ZERO_ENTRY;
       
      lac_InvertMatrix(a_ij, block_n, a_kk_i);

      // This computes the terms in the lower diagonal portion that gives the
      // cancel
      lac_MatRowMatRowMult(a_ik, a_kk_i, block_n, block_n, block_n, working);
      std::memcpy(a_ik, working, n_block * sizeof(double));

      // Then we apply this transformation to other entries in this row before
      // doing more cancels
      for (int jj_i = kk_i + 1; jj_i < n_col[ii + 1]; ++jj_i){
        const int jj = col_index[jj_i];

        if (lac_BlockCRSGetData(p_LU, ii, jj, &a_ij) == lac_INDEX_ZERO_ENTRY)
          return lac_INDEX_ZERO_ENTRY;
        if (lac_BlockCRSGetData(p_LU, kk, jj, &a_kj) == lac_INDEX_ZERO_ENTRY)
          continue;

        lac_MatRowMatRowMult(a_ik, a_kj, block_n, block_n, block_n, working);
        for (int entry = 0; entry < block_n * block_n; ++entry)
            a_ij[entry] -= working[entry];
      } // for
    } // for
  } // for

  delete[] a_kk_i;
  delete[] working;

  return lac_OK;

} // lac_BlockCRSILU0
  
int lac_BlockCRSPrint(lac_BlockCRSMatrix *p_Mat, const int block_n, const int block_m){
  const int n = p_Mat->n;
  const int m = p_Mat->m;

  double **data_pointers = new double*[m];

  for (int row = 0; row < n; ++row){
    // Get all of the data from all the rows
    for (int col = 0; col < m; ++col){
      int ierr = lac_BlockCRSGetData(p_Mat, row, col, data_pointers + col);
      if (ierr == lac_INDEX_ZERO_ENTRY) data_pointers[col] = nullptr;
    } // for

    // Loop over the block rows
    printf("  ");
    for (int col = 0; col < m; ++ col){
      printf("  ");
      for (int block_col = 0; block_col < block_m; ++block_col)
        printf("---------");
    } // for
    printf("  \n");
    for (int block_row = 0; block_row < block_n; ++block_row){
      printf("| ");
      for (int col = 0; col < m; ++col){
        printf("| ");
        if (data_pointers[col]){
          for (int block_col = 0; block_col < block_m; ++block_col)
            printf("%8.1e ", data_pointers[col][block_col]);
        } // if
        else{
          for (int block_col = 0; block_col < block_m; ++block_col)
            printf("         ");
        } // else
      } // for
      printf("|\n");
    } // for
  } // for
  // Loop over the block rows
  printf("  ");
  for (int col = 0; col < m; ++ col){
    printf("  ");
    for (int block_col = 0; block_col < block_m; ++block_col)
      printf("---------");
  } // for
  printf("  \n");

  return lac_OK;

} // lac_BlockCRSPrint

