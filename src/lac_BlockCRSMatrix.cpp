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
#include "lac_Default.hpp"
#include "lac_Error.hpp"

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


