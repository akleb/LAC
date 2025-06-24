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

#include <algorithm>
#include <cstring>
#include <cstdio>

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

int lac_BlockCRSTranspose(lac_BlockCRSMatrix *p_Mat, const int block_n, const int block_m){
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
      double *T_data_local = T_data + n_block * cur_index[row];
      const double *data_local = data + n_block * jj;
      for (int block_row = 0; block_row < block_n; ++block_row){
        for (int block_col = 0; block_col < block_m; ++block_col){
          T_data_local[block_n * block_col + block_row] = data_local[block_m * block_row + block_col];
        } // for
      } // for
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

  
//int lac_BlockCRS_ILU0(const lac_BlockCRSMatrix *p_A, const int block_n, lac_BlockCRSMatrix *p_LU){
//  const int n_block = p_A->n_block;
//  const int n = p_A->n;
//  const int *n_col = p_A->n_col;
//  const int *col_index = p_A->col_index;
//
//  std::memcpy(p_LU->data, p_A->data, sizeof(double) * p_A->n_block * p_A->n_nzero);
//  double *a_ik= nullptr;
//  double *a_ij = nullptr;
//  double *a_kj = nullptr;
//  double *a_kk_i = new double[block_n * block_n];
//  double *working = new double[block_n * block_n];
//
//  // For each block row
//  for (int ii = 1; ii < n; ++ii){
//    for (int kk_i = n_col[ii], kk = col_index[kk_i]; kk < ii && kk_i < n_col[ii + 1]; kk = col_index[++kk_i]){
//
//      // skip if (i, k) is zero
//      if (lac_BlockCRSGetData(p_LU, ii, kk, &a_ik) == lac_INDEX_ZERO_ENTRY)
//        return lac_INDEX_ZERO_ENTRY;
//
//      if (lac_BlockCRSGetData(p_LU, kk, kk, &a_ij) == lac_INDEX_ZERO_ENTRY)
//        return lac_INDEX_ZERO_ENTRY;
//       
//      lac_InvertMatrix(a_ij, block_n, a_kk_i);
//
//      // This computes the terms in the lower diagonal portion that gives the
//      // cancel
//      lac_MatRowMatRowMult(a_ik, a_kk_i, block_n, block_n, block_n, working);
//      std::memcpy(a_ik, working, n_block * sizeof(double));
//
//      // Then we apply this transformation to other entries in this row before
//      // doing more cancels
//      for (int jj_i = kk_i + 1; jj_i < n_col[ii + 1]; ++jj_i){
//        const int jj = col_index[jj_i];
//
//        if (lac_BlockCRSGetData(p_LU, kk, jj, &a_kj) == lac_INDEX_ZERO_ENTRY)
//          continue;
//        if (lac_BlockCRSGetData(p_LU, ii, jj, &a_ij) == lac_INDEX_ZERO_ENTRY)
//          return lac_INDEX_ZERO_ENTRY;
//
//        lac_MatRowMatRowMult(a_ik, a_kj, block_n, block_n, block_n, working);
//        for (int entry = 0; entry < block_n * block_n; ++entry)
//            a_ij[entry] -= working[entry];
//      } // for
//    } // for
//  } // for
//
//  delete[] a_kk_i;
//  delete[] working;
//
//  return lac_OK;
//
//} // lac_BlockCRS_ILU0
  
int lac_BlockCRS_ILU0(const lac_BlockCRSMatrix *p_A, const int block_n, lac_BlockCRSMatrix *p_LU){
  const int n = p_A->n;
  const int *n_col = p_A->n_col;
  const int *col_index = p_A->col_index;

  std::memcpy(p_LU->data, p_A->data, sizeof(double) * p_A->n_block * p_A->n_nzero);
  double *a_ik = nullptr;
  double *a_kk= nullptr;
  double *a_kj = nullptr;
  double *a_ij = nullptr;

  // For each block row
  for (int ii = 0; ii < n; ++ii){
    // Do each row in the block infividually
    int kk_i = n_col[ii];
    for (int kk = col_index[kk_i]; kk < ii && kk_i < n_col[ii + 1]; kk = col_index[++kk_i]){
      // This should not be zero because we are iterating on the col_index
      // above.
      LAC_BLOCKCRSGETDATA(p_LU, ii, kk, a_ik);
      // This should also never be zero because it is the diagonal entry
      LAC_BLOCKCRSGETDATA(p_LU, kk, kk, a_kk);
      // Need to iterate through the actual elements in this specific block
      // NEVER HAVE TO WORRY ABOUT THE FIRST ROW BECAUSE THERE WILL BE NO FULL K
      // BLOCKS, THE FIRST BLOCK IN THE ROW WILL BE DIAGONAL
      for (int row = 0; row < block_n; ++row){
        // The index of k in this block
        for (int col = 0; col < block_n; ++col){
          a_ik[block_n * row + col] /= a_kk[block_n * col + col];
          // Need to do the j update for this block while we are doing the k
          // stuff
          for (int jj = col + 1; jj < block_n; ++jj){
            a_ik[block_n * row + jj] -= a_ik[block_n * row + col] * a_kk[block_n * col + jj];
          } // for
        } // for
      } // for

      // Then we apply this transformation to other entries in this row before
      // doing more cancels
      for (int jj_i = kk_i + 1; jj_i < n_col[ii + 1]; ++jj_i){
        const int jj = col_index[jj_i];

        if (lac_BlockCRSGetData(p_LU, kk, jj, &a_kj) == lac_INDEX_ZERO_ENTRY)
          continue;
        // This should always hit because we are iterating over the col_index
        // above
        LAC_BLOCKCRSGETDATA(p_LU, ii, jj, a_ij);
        // operate on each element in the block individually
        for (int row = 0; row < block_n; row++){
          // col_1 is the columns in the ik block each column in this needs to
          // operate on every element in the current row
          for (int col_1 = 0; col_1 < block_n; col_1++){
            // cols  of a_ij
            for (int col = 0; col < block_n; ++col){
              a_ij[block_n * row + col] -= a_ik[block_n * row + col_1] * a_kj[block_n * col_1 + col];
            } // for
          } // for
        } // for
      } // for
    } // for
    
    if (col_index[kk_i] != ii){
      LAC_ERR("Missing diagonal element in BlockCRSMatrix during ILU0 row %d\n", ii);
      return lac_INDEX_ZERO_ENTRY;
    } // if
      
    LAC_BLOCKCRSGETDATA(p_LU, ii, ii, a_kk);
    // Need to do the "diagonal block" specially 
    for (int row = 1; row < block_n; ++row){
      for (int kk = 0; kk < row; ++kk){
        a_kk[block_n * row + kk] /= a_kk[block_n * kk + kk];
        // propagating changes within this block
        for (int jj = kk + 1; jj < block_n; ++jj){
          a_kk[block_n * row + jj] -= a_kk[block_n * row + kk] * a_kk[block_n * kk + jj];
        } // for
      } // for
    } // for
    // Need to propogate the changes we just made down the rest of the rows
    for (int jj_i = kk_i + 1; jj_i < n_col[ii + 1]; ++jj_i){
      int jj = col_index[jj_i];
      LAC_BLOCKCRSGETDATA(p_LU, ii, jj, a_ij);
      a_ik = a_kk;
      for (int row = 1; row < block_n; ++row){
        // This is the index in k in the diagonal block
        for (int col_1 = 0; col_1 < row; ++col_1){
          // This is the index across all columns in j block because they all
          // need to be updated
          for (int col = 0; col < block_n; ++col){
            a_ij[block_n * row + col] -= a_kk[block_n * row + col_1] * a_ij[block_n * col_1 + col];
          } // for
        } // for
      } // for

    } // for
  } // for

  return lac_OK;

} // lac_BlockCRS_ILU0

int lac_BlockCRS_LUForwardBackwardSub(lac_BlockCRSMatrix *p_LU, const int block_n, const double *b, double *x){
  const int n = p_LU->n;
  const int *n_col = p_LU->n_col;
  const int *col_index = p_LU->col_index;
  double *working = new double[block_n];
  double *data;
  int *P = new int[block_n];
  for (int ii = 0; ii < block_n; ++ii)
    P[ii] = ii;

  //Solving LUx = b
  // Ly = b
  std::memcpy(x, b, sizeof(double) * n * block_n);
  for (int row = 0; row < n; ++row){
    for (int col_i = n_col[row], col = col_index[col_i]; col < row && col_i < n_col[row + 1]; col = col_index[++col_i]){
      if (lac_BlockCRSGetData(p_LU, row, col, &data) == lac_INDEX_ZERO_ENTRY)
        return lac_INDEX_ZERO_ENTRY;
      // matrix multiply 
      lac_MatVecMultRow(data, x + block_n * col, block_n, block_n, working);
      for (int ii = 0; ii < block_n; ++ii)
        x[block_n * row + ii] -= working[ii]; 
    } // for
      
    // Should always have something on the diagonal
    LAC_BLOCKCRSGETDATA(p_LU, row, row, data);
    lac_PLForwardSub(data, P, block_n, x + block_n * row, working);
    std::memcpy(x + block_n * row, working, sizeof(double) * block_n);
  } // for

  //Ux = y
  for (int row = n - 1; row >= 0; --row){
    const int start = std::upper_bound(col_index + n_col[row], col_index + n_col[row + 1], row) - col_index; 
    for (int col_i = start; col_i < n_col[row + 1]; ++col_i){
      const int col = col_index[col_i];
      if (lac_BlockCRSGetData(p_LU, row, col, &data) == lac_INDEX_ZERO_ENTRY)
        return lac_INDEX_ZERO_ENTRY;
      lac_MatVecMultRow(data, x + block_n * col, block_n, block_n, working);
      for (int ii = 0; ii < block_n; ++ii)
        x[block_n * row + ii] -= working[ii];
    } // for
    if (lac_BlockCRSGetData(p_LU, row, row, &data) == lac_INDEX_ZERO_ENTRY)
      return lac_INDEX_ZERO_ENTRY;
    lac_UBackwardSub(data, block_n, x + block_n * row);
  } // for

  delete[] working;
  delete[] P;
  return lac_OK;

} // lac_BlockCRSLUForwardBackwardSub
  
int lac_BlockCRSDense(lac_BlockCRSMatrix *p_Mat, const int block_n, const int block_m, double **A){

  if (block_n * block_m != p_Mat->n_block) {
    LAC_ERR("block_n * block_m != p_Mat->n_block, %d * %d != %d", block_n, block_m, p_Mat->n_block);
    return lac_VALUE_ERROR;
  } // if

  *A = new double[p_Mat->n_block * p_Mat->n * p_Mat->m];
  std::memset(*A, 0, sizeof(double) * p_Mat->n_block * p_Mat->n * p_Mat->m);

  const int tot_cols = block_m * p_Mat->m;
  // Loop through non-zero entries in the block matrix
  for (int row = 0; row < p_Mat->n; ++row){
    for (int ii = p_Mat->n_col[row]; ii < p_Mat->n_col[row + 1]; ++ii){
      const int col = p_Mat->col_index[ii];
      // loop through the current non-zero block and add the entries to the
      // dense matrix
      for (int block_row = 0; block_row < block_n; ++block_row){
        for(int block_col = 0; block_col < block_m; ++block_col){
          (*A)[tot_cols * (row * block_n + block_row) + block_m * col + block_col] 
              = p_Mat->data[p_Mat->n_block * ii + block_m * block_row + block_col];
        } // for
      } // for
    } // for
  } // for

  return lac_OK;

} // lac_BlockCRSDense
  
int lac_BlockCRSPrint(lac_BlockCRSMatrix *p_Mat, const int block_n, const int block_m){
  const int n = p_Mat->n;
  const int m = p_Mat->m;

  double **data_pointers = new double*[m];

  for (int row = 0; row < n; ++row){
    // Get all of the data from all the columns
    for (int col = 0; col < m; ++col){
      int ierr = lac_BlockCRSGetData(p_Mat, row, col, data_pointers + col);
      if (ierr == lac_INDEX_ZERO_ENTRY) data_pointers[col] = nullptr;
    } // for

    // Loop over the block rows
    for (int col = 0; col < m; ++ col){
      printf("  ");
      for (int block_col = 0; block_col < block_m; ++block_col)
        printf("---------");
    } // for
    printf("  \n");
    for (int block_row = 0; block_row < block_n; ++block_row){
      for (int col = 0; col < m; ++col){
        printf("| ");
        if (data_pointers[col]){
          for (int block_col = 0; block_col < block_m; ++block_col)
            printf("%8.1e ", data_pointers[col][block_m * block_row + block_col]);
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
  for (int col = 0; col < m; ++ col){
    printf("  ");
    for (int block_col = 0; block_col < block_m; ++block_col)
      printf("---------");
  } // for
  printf("  \n");

  delete[] data_pointers;

  return lac_OK;

} // lac_BlockCRSPrint

int lac_BlockCRSFileWrite(lac_BlockCRSMatrix *p_Mat, const std::string filename){

  FILE *f_out = fopen(filename.c_str(), "w");
  if (f_out == nullptr)
    return lac_FILE_OPEN_ERROR;

  fprintf(f_out, "%d %d %d %d\n", p_Mat->n, p_Mat->m, p_Mat->n_block, p_Mat->n_nzero);

  // first print out the number of entries in each row
  for (int ii = 0; ii < p_Mat->n + 1; ++ii)
    fprintf(f_out, "%d ", p_Mat->n_col[ii]);
  fprintf(f_out, "\n");

  // Next print out the column index followed by data for each entry
  for (int ii = 0; ii < p_Mat->n_nzero; ++ii){
      fprintf(f_out, "%d ", p_Mat->col_index[ii]);
      for (int state = 0; state < p_Mat->n_block; ++state)
        fprintf(f_out, "%25.16e ", p_Mat->data[p_Mat->n_block * ii + state]);
      fprintf(f_out, "\n");
  } // for

  fclose(f_out);

  return lac_OK;

} // lac_BlockCRSFPrint
 
int lac_BlockCRSInitFromFile(lac_BlockCRSMatrix **pp_Mat, const std::string filename){
  *pp_Mat = new lac_BlockCRSMatrix;

  FILE *f_in = fopen(filename.c_str(), "r");
  if (f_in == nullptr)
    return lac_FILE_OPEN_ERROR;

  int read = fscanf(f_in, "%d %d %d %d\n", &(*pp_Mat)->n, &(*pp_Mat)->m, &(*pp_Mat)->n_block, &(*pp_Mat)->n_nzero);
  if (read != 4){
    LAC_ERR("Failed to read the hearder information in %s.\n", filename.c_str());
    return lac_FILE_READ_ERROR;
  } // if 

  (*pp_Mat)->n_col = new int[(*pp_Mat)->n + 1];
  (*pp_Mat)->col_index = new int[(*pp_Mat)->n_nzero];
  (*pp_Mat)->data = new double[(*pp_Mat)->n_block * (*pp_Mat)->n_nzero];

  // first read in the number of entries in each row
  for (int ii = 0; ii < (*pp_Mat)->n + 1; ++ii){
    read = fscanf(f_in, "%d ", (*pp_Mat)->n_col + ii);
    if (read != 1){
      LAC_ERR("Failed to read n_col for index %d in file %s.\n", ii, filename.c_str());
      return lac_FILE_READ_ERROR;
    } // if
  } // for

  // Next read in the column index followed by data for each entry
  for (int ii = 0; ii < (*pp_Mat)->n_nzero; ++ii){
      read = fscanf(f_in, "%d", &(*pp_Mat)->col_index[ii]);
      if (read != 1){
        LAC_ERR("Failed to read col index for entry %d in file %s.\n", ii, filename.c_str());
        return lac_FILE_READ_ERROR;
      } // if
      for (int state = 0; state < (*pp_Mat)->n_block; ++state){
        read = fscanf(f_in, "%lf", (*pp_Mat)->data + (*pp_Mat)->n_block * ii + state);
        if (read != 1){
          LAC_ERR("Failed to read block data at state %d for entry %d in file %s.\n", state, ii, filename.c_str());
          return lac_FILE_READ_ERROR;
        } // if
      } // for
  } // for

  fclose(f_in);

  return lac_OK;

} // lac_BlockCRSFPrint
 
