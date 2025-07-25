/**
 * @file test_BlockCRSMatrix.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test the Block CRS Matrix format
 * @version 0.1
 * @date 2023-12-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "unit_test_framework.h"
#include "lac_Error.hpp"
#include "lac_BlockCRSMatrix.hpp"
#include "lac_MatrixMath.hpp"
#include "lac_Direct.hpp"
#include "lac_Norms.hpp"
#include <cstring>

TEST(test_base){

  const int n = 5;
  const int m = 5;
  const int n_nzero = 9;
  const int n_block = 3;
  // MATRIX:
  // | X 0 0 0 X | 
  // | 0 0 0 0 0 |
  // | 0 0 X X X |
  // | 0 X 0 0 0 |
  // | X 0 X 0 X |
  int *n_col = new int[n + 1];
  n_col[0] = 0;
  n_col[1] = 2;
  n_col[2] = 2;
  n_col[3] = 5;
  n_col[4] = 6;
  n_col[5] = 9;

  int *col_index = new int[n_nzero];
  col_index[0] = 0;
  col_index[1] = 4;
  col_index[2] = 2;
  col_index[3] = 3;
  col_index[4] = 4;
  col_index[5] = 1;
  col_index[6] = 0;
  col_index[7] = 2;
  col_index[8] = 4;

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, n_block, n, m, col_index, n_col);

  // ensure that we copied over the appropriate things over
  ASSERT_EQUAL(p_Mat->n, n);
  ASSERT_EQUAL(p_Mat->m, m);
  ASSERT_EQUAL(p_Mat->n_block, n_block);
  ASSERT_EQUAL(p_Mat->n_nzero, n_nzero);

  for (int ii = 0; ii < n + 1; ++ii){
    ASSERT_EQUAL(p_Mat->n_col[ii], n_col[ii]);
  } // for
  for (int ii = 0; ii < n_nzero; ++ii){
    ASSERT_EQUAL(p_Mat->col_index[ii], col_index[ii]);
  } // for
  for (int ii = 0; ii < n_block * n_nzero; ++ii){
    ASSERT_EQUAL(p_Mat->data[ii], 0);
  } // for

  double *data;
  int ierr = lac_BlockCRSGetData(p_Mat, 5, 0, &data);
  ASSERT_EQUAL(ierr, lac_INDEX_OUT_OF_BOUNDS);
  ierr = lac_BlockCRSGetData(p_Mat, 2, 6, &data);
  ASSERT_EQUAL(ierr, lac_INDEX_OUT_OF_BOUNDS);
  ierr = lac_BlockCRSGetData(p_Mat, 1, 3, &data);
  ASSERT_EQUAL(ierr, lac_INDEX_ZERO_ENTRY);
  ierr = lac_BlockCRSGetData(p_Mat, 2, 1, &data);
  ASSERT_EQUAL(ierr, lac_INDEX_ZERO_ENTRY);
  ierr = lac_BlockCRSGetData(p_Mat, 0, 1, &data);
  ASSERT_EQUAL(ierr, lac_INDEX_ZERO_ENTRY);

  ierr = lac_BlockCRSGetData(p_Mat, 3, 1, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  ASSERT_EQUAL(data, p_Mat->data + 5 * n_block);
  ierr = lac_BlockCRSGetData(p_Mat, 2, 2, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  ASSERT_EQUAL(data, p_Mat->data + 2 * n_block);
  ierr = lac_BlockCRSGetData(p_Mat, 2, 3, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  ASSERT_EQUAL(data, p_Mat->data + 3 * n_block);
  ierr = lac_BlockCRSGetData(p_Mat, 2, 4, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  ASSERT_EQUAL(data, p_Mat->data + 4 * n_block);
  ierr = lac_BlockCRSGetData(p_Mat, 0, 0, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  ASSERT_EQUAL(data, p_Mat->data);
  ierr = lac_BlockCRSGetData(p_Mat, 0, 4, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  ASSERT_EQUAL(data, p_Mat->data + 1 * n_block);

  double data_place[3] = {1, 2, 3};
  ierr = lac_BlockCRSAddData(p_Mat, 0, 0, data_place);
  ASSERT_EQUAL(ierr, lac_OK);
  ierr = lac_BlockCRSAddData(p_Mat, 0, 0, data_place);
  ASSERT_EQUAL(ierr, lac_OK);
  ierr = lac_BlockCRSGetData(p_Mat, 0, 0, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  for (int ii = 0; ii < 3; ++ii){
    ASSERT_ALMOST_EQUAL(data[ii], data_place[ii] * 2, 1e-10);
  } // for

  ierr = lac_BlockCRSSetData(p_Mat, 0, 0, data_place);
  ASSERT_EQUAL(ierr, lac_OK);
  ierr = lac_BlockCRSGetData(p_Mat, 0, 0, &data);
  ASSERT_EQUAL(ierr, lac_OK);
  for (int ii = 0; ii < 3; ++ii){
    ASSERT_EQUAL(data[ii], data_place[ii]);
  } // for

  delete[] n_col;
  delete[] col_index;
  lac_BlockCRSFree(&p_Mat);
  ASSERT_EQUAL(p_Mat, nullptr);

  return;

} // test_Init

TEST(test_InitFromFile){

  lac_BlockCRSMatrix *p_Mat;
  int ierr = lac_BlockCRSInitFromFile(&p_Mat, "test_matrix.dat");
  ASSERT_EQUAL(ierr, lac_OK);

  // MATRIX:
  // | X 0 0 0 X 0 0 | 
  // | 0 0 0 0 0 0 X |
  // | 0 0 X X X 0 0 |
  // | 0 X 0 0 0 X 0 |
  // | X 0 X 0 X 0 X |

  ASSERT_EQUAL(p_Mat->n, 5);
  ASSERT_EQUAL(p_Mat->m, 7);
  ASSERT_EQUAL(p_Mat->n_nzero, 12);
  ASSERT_EQUAL(p_Mat->n_block, 2);

  int correct_ncol[6] = {0, 2, 3, 6, 8, 12};
  for (int ii = 0; ii < p_Mat->n + 1; ++ii){
    ASSERT_EQUAL(p_Mat->n_col[ii], correct_ncol[ii]);
  } // for

  int correct_col_index[12] = {0, 4, 6, 2, 3, 4, 1, 5, 0, 2, 4, 6};
  for (int ii = 0; ii < p_Mat->n_nzero; ++ii){
    ASSERT_EQUAL(p_Mat->col_index[ii], correct_col_index[ii]);
  } // for

  for (int row = 0; row < p_Mat->n; ++row){
    for (int ii = p_Mat->n_col[row]; ii < p_Mat->n_col[row + 1]; ++ii){
      double *data;
      ierr = lac_BlockCRSGetData(p_Mat, row, p_Mat->col_index[ii], &data);
      ASSERT_EQUAL(ierr, lac_OK);
      ASSERT_EQUAL(data[0], row);
      ASSERT_EQUAL(data[1], p_Mat->col_index[ii]);
    } // for
  } // for
    //
  lac_BlockCRSFree(&p_Mat);

  return;

} // test_InitFromFile

TEST(test_InitWriteToFile){

  lac_BlockCRSMatrix *p_Mat;
  int ierr = lac_BlockCRSInitFromFile(&p_Mat, "test_matrix.dat");
  ASSERT_EQUAL(ierr, lac_OK);
  ierr = lac_BlockCRSFileWrite(p_Mat, "_test_BlockCRSWrite.dat");
  ASSERT_EQUAL(ierr, lac_OK);
  lac_BlockCRSMatrix *p_Write;
  ierr = lac_BlockCRSInitFromFile(&p_Write, "_test_BlockCRSWrite.dat");
  ASSERT_EQUAL(ierr, lac_OK);


  // MATRIX:
  // | X 0 0 0 X 0 0 | 
  // | 0 0 0 0 0 0 X |
  // | 0 0 X X X 0 0 |
  // | 0 X 0 0 0 X 0 |
  // | X 0 X 0 X 0 X |

  ASSERT_EQUAL(p_Write->n, p_Mat->n);
  ASSERT_EQUAL(p_Write->m, p_Mat->m);
  ASSERT_EQUAL(p_Write->n_nzero, p_Mat->n_nzero);
  ASSERT_EQUAL(p_Write->n_block, p_Mat->n_block);

  for (int ii = 0; ii < p_Write->n + 1; ++ii){
    ASSERT_EQUAL(p_Write->n_col[ii], p_Mat->n_col[ii]);
  } // for

  for (int ii = 0; ii < p_Write->n_nzero; ++ii){
    ASSERT_EQUAL(p_Write->col_index[ii], p_Mat->col_index[ii]);
  } // for

  for (int ii = 0; ii < p_Write->n_nzero * p_Write->n_block; ++ii){
    ASSERT_EQUAL(p_Write->data[ii], p_Mat->data[ii]);
  } // for

  lac_BlockCRSFree(&p_Mat);
  lac_BlockCRSFree(&p_Write);

  return;

} // test_InitFromFile

TEST(test_transpose_blocks){
  const int n = 4;
  const int m = 5;
  const int n_nzero = 6;
  const int n_block = 2;
  const int block_n = 2;
  const int block_m = 1;

  // MATRIX:
  // | X 0 0 0 X | 
  // | 0 0 0 0 0 |
  // | 0 0 X X X |
  // | 0 X 0 0 0 |
  int n_col[5];
  n_col[0] = 0;
  n_col[1] = 2;
  n_col[2] = 2;
  n_col[3] = 5;
  n_col[4] = 6;

  int col_index[6];
  col_index[0] = 0;
  col_index[1] = 4;
  col_index[2] = 2;
  col_index[3] = 3;
  col_index[4] = 4;
  col_index[5] = 1;

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, n_block, n, m, col_index, n_col);
  double data[2];
  for (int row = 0; row < n; ++row){
    for (int jj = n_col[row]; jj < n_col[row+1]; ++jj){
      const int col = col_index[jj];
      data[0] = row;
      data[1] = col;

      int ierr = lac_BlockCRSSetData(p_Mat, row, col, data);
      if (ierr != lac_OK) {ASSERT_TRUE(false);}
    } // for
  } // for

  lac_BlockCRSTranspose(p_Mat, block_n, block_m);
  // MATRIX TRANSPOSED:
  // | X 0 0 0 | 
  // | 0 0 0 X |
  // | 0 0 X 0 |
  // | 0 0 X 0 |
  // | X 0 X 0 |
  ASSERT_EQUAL(p_Mat->n, m);
  ASSERT_EQUAL(p_Mat->m, n);
  ASSERT_EQUAL(p_Mat->n_nzero, 6);
  ASSERT_EQUAL(p_Mat->n_block, 2);

  int correct_n_col[6];
  correct_n_col[0] = 0;
  correct_n_col[1] = 1;
  correct_n_col[2] = 2;
  correct_n_col[3] = 3;
  correct_n_col[4] = 4;
  correct_n_col[5] = 6;
  for (int ii = 0; ii < 6; ++ii){
    ASSERT_EQUAL(p_Mat->n_col[ii], correct_n_col[ii]);
  } // for
    
  int correct_col_index[6];
  correct_col_index[0] = 0;
  correct_col_index[1] = 3;
  correct_col_index[2] = 2;
  correct_col_index[3] = 2;
  correct_col_index[4] = 0;
  correct_col_index[5] = 2;
  for (int ii = 0; ii < 6; ++ii){
    ASSERT_EQUAL(p_Mat->col_index[ii], correct_col_index[ii]);
  } // for

  // check that the data was tranposed
  for (int row = 0; row < n; ++row){
    for (int jj = p_Mat->n_col[row]; jj < p_Mat->n_col[row+1]; ++jj){
      const int col = p_Mat->col_index[jj];

      double *t_data = nullptr;
      int ierr = lac_BlockCRSGetData(p_Mat, row, col, &t_data);
      if (ierr != lac_OK) {ASSERT_TRUE(false);}

      ASSERT_EQUAL(t_data[0], col);
      ASSERT_EQUAL(t_data[1], row);
    } // for
  } // for

  lac_BlockCRSFree(&p_Mat);

  return;

} // test_transpose

TEST(test_transpose_data){
  const int n = 4;
  const int m = 5;
  const int n_nzero = 6;
  const int n_block = 6;
  const int block_n = 2;
  const int block_m = 3;

  // MATRIX:
  // | X 0 0 0 X | 
  // | 0 0 0 0 0 |
  // | 0 0 X X X |
  // | 0 X 0 0 0 |
  int n_col[5];
  n_col[0] = 0;
  n_col[1] = 2;
  n_col[2] = 2;
  n_col[3] = 5;
  n_col[4] = 6;

  int col_index[6];
  col_index[0] = 0;
  col_index[1] = 4;
  col_index[2] = 2;
  col_index[3] = 3;
  col_index[4] = 4;
  col_index[5] = 1;

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, n_block, n, m, col_index, n_col);
  double data[6] = {0, 1, 2, 3, 4, 5};
  for (int row = 0; row < n; ++row){
    for (int jj = n_col[row]; jj < n_col[row+1]; ++jj){
      const int col = col_index[jj];
      int ierr = lac_BlockCRSSetData(p_Mat, row, col, data);
      if (ierr != lac_OK) {ASSERT_TRUE(false);}
    } // for
  } // for

  lac_BlockCRSTranspose(p_Mat, block_n, block_m);
  // MATRIX TRANSPOSED:
  // | X 0 0 0 | 
  // | 0 0 0 X |
  // | 0 0 X 0 |
  // | 0 0 X 0 |
  // | X 0 X 0 |
  ASSERT_EQUAL(p_Mat->n, m);
  ASSERT_EQUAL(p_Mat->m, n);
  ASSERT_EQUAL(p_Mat->n_nzero, 6);
  ASSERT_EQUAL(p_Mat->n_block, 6);

  int correct_n_col[6];
  correct_n_col[0] = 0;
  correct_n_col[1] = 1;
  correct_n_col[2] = 2;
  correct_n_col[3] = 3;
  correct_n_col[4] = 4;
  correct_n_col[5] = 6;
  for (int ii = 0; ii < 6; ++ii){
    ASSERT_EQUAL(p_Mat->n_col[ii], correct_n_col[ii]);
  } // for
    
  int correct_col_index[6];
  correct_col_index[0] = 0;
  correct_col_index[1] = 3;
  correct_col_index[2] = 2;
  correct_col_index[3] = 2;
  correct_col_index[4] = 0;
  correct_col_index[5] = 2;
  for (int ii = 0; ii < 6; ++ii){
    ASSERT_EQUAL(p_Mat->col_index[ii], correct_col_index[ii]);
  } // for

  // check that the data was tranposed
  double data_T[6] = {0, 3, 1, 4, 2, 5};
  for (int row = 0; row < n; ++row){
    for (int jj = p_Mat->n_col[row]; jj < p_Mat->n_col[row+1]; ++jj){
      const int col = p_Mat->col_index[jj];

      double *t_data = nullptr;
      int ierr = lac_BlockCRSGetData(p_Mat, row, col, &t_data);
      if (ierr != lac_OK) {ASSERT_TRUE(false);}

      ASSERT_EQUAL(t_data[0], data_T[0]);
      ASSERT_EQUAL(t_data[1], data_T[1]);
      ASSERT_EQUAL(t_data[2], data_T[2]);
      ASSERT_EQUAL(t_data[3], data_T[3]);
      ASSERT_EQUAL(t_data[4], data_T[4]);
      ASSERT_EQUAL(t_data[5], data_T[5]);
    } // for
  } // for

  lac_BlockCRSFree(&p_Mat);

  return;

} // test_transpose_data
  
TEST(test_ILU0){
  const int n = 4;
  // MATRIX:
  // |  2 -1  0  0 | 
  // | -1  2 -1  0 |
  // |  0 -1  2 -1 |
  // |  0  0 -1  2 |
  
  int n_nzero = 10;
  int col_index[10] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
  int n_col[5] = {0, 2, 5, 8, 10};

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, 4, 4, 4, col_index, n_col);
  double *data;
  for (int row = 0; row < 4; ++row){
    for (int col = 0; col < 4; ++col){
      if (lac_BlockCRSGetData(p_Mat, row, col, &data) == lac_INDEX_ZERO_ENTRY)
        continue;
      std::memset(data, 0, sizeof(double)*4);
      double val = (row == col) ? 2 : -1;
      data[0] = val;
      data[3] = val;
    } // for
  } // for

  double correct[16] = {   2,    -1,     0,    0,
                        -0.5,   1.5,    -1,    0,
                           0, -2./3,  4./3, -1.0,
                           0,     0, -3./4, 5./4};
  lac_BlockCRSMatrix *p_LU;
  lac_BlockCRSInit(&p_LU, 4, 4, 4, col_index, n_col);
  lac_BlockCRS_ILU0(p_Mat, 2, p_LU);

  ASSERT_EQUAL(p_LU->n_nzero, 10);
  ASSERT_EQUAL(p_LU->n_block, 4);
  ASSERT_EQUAL(p_LU->n, 4);
  ASSERT_EQUAL(p_LU->m, 4);
  for (int ii = 0; ii < 5; ++ii){
    ASSERT_EQUAL(p_LU->n_col[ii], n_col[ii]);
  } // for
  for (int ii = 0; ii < 10; ++ii){
    ASSERT_EQUAL(p_LU->col_index[ii], col_index[ii]);
  } // for

  // Check the data matches
  for (int row = 0; row < 4; ++row){
    for (int col_i = n_col[row]; col_i < n_col[row + 1]; ++col_i){
      const int col = col_index[col_i];
      ASSERT_TRUE(lac_BlockCRSGetData(p_LU, row, col, &data) == lac_OK);
      ASSERT_ALMOST_EQUAL(data[0], correct[4 * row + col], 1e-12);
      ASSERT_ALMOST_EQUAL(data[1], 0, 1e-12);
      ASSERT_ALMOST_EQUAL(data[2], 0, 1e-12);
      ASSERT_ALMOST_EQUAL(data[3], correct[4 * row + col], 1e-12);
    } // for
  } // for
  
  //lac_BlockCRSPrint(p_LU, 2, 2);
  double A_col[64] = { 2,  0, -1,  0,  0,  0,  0,  0, 
                       0,  2,  0, -1,  0,  0,  0,  0, 
                      -1,  0,  2,  0, -1,  0,  0,  0,
                       0, -1,  0,  2,  0, -1,  0,  0,
                       0,  0, -1,  0,  2,  0, -1,  0,
                       0,  0,  0, -1,  0,  2,  0, -1,
                       0,  0,  0,  0, -1,  0,  2,  0,
                       0,  0,  0,  0,  0, -1,  0,  2};
  double x[8];
  for (int col = 0; col < 4; ++col){
    ASSERT_EQUAL(lac_BlockCRS_LUForwardBackwardSub(p_LU, 2, A_col + 8 * col, x), lac_OK);
    for (int row = 0; row < 4; ++row){
      if (row == col) {
        ASSERT_ALMOST_EQUAL(x[row], 1, 1e-12);
      } // if
      else{
        ASSERT_ALMOST_EQUAL(x[row], 0, 1e-12);
      } // else
    } // for
  } // for

  lac_BlockCRSFree(&p_Mat);
  lac_BlockCRSFree(&p_LU);
  return;

} // test_ILU0

TEST(test_ILU0_rand){
  const int n = 4;
  // MATRIX:
  // |  2 -1  0  0 | 
  // | -1  2 -1  0 |
  // |  0 -1  2 -1 |
  // |  0  0 -1  2 |
  
  int n_nzero = 10;
  int col_index[10] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
  int n_col[5] = {0, 2, 5, 8, 10};

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, 4, 4, 4, col_index, n_col);
  double *data;
  srand(0);
  for (int row = 0; row < 4; ++row){
    for (int col = 0; col < 4; ++col){
      if (lac_BlockCRSGetData(p_Mat, row, col, &data) == lac_INDEX_ZERO_ENTRY)
        continue;
      double val = (row == col) ? 2 : -1;
      data[0] = val;
      data[3] = val;
      data[1] = ((double) rand() / RAND_MAX) * 0.5;
      data[2] = -data[1];

    } // for
  } // for

  lac_BlockCRSMatrix *p_LU;
  lac_BlockCRSInit(&p_LU, 4, 4, 4, col_index, n_col);
  lac_BlockCRS_ILU0(p_Mat, 2, p_LU);

  ASSERT_EQUAL(p_LU->n_nzero, 10);
  ASSERT_EQUAL(p_LU->n_block, 4);
  ASSERT_EQUAL(p_LU->n, 4);
  ASSERT_EQUAL(p_LU->m, 4);
  for (int ii = 0; ii < 5; ++ii){
    ASSERT_EQUAL(p_LU->n_col[ii], n_col[ii]);
  } // for
  for (int ii = 0; ii < 10; ++ii){
    ASSERT_EQUAL(p_LU->col_index[ii], col_index[ii]);
  } // for

  double x[8];
  for (int ii = 0; ii < 8; ++ii)
    x[ii] = (double) rand() / RAND_MAX - 0.5;

  double Ax[8];
  double AiAx[8];
  lac_MatVecMultBlockCRS(p_Mat, 2, 2, x, Ax);
  lac_BlockCRS_LUForwardBackwardSub(p_LU, 2, Ax, AiAx);
  double R, R_tilde;
  lac_L2Norm(Ax, 8, &R);
  lac_L2Norm(AiAx, 8, &R_tilde);

  //lac_BlockCRSPrint(p_Mat, 2, 2);
  //lac_BlockCRSPrint(p_LU, 2, 2);
  printf("R       = %.3e\nR_tilde = %.3e\n", R, R_tilde);
  ASSERT_TRUE(R_tilde < R);

  lac_BlockCRSFree(&p_Mat);
  lac_BlockCRSFree(&p_LU);
  return;

} // test_ILU0_rand

TEST(test_dense_ILU0){
  const int n = 4;
  // MATRIX:
  // |  2 -1  0  0 | 
  // | -1  2 -1  0 |
  // |  0 -1  2 -1 |
  // |  0  0 -1  2 |
  
  int n_nzero = 16;
  int col_index[16] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  int n_col[5] = {0, 4, 8, 12, 16};

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, 9, 4, 4, col_index, n_col);
  double A[12 * 12];
  std::memset(A, 0, sizeof(double) * 12 * 12);

  double *data;
  srand(0);
  for (int row = 0; row < 4; ++row){
    for (int col = 0; col < 4; ++col){
      if (abs(row - col) > 1) continue;
      LAC_BLOCKCRSGETDATA(p_Mat, row, col, data);

      double val = (row == col) ? 2 : -1;
      data[0] = val;
      data[4] = val;
      data[1] = ((double) rand() / RAND_MAX) * 0.5;
      data[3] = -data[1];
      data[2] = ((double) rand() / RAND_MAX) * 0.5;
      data[6] = -data[2];
      data[5] = ((double) rand() / RAND_MAX) * 0.5;
      data[7] = -data[5];
      data[8] = val;

      A[12 * (3 *row + 0) + 3 * col + 0] = data[0];
      A[12 * (3 *row + 0) + 3 * col + 1] = data[1];
      A[12 * (3 *row + 0) + 3 * col + 2] = data[2];
      A[12 * (3 *row + 1) + 3 * col + 0] = data[3];
      A[12 * (3 *row + 1) + 3 * col + 1] = data[4];
      A[12 * (3 *row + 1) + 3 * col + 2] = data[5];
      A[12 * (3 *row + 2) + 3 * col + 0] = data[6];
      A[12 * (3 *row + 2) + 3 * col + 1] = data[7];
      A[12 * (3 *row + 2) + 3 * col + 2] = data[8];

    } // for
  } // for

  int n_col_nb[13] = {0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144};
  int col_index_nb[12*12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                             0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  lac_BlockCRSMatrix *p_no_block;
  lac_BlockCRSInit(&p_no_block, 1, 12, 12, col_index_nb, n_col_nb);
  std::memcpy(p_no_block->data, A, sizeof(double) * 12 * 12);

  lac_BlockCRSMatrix *p_LU, *p_LU_nb;
  lac_BlockCRSInit(&p_LU, 9, 4, 4, col_index, n_col);
  lac_BlockCRSInit(&p_LU_nb, 1, 12, 12, col_index_nb, n_col_nb);

  lac_BlockCRS_ILU0(p_Mat, 3, p_LU);
  lac_BlockCRS_ILU0(p_no_block, 1, p_LU_nb);
  double A_inv[12 * 12];
  lac_InvertMatrix(A, 12, A_inv);

  //lac_BlockCRSPrint(p_LU, 3, 3);
  //lac_BlockCRSPrint(p_LU_nb, 1, 1);
  
  double col[12];
  for (int ii = 0; ii < 12; ++ii){
    double x[12];
    std::memset(x, 0, sizeof(double) * 12);
    x[ii] = 1;

    lac_BlockCRS_LUForwardBackwardSub(p_LU_nb, 1, x, col);
    for (int check = 0; check < 12; ++check){
      ASSERT_ALMOST_EQUAL(col[check], A_inv[12 * check + ii], 1e-10);
    } // for
  } // for
    
  for (int ii = 0; ii < 12; ++ii){
    double x[12];
    std::memset(x, 0, sizeof(double) * 12);
    x[ii] = 1;

    lac_BlockCRS_LUForwardBackwardSub(p_LU, 3, x, col);
    for (int check = 0; check < 12; ++check){
      ASSERT_ALMOST_EQUAL(col[check], A_inv[12 * check + ii], 1e-10);
    } // for
  } // for

  lac_BlockCRSFree(&p_Mat);
  lac_BlockCRSFree(&p_no_block);
  lac_BlockCRSFree(&p_LU);
  lac_BlockCRSFree(&p_LU_nb);

  return;

} // test_dense_ILU0

TEST(test_ILU0_BlockSize){
  // MATRIX:
  // |  2 -1  0  0 | 
  // | -1  2 -1  0 |
  // |  0 -1  2 -1 |
  // |  0  0 -1  2 |
  
  int n_nzero = 10;
  int col_index[10] = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
  int n_col[5] = {0, 2, 5, 8, 10};

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, 9, 4, 4, col_index, n_col);

  int n_zero = 10 * 9;
  int col_index_nb[10 * 9] = {0, 1, 2, 3, 4, 5, 
                              0, 1, 2, 3, 4, 5, 
                              0, 1, 2, 3, 4, 5, 
                              0, 1, 2, 3, 4, 5, 6, 7, 8,
                              0, 1, 2, 3, 4, 5, 6, 7, 8,
                              0, 1, 2, 3, 4, 5, 6, 7, 8,
                              3, 4, 5, 6, 7, 8, 9, 10, 11,
                              3, 4, 5, 6, 7, 8, 9, 10, 11,
                              3, 4, 5, 6, 7, 8, 9, 10, 11,
                              6, 7, 8, 9, 10, 11,
                              6, 7, 8, 9, 10, 11,
                              6, 7, 8, 9, 10, 11};
  int n_col_nb[13] = {0, 6, 12, 18, 27, 36, 45, 54, 63, 72, 78, 84, 90};
  lac_BlockCRSMatrix * p_Mat_nb;
  lac_BlockCRSInit(&p_Mat_nb, 1, 12, 12, col_index_nb, n_col_nb);


  double *data, *data_nb;
  srand(0);
  for (int row = 0; row < 4; ++row){
    for (int col = 0; col < 4; ++col){
      if (abs(row - col) > 1) continue;
      LAC_BLOCKCRSGETDATA(p_Mat, row, col, data);

      double val = (row == col) ? 2 : -1;
      data[0] = val;
      data[4] = val;
      data[1] = ((double) rand() / RAND_MAX) * 0.5;
      data[3] = -data[1];
      data[2] = ((double) rand() / RAND_MAX) * 0.5;
      data[6] = -data[2];
      data[5] = ((double) rand() / RAND_MAX) * 0.5;
      data[7] = -data[5];
      data[8] = val;

      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 0, 3 * col + 0, data_nb);
      data_nb[0] = data[0];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 0, 3 * col + 1, data_nb);
      data_nb[0] = data[1];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 0, 3 * col + 2, data_nb);
      data_nb[0] = data[2];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 1, 3 * col + 0, data_nb);
      data_nb[0] = data[3];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 1, 3 * col + 1, data_nb);
      data_nb[0] = data[4];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 1, 3 * col + 2, data_nb);
      data_nb[0] = data[5];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 2, 3 * col + 0, data_nb);
      data_nb[0] = data[6];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 2, 3 * col + 1, data_nb);
      data_nb[0] = data[7];
      LAC_BLOCKCRSGETDATA(p_Mat_nb, 3 * row + 2, 3 * col + 2, data_nb);
      data_nb[0] = data[8];
    } // for
  } // for

  double b[12];
  for (int ii = 0; ii < 12; ++ii)
    b[ii] = (double) rand() / RAND_MAX - 0.5;

  lac_BlockCRSMatrix *p_LU;
  lac_BlockCRSInit(&p_LU, 9, 4, 4, col_index, n_col);
  lac_BlockCRS_ILU0(p_Mat, 3, p_LU);
  double x[12];
  lac_BlockCRS_LUForwardBackwardSub(p_LU, 3, b, x);


  lac_BlockCRSMatrix *p_LU_nb;
  lac_BlockCRSInit(&p_LU_nb, 1, 12, 12, col_index_nb, n_col_nb);
  lac_BlockCRS_ILU0(p_Mat_nb, 1, p_LU_nb);
  double x_nb[12];
  lac_BlockCRS_LUForwardBackwardSub(p_LU, 3, b, x_nb);

  for (int ii = 0; ii < 12; ++ii){
    ASSERT_ALMOST_EQUAL(x[ii], x_nb[ii], 1e-10);
  } // for

  lac_BlockCRSFree(&p_Mat);
  lac_BlockCRSFree(&p_Mat_nb);
  lac_BlockCRSFree(&p_LU);
  lac_BlockCRSFree(&p_LU_nb);

  return;

} // test_ILU0_BlockSize
  
TEST(test_Dense){

  // SPARSITY PATTERN
  // | X 0 X |
  // | 0 0 0 |
  // | 0 X 0 |
  int n_col[4] = {0, 2, 2, 3};
  int col_index[3] = {0, 2, 1};
  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, 6, 3, 3, col_index, n_col);
  for (int ii = 0; ii < 18; ++ii)
    p_Mat->data[ii] = ii + 1;
  
  double A_23[54] = { 1,  2,  3,  0,  0,  0,  7,  8,  9,
                      4,  5,  6,  0,  0,  0, 10, 11, 12,
                      0,  0,  0,  0,  0,  0,  0,  0,  0,
                      0,  0,  0,  0,  0,  0,  0,  0,  0,
                      0,  0,  0, 13, 14, 15,  0,  0,  0,
                      0,  0,  0, 16, 17, 18,  0,  0,  0};
  double A_32[54] = { 1,  2,  0,  0,  7,  8,
                      3,  4,  0,  0,  9, 10,
                      5,  6,  0,  0, 11, 12,
                      0,  0,  0,  0,  0,  0, 
                      0,  0,  0,  0,  0,  0, 
                      0,  0,  0,  0,  0,  0, 
                      0,  0, 13, 14,  0,  0, 
                      0,  0, 15, 16,  0,  0, 
                      0,  0, 17, 18,  0,  0};

  double *A;
  lac_BlockCRSDense(p_Mat, 2, 3, &A);
  for (int ii = 0; ii < 54; ++ii){
    ASSERT_EQUAL(A[ii], A_23[ii]);
  } // for
  delete[] A;

  lac_BlockCRSDense(p_Mat, 3, 2, &A);
  for (int ii = 0; ii < 54; ++ii){
    ASSERT_EQUAL(A[ii], A_32[ii]);
  } // for
    
  delete[] A;
  lac_BlockCRSFree(&p_Mat);

  return;

} // test_Dense
  
TEST_MAIN()

