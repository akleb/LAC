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
  const int n_block = 12;
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
  ASSERT_EQUAL(p_Mat->n_block, 12);

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

TEST_MAIN()

