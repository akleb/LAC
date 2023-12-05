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
  lac_BlockCRSInit(&p_Mat, n_block, n, m, n_nzero, col_index, n_col);

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

TEST_MAIN()

