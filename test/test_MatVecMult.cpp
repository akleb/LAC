/**
 * @file test_MatVecMult.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test the Matrix vector multiply functions
 * @version 0.1
 * @date 2023-10-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "unit_test_framework.h"
#include "lac_MatrixMath.hpp"
#include "lac_BlockCRSMatrix.hpp"
#include "lac_Error.hpp"

TEST(test_identity){

  double I[16] = {1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1};
  double b[4];
  srand(0);
  for (int ii = 0; ii < 4; ++ii)
    b[ii] = (double) rand() / RAND_MAX;

  double calc_b[4];

  int ierr = lac_MatVecMultCol(I, b, 4, 4, calc_b);
  ASSERT_EQUAL(ierr, lac_OK);

  for (int ii = 0; ii < 4; ++ii){
    ASSERT_ALMOST_EQUAL(calc_b[ii], b[ii], 1e-10);
  } // for
  
  ierr = lac_MatVecMultRow(I, b, 4, 4, calc_b);
  ASSERT_EQUAL(ierr, lac_OK);

  for (int ii = 0; ii < 4; ++ii){
    ASSERT_ALMOST_EQUAL(calc_b[ii], b[ii], 1e-10);
  } // for

  return;

} // test_identity
  
TEST(test_rand){

  const double A_row[6] = {0, 0.5, 2,
                          -1, 2.4, 1};
  const double A_col[6] = {0, -1,
                           0.5, 2.4, 
                           2, 1};

  const double x[3] = {2.4, 0.2, -1};
  const double b[2] = {-1.9, -2.92};
  double calc_b[2];

  int ierr = lac_MatVecMultRow(A_row, x, 2, 3, calc_b);
  ASSERT_EQUAL(ierr, lac_OK);

  for (int ii = 0; ii < 2; ++ii){
    ASSERT_ALMOST_EQUAL(calc_b[ii], b[ii], 1e-10);
  } // for
  
  ierr = lac_MatVecMultCol(A_col, x, 2, 3, calc_b);
  ASSERT_EQUAL(ierr, lac_OK);

  for (int ii = 0; ii < 2; ++ii){
    ASSERT_ALMOST_EQUAL(calc_b[ii], b[ii], 1e-10);
  } // for

  return;

} // test_identity

TEST(test_BlockCRS){
  // MATRIX
  // | 2 3 0 0 0 0 2 2 |
  // | 0 2 0 0 0 0 1 1 |
  // | 1 0 0 0 0 0 7 9 |
  // | 0 0 9 5 0 0 0 0 |
  // | 0 0 7 0 0 0 0 0 |
  // | 0 0 0 2 0 0 0 0 |
  // | 0 0 0 0 0 0 8 8 |
  // | 0 0 0 0 0 0 0 2 |
  // | 0 0 0 0 0 0 0 9 |
  // | X 0 0 X |
  // | 0 X 0 0 |
  // | 0 0 0 X |
  int n = 3;
  int m = 4;
  int n_col[4] = {0, 2, 3, 4};
  int col_index[4] = {0, 3, 1, 3};

  lac_BlockCRSMatrix *p_Mat;
  lac_BlockCRSInit(&p_Mat, 6, n, m, col_index, n_col);

  double data[6] = {2, 3, 0, 2, 1, 0};
  lac_BlockCRSSetData(p_Mat, 0, 0, data);

  data[0] = 2;
  data[1] = 2;
  data[2] = 1;
  data[3] = 1;
  data[4] = 7;
  data[5] = 9;
  lac_BlockCRSSetData(p_Mat, 0, 3, data);

  data[0] = 9;
  data[1] = 5;
  data[2] = 7;
  data[3] = 0;
  data[4] = 0;
  data[5] = 2;
  lac_BlockCRSSetData(p_Mat, 1, 1, data);

  data[0] = 8;
  data[1] = 8;
  data[2] = 0;
  data[3] = 2;
  data[4] = 0;
  data[5] = 9;
  lac_BlockCRSSetData(p_Mat, 2, 3, data);

  double x[8] = {2, 3, 3, 4, 8, 2, 4, 7};
  double b[9];
  double b_correct[9] = {35, 17, 93, 47, 21, 8, 88, 14, 63};

  lac_MatVecMultBlockCRS(p_Mat, 3, 2, x, b);

  for (int ii = 0; ii < 9; ++ii){
    ASSERT_ALMOST_EQUAL(b[ii], b_correct[ii], 1e-12);
  } // for

  lac_BlockCRSFree(&p_Mat);

} // test_BlockCRS

TEST_MAIN()
