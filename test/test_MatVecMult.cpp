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

TEST_MAIN()
