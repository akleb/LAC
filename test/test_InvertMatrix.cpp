/**
 * @file test_InvertMatrix.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test that we can invert matrices 
 * @version 0.1
 * @date 2023-07-26
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "unit_test_framework.h"
#include "lac_Direct.hpp"
#include "lac_Error.hpp"

TEST(test_identity) {
  double A[16] = {1, 0, 0, 0, 
                  0, 1, 0, 0, 
                  0, 0, 1, 0,
                  0, 0, 0, 1};
  double A_inv[16];

  int ierr = lac_InvertMatrix(A, 4, A_inv);
  if (ierr != lac_OK){ASSERT_TRUE(false)};

  for (int ii = 0; ii < 16; ++ii) {
    ASSERT_ALMOST_EQUAL(A_inv[ii], A[ii], 1e-10);
  } // for

  return;

} // test_no_perm

TEST(test_2x2) {
  double A[4] = {1, -3, 2, -1};
  double correct[4] = {-1./5, 3./5, -2./5, 1./5};

  double A_inv[4];

  int ierr = lac_InvertMatrix(A, 2, A_inv);
  if (ierr != lac_OK){ASSERT_TRUE(false)};

  for (int ii = 0; ii < 4; ++ii) {
    ASSERT_ALMOST_EQUAL(A_inv[ii], correct[ii], 1e-10);
  } // for

  return;

} // test_2x2
 
TEST(test_3x3) {
  double A[9] = {-1, 7, 8, 
                  4, 2, 3, 
                  3, 0, 1};
  double correct[9] = {-2./15, 7./15, -1./3,
                       -1./3, 5./3, -7./3,
                        2./5, -7./5, 2};
  double x[3];
  double A_inv[9];

  int ierr = lac_InvertMatrix(A, 3, A_inv);
  if (ierr != lac_OK){ASSERT_TRUE(false)};

  for (int ii = 0; ii < 9; ++ii) {
    ASSERT_ALMOST_EQUAL(A_inv[ii], correct[ii], 1e-10);
  } // for

  return;

} // test_no_perm

TEST(test_singular){
  double A[9] = {-1, 7, 8, 4, 2, 3, 3, -21, -24};
  double A_inv[9];

  int ierr = lac_InvertMatrix(A, 3, A_inv);
  ASSERT_EQUAL(ierr, lac_SINGULAR_MATRIX);

  return;

} // test_singular

TEST_MAIN()

