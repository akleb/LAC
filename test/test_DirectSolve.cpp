/**
 * @file test_DirectSolve.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test that we can solve linear systems
 * @version 0.1
 * @date 2023-07-26
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "unit_test_framework.h"
#include "lac_Direct.hpp"
#include "lac_Error.hpp"

TEST(test_no_perm) {
  double A[9] = {4, 2, 3, -1, 7, 8, 3, 0, 1};
  double b[3] = {1, -4, 3};
  double correct[3] = {0, -4, 3};
  double x[3];

  int ierr = lac_DirectSolve(A, b, 3, x);
  if (ierr != lac_OK){ASSERT_TRUE(false)};

  for (int ii = 0; ii < 3; ++ii) {
    ASSERT_ALMOST_EQUAL(x[ii], correct[ii], 1e-10);
  } // for

  return;

} // test_no_perm

TEST(test_perm) {
  double A[9] = {-1, 7, 8, 4, 2, 3, 3, 0, 1};
  double b[3] = {-4, 1, 3};
  double correct[3] = {0, -4, 3};
  double x[3];

  int ierr = lac_DirectSolve(A, b, 3, x);
  if (ierr != lac_OK){ASSERT_TRUE(false)};

  for (int ii = 0; ii < 3; ++ii) {
    ASSERT_ALMOST_EQUAL(x[ii], correct[ii], 1e-10);
  } // for

  return;

} // test_no_perm

TEST(test_singular){
  double A[9] = {-1, 4, 8, 4, 2, 4, 2, -8, -16};
  double b[3] = {-4, 1, 3};
  double x[3];

  int ierr = lac_DirectSolve(A, b, 3, x);
  ASSERT_EQUAL(ierr, lac_SINGULAR_MATRIX);

  return;

} // test_singular

TEST_MAIN()

