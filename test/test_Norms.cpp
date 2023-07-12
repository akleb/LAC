/**
 * @file test_Norms.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test the Norm functions
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "unit_test_framework.h"
#include "lac_Norms.hpp"
#include "lac_Error.hpp"
#include <cmath>
#include <cstring>

TEST(test_L1){
  double a[5] = {0.1, -0.1, 2, 0.8, -2.4};
  double l1_norm;
  int ierr = lac_L1Norm(a, 5, &l1_norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(l1_norm, 5.4, 1e-12);

  std::memset(a, 0, sizeof(double)*5);
  ierr = lac_L1Norm(a, 5, &l1_norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(l1_norm, 0, 1e-12);


} // test_L1

TEST(test_L2){
  double a[5] = {0.1, -0.1, 2, 0.8, -2.4};
  double l2_norm;
  int ierr = lac_L2Norm(a, 5, &l2_norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(l2_norm, sqrt(10.42), 1e-12);

  std::memset(a, 0, sizeof(double)*5);
  ierr = lac_L2Norm(a, 5, &l2_norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(l2_norm, 0, 1e-12);


} // test_L1

TEST_MAIN()
