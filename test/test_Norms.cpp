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

  return;

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

  return;

} // test_L2

TEST(test_Parallel){
  int size, rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double a[5];
  for (int ii = 0; ii < 5; ++ii){
    a[ii] = 5 * rank + ii;
    if (ii % 2 == 0) a[ii] *= -1;
  } // for

  double l1 = 0;
  double l2 = 0;
  for (int ii = 0; ii < 5 * size; ++ii){
    l1 += ii;
    l2 += ii * ii;
  } // for
  l2 = sqrt(l2);

  double norm;
  int ierr = lac_L1NormAllReduce(a, 5, &norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_EQUAL(l1, norm);
  ierr = lac_DeterministicL1NormAllReduce(a, 5, &norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_EQUAL(l1, norm);
  ierr = lac_L2NormAllReduce(a, 5, &norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_ALMOST_EQUAL(l2, norm, 1e-13);
  ierr = lac_DeterministicL2NormAllReduce(a, 5, &norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_ALMOST_EQUAL(l2, norm, 1e-13);

  return;

} // test_Parallel

TEST_PARALLEL_MAIN()
