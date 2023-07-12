/**
 * @file test_DotProduct.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief test the Dot product function
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "unit_test_framework.h"
#include "lac_MatrixMath.hpp"
#include "lac_Error.hpp"
#include <cstring>

TEST(test_zero){

  const int N = 100;
  double *a = new double[N];
  double *b = new double[N];

  srand(0);
  for (int ii = 0; ii < N; ++ii)
    a[ii] = (double)rand() / RAND_MAX;
  std::memset(b, 0, sizeof(double)*N);

  double dot;

  int ierr = lac_DotProduct(a, b, N, &dot);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_ALMOST_EQUAL(dot, 0, 1e-12);

  ierr = lac_DotProduct(b, a, N, &dot);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_ALMOST_EQUAL(dot, 0, 1e-12);

  delete[] a;
  delete[] b;
  return;

} // test_zero

TEST(test_transpose){

  const int N = 100;
  double *a = new double[N];
  double *b = new double[N];

  srand(0);
  for (int ii = 0; ii < N; ++ii){
    a[ii] = (double)rand() / RAND_MAX;
    b[ii] = (double)rand() / RAND_MAX;
  } // for

  double dot, dot_T;

  int ierr = lac_DotProduct(a, b, N, &dot);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ierr = lac_DotProduct(b, a, N, &dot_T);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(dot, dot_T, 1e-12);

  delete[] a;
  delete[] b;
  return;

} // test_transpose

TEST(test_basic){
  
  double a[4] = {-1, 2, 0, 3};
  double b[4] = {0.5, 1, 0.3, -0.1};

  double dot;
  int ierr = lac_DotProduct(a, b, 4, &dot);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_ALMOST_EQUAL(dot, 1.2, 1e-12);

  return;

} // test_basic

TEST_MAIN()

