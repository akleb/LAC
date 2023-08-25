/**
 * @file test_MPI.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test the functions that utilize MPI calls
 * @version 0.1
 * @date 2023-08-25
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_Error.hpp"
#include "lac_MPI.hpp"
#include "lac_MatrixMath.hpp"
#include "lac_Norms.hpp"
#include "unit_test_framework.h"

#include <cmath>
#include <cstring>

TEST(test_AllReduce){
  int size; 
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double a[5] = {0.1, -0.1, 2, 0.8, -2.4};

  double reduce[5];

  int ierr = lac_Allreduce(a, reduce, 1, MPI_SUM);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
    
  ASSERT_ALMOST_EQUAL(reduce[0], a[0]*size, 1e-12);

  ierr = lac_Allreduce(a, reduce, 5, MPI_SUM);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  for (int ii = 0; ii < 5; ++ii){
    ASSERT_ALMOST_EQUAL(reduce[ii], a[ii]*size, 1e-12);
  } // for

  return;

} // test_Allreduce
  
TEST(test_L1NormAllReduce){
  int size; 
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double a[5] = {0.1, -0.1, 2, 0.8, -2.4};

  double norm;
  int ierr = lac_L1NormAllReduce(a, 5, &norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
    
  ASSERT_ALMOST_EQUAL(norm, 5.4*size, 1e-12);

  std::memset(a, 0, sizeof(double)*5);
  ierr = lac_L1Norm(a, 5, &norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(norm, 0, 1e-12);

  return;

} // test_L1NormAllReduce

TEST(test_L2NormAllReduce){
  int size; 
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double a[5] = {0.1, -0.1, 2, 0.8, -2.4};
  double l2_norm;
  int ierr = lac_L2NormAllReduce(a, 5, &l2_norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(l2_norm, sqrt(10.42 * size), 1e-12);

  std::memset(a, 0, sizeof(double)*5);
  ierr = lac_L2NormAllReduce(a, 5, &l2_norm);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if

  ASSERT_ALMOST_EQUAL(l2_norm, 0, 1e-12);

  return;

} // test_L2NormAllReduce
  
TEST(test_DotProductAllReduce){
  int size; 
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  double a[4] = {-1, 2, 0, 3};
  double b[4] = {0.5, 1, 0.3, -0.1};

  double dot;
  int ierr = lac_DotProductAllReduce(a, b, 4, &dot);
  if (ierr != lac_OK){
    ASSERT_FALSE(true);
  } // if
  ASSERT_ALMOST_EQUAL(dot, 1.2*size, 1e-12);

  return;

} // test_basic
  
TEST_PARALLEL_MAIN();

