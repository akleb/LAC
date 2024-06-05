/**
 * @file test_MatMatMult.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Test the Matrix Matrix multiply functions
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

TEST(test_rand){
  double A[6] = {   0,    1,
                 -0.5,    3,
                  5.1, -0.1};

  double B[8] = { 1, -0.2, 3.4, 9,
                 -2, -0.1, 2.5, 8};

  double C[12];
  int ierr = lac_MatRowMatRowMult(A, B, 3, 4, 2, C);
  ASSERT_EQUAL(ierr, lac_OK);

  double corrcet_C[12] = {  -2,  -0.1,   2.5,    8,
                          -6.5,  -0.2,   5.8, 19.5,
                           5.3, -1.01, 17.09, 45.1};

  for (int ii = 0; ii < 12; ++ii){
    ASSERT_ALMOST_EQUAL(C[ii], corrcet_C[ii], 1e-12);
  } // for
  
  return;

} // test_identity

TEST_MAIN()
