/**
 * @file lac_Norms.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Functions for computing matrix/vector norms
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_Error.hpp"
#include "lac_Norms.hpp"
#include <cmath>

int lac_L2Norm(const double *a, const int n, double *norm){
  *norm = 0;
  for (int ii = 0; ii < n; ++ii)
    (*norm) += a[ii]*a[ii];
  (*norm) = sqrt(*norm);

  return lac_OK;

} // lac_L2Norm
 
int lac_L1Norm(const double *a, const int n, double *norm){
  *norm = 0;
  for (int ii = 0; ii < n; ++ii)
    (*norm) += fabs(a[ii]);

  return lac_OK;

} // lac_L1Norm

