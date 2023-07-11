/**
 * @file lac_GMRES.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Generalized Minimum Residual Linear Sovler
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_GMRES.hpp"
#include "lac_Error.hpp"
#include <cmath>

int lac_GivensRotation(const int col, double *H_col, double *e, double *F){
  // perform all the rotations already stored on the new column of H
  for (int ii = 0; ii < col; ++ii){
    const double c = F[2*ii];
    const double s = F[2*ii + 1];
    const double u = H_col[ii];
    const double v = H_col[ii+1];
    H_col[ii] = c*u + s*v;
    H_col[ii + 1] = -s*u + c*v;
  } // for

  // compute the Givens rotation to remove the final off diagonal piece
  const double u = H_col[col];
  const double v = H_col[col + 1];
  const double r = hypot(u, v);
  const double c = F[2*col] = u / r;
  const double s = F[2*col + 1] = v / r;

  // rotate the final part of H
  H_col[col] = c*u + s*v;
  H_col[col + 1] = -s*u + c*v;

  // rotate e as well
  const double e_u = e[col];
  const double e_v = e[col + 1];
  e[col] = c*e_u + s*e_v;
  e[col + 1] = -s*e_u + c*e_v;

  return lac_OK;

} // lac_GivensRotation
  
