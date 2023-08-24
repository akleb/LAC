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
#include "lac_Norms.hpp"
#include <mpi.h>
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
  
#ifndef _LAC_GMRES_CLEANUP
#define _LAC_GMRES_CLEANUP { \
        delete[] r;           \
        delete[] K;           \
        delete[] H;           \
        delete[] e;           \
        delete[] F;           \
        delete[] y;           \
        delete[] temp;        \
        return ierr;          \
    }
int lac_GMRES(lac_MatrixFreeLinearSystem *obj, const double *b, const int n, const int nRst, 
               const double tol, const bool precondition, double *x,
               const int size, const int rank, bool verbose, int *iters){

  if (verbose && rank != 0) verbose = false;
  if (iters) *iters = 0;

  // Used to update the value of x
  double *r = new double[n];
  // K and H are both column major
  // initialize the Krylov subspace
  double *K = new double[n*(nRst+1)];
  // initialize the Hessenberg matrix
  double *H = new double[nRst*(nRst+1)];
  double *e = new double[nRst + 1];
  // Stores the previous Givens rotation
  double *F = new double[2*nRst];
  // used to compute the LS operation after building the subspace
  double *y = new double[nRst];

  double *temp = new double[n];

  int ierr;
  // compute the initial linear residual norm
  obj->MatVecProd(x, temp);
  for (int ii = 0; ii < n; ++ii)
    temp[ii] -= b[ii];

  double r_norm, init_r_norm;
  ierr = (size > 1) ? lac_L2NormAllReduce(temp, n, &init_r_norm) : 
                               lac_L2Norm(temp, n, &init_r_norm);
  if (ierr != lac_OK) _LAC_GMRES_CLEANUP; 
  r_norm = init_r_norm;

  // Set the convergence criteria relative to initial linear residual
  const double target_norm = init_r_norm * tol;

  // Start the GMRES iteration
  int nOuter = 0;
  while (r_norm > target_norm){
    std::memset(K, 0, sizeof(double)*n*(nRst+1));
    std::memset(H, 0, sizeof(double)*nRst*(nRst+1));
    std::memset(e, 0, sizeof(double)*(nRst + 1));
    std::memset(F, 0, sizeof(double)*2*nRst);

    // Arnoldi iteration, both K and H are column major
    // need the first vector as Ax - b or AP^-1 x - b
    if (precondition){
      ierr = obj->RightPrecondition(x, temp);
      if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
      ierr = obj->MatVecProd(temp, K);
      if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
    } // if
    else{
      ierr = obj->MatVecProd(x, K);
      if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
    } // else
    for (int ii = 0; ii < n; ++ii)
      K[ii] = b[ii] - K[ii];
    // compute the current residual from this and initialize e with it
    ierr = (size > 1) ? lac_L2NormAllReduce(K, n, &r_norm) :
                                 lac_L2Norm(K, n, &r_norm);
    if (ierr != lac_OK) _LAC_GMRES_CLEANUP;

    if (r_norm <= 1e-16)
      break;
    for (int ii = 0; ii < n; ++ii)
      K[ii] /= r_norm;
    e[0] = r_norm;

    // loop over the Arnoldi vectors to add
    int actual_nRst = nRst;
    for (int col = 0; col < nRst; ++col){
      // column starts with the matrix multiply
      if (precondition){
        ierr = obj->RightPrecondition(K + n*col, temp);
        if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
        ierr = obj->MatVecProd(temp, K + n*(col+1));
        if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
      } // if
      else{
        ierr = obj->MatVecProd(K + n*col, K + n*(col+1));
        if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
      } // else

      // MGS to make it orthogonal
      for (int row = 0; row < col + 1; ++row){
        double dot;
        ierr = lac_DotProduct(K + n*row, K + n*(col+1), n, &dot);
        if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
        if (size > 1) 
          MPI_Allreduce(&dot, H + (nRst+1)*col + row, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        else
          H[(nRst+1)*col + row] = dot;

        for (int ii = 0; ii < n; ++ii)
          K[n*(col+1) + ii] -= H[(nRst+1)*col + row] * K[n*row + ii];
      } // for
      ierr = (size > 1) ? lac_L2NormAllReduce(K + n*(col+1), n, H + (nRst+1)*col + col + 1) :
                                   lac_L2Norm(K + n*(col+1), n, H + (nRst+1)*col + col + 1);
      if (ierr != lac_OK) _LAC_GMRES_CLEANUP;

      // normalize the Krylov subspace column
      for (int row = 0; row < n; ++row){
        K[n*(col+1) + row] /= H[(nRst+1)*col + col + 1];
      } // for

      // Perform the Givens rotation on the last column we added to H for Least
      // Squares
      ierr = lac_GivensRotation(col, H + (nRst+1)*col, e, F);
      if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
      if (iters) (*iters)++;

      // check if we have already converged
      r_norm = fabs(e[col+1]);
      if (verbose){
        CONT("Outer: %6d Inner: %6d Residual: %16.10e\n", nOuter, col, r_norm);
      } // if
      if (r_norm <= target_norm){
        actual_nRst = col + 1;
        break;
      } // if
    } // for
    
    //perform the back sub on H and e to get the solution to the LS
    for (int col = actual_nRst - 1; col >= 0; col--){
      const double *H_col = H + (nRst + 1)*col;
      // compute this column's variable value
      y[col] = e[col] / H_col[col];
      // subtract it from all the other variable's RHS
      for (int row = col - 1; row >= 0; row--){
        e[row] -= y[col] * H_col[row];
      } // for rows
    } // for columns

    // Update the x vector with the computed update
    ierr = lac_MatVecMultCol(K, y, n, actual_nRst, r);
    if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
    for (int row = 0; row < n; ++row)
      x[row] += r[row];

    nOuter++;

  } // while not converged
    
  if (precondition){
    ierr = obj->RightPrecondition(x, temp);
    if (ierr != lac_OK) _LAC_GMRES_CLEANUP;
    std::memcpy(x, temp, sizeof(double)*n);
  } // if

  delete[] temp;
  delete[] K;
  delete[] H;
  delete[] r;
  delete[] y;
  delete[] e;
  delete[] F;
  if (verbose){
    NOTE("GMRES converged residual: %16.10e\n", r_norm);
  } // if

  return lac_OK;

} // lac_GMRES
#undef _LAC_GMRES_CLEANUP
#else
#error "_LAC_GMRES_CLEANUP already defined somewhere!"
#endif
  
