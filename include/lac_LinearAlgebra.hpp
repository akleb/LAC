/**
 * @file lac_LinearAlgebra.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Linear Algebra package in LAC
 * @version 0.1
 * @date 2023-07-06
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_LINEAR_ALGEBRA_HPP_
#define _LAC_LINEAR_ALGEBRA_HPP_

#include <cstring>
#include <cmath>
#include "lac_Default.hpp"
#include "lac_Error.hpp"

/**
 * @brief Computes the matrix vector product Ax -> b
 *
 * @param A an n x m column major matrix
 * @param x a length m vector
 * @param n the number of rows in A
 * @param m the number of columns in A
 * @param b a length n vector set to the product of Ax
 */
int lac_MatVecMultCol(const double *A, const double *x, const int n, 
                        const int m, double *b);
/**
 * @brief computes the dot product of a and b
 *
 * @param a first vector
 * @param b second vector
 * @param n the size of the vectors
 * @param dot the dot product of a and b
 */
int lac_DotProduct(const double *a, const double *b, const int n, double *dot);

/**
 * @brief computes the L2 norm of a
 *
 * @param a the vector
 * @param n the size of the vector
 * @param norm the computed L2 norm
 */
int lac_L2Norm(const double *a, const int n, double *norm);

/**
 * @brief computes the L2 norm of a
 *
 * @param a the vector
 * @param n the size of the vector
 * @param norm the computed L1 norm
 */
int lac_L1Norm(const double *a, const int n, double *norm);

/**
 * @brief Applies previous Givens rotations to H_col and comptues the new
 * rotation to remove the last entry in the column
 *
 * @param col the column index we are acting on inside of H
 * @param H_col the col column in the upper Hessenberg without Givens rotations
 * applied
 * @param e the RHS for the least-squares that needs the new Givens rotation
 * applied to it as well
 * @param F stores all the previous Givens rotations that need to be applied to
 * H_col before computing the new rotation, the new rotation is stored in
 * F[2*col] <= c and F[2*col + 1] <= s
 */
int lac_GivensRotation(const int col, double *H_col, double *e, double *F);

/**
 * @brief computes the solution to Ax = b using a matrix free GMRES approach
 *
 * @param obj an object that contains preconditioner and matrix-vector product
 * functions. 
 *      void T::MatVecProd(const double *v, double *Av)
 *      void T::RightPrecondition(const double *x, double *Mx)
 * @param b the right hand side we are solving for
 * @param n the number of unkowns
 * @param nRst the number of Krylov vectors to build before reseting
 * @param tol the tolerance to converge the residual to
 * @param precondition true to right precondition, false otherwise
 * @param x the solution vector
 * @param verbose true for extra print out, false for less
 */
template <typename T>
int lac_GMRES(T *obj, const double *b, const int n, const int nRst, 
               const double tol, const bool precondition, double *x, 
               bool verbose = false){

  int ierr;
  double r_norm = HUGE_VAL;

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
  int nOuter = 0;
  while (r_norm > tol){
    std::memset(K, 0, sizeof(double)*n*(nRst+1));
    std::memset(H, 0, sizeof(double)*nRst*(nRst+1));
    std::memset(e, 0, sizeof(double)*(nRst + 1));
    std::memset(F, 0, sizeof(double)*2*nRst);
    // Previous Givens rotation set to Identity, so it doesn't rotate

    // Arnoldi iteration, both K and H are column major
    // need the first vector as Ax - b
    obj->MatVecProd(x, K);
    for (int ii = 0; ii < n; ++ii)
      K[ii] = b[ii] - K[ii];
    // compute the current residual from this and initialize e with it
    ierr = lac_L2Norm(K, n, &r_norm);
    if (ierr != lac_OK){
      delete[] r;
      delete[] K;
      delete[] H;
      delete[] e;
      delete[] F;
      delete[] y;
      delete[] temp;
      return ierr;
    } // if
     
    if (r_norm <= tol)
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
        if (ierr != lac_OK){
          delete[] r;
          delete[] K;
          delete[] H;
          delete[] e;
          delete[] F;
          delete[] y;
          delete[] temp;
          return ierr;
        } // if
        ierr = obj->MatVecProd(temp, K + n*(col+1));
        if (ierr != lac_OK){
          delete[] r;
          delete[] K;
          delete[] H;
          delete[] e;
          delete[] F;
          delete[] y;
          delete[] temp;
          return ierr;
        } // if
      } // if
      else{
        ierr = obj->MatVecProd(K + n*col, K + n*(col+1));
        if (ierr != lac_OK){
          delete[] r;
          delete[] K;
          delete[] H;
          delete[] e;
          delete[] F;
          delete[] y;
          delete[] temp;
          return ierr;
        } // if
      } // else
      // MGS to make it orthogonal
      for (int row = 0; row < col + 1; ++row){
        ierr = lac_DotProduct(K + n*row, K + n*(col+1), n, H + (nRst+1)*col + row);
        if (ierr != lac_OK){
          delete[] r;
          delete[] K;
          delete[] H;
          delete[] e;
          delete[] F;
          delete[] y;
          delete[] temp;
          return ierr;
        } // if

        for (int ii = 0; ii < n; ++ii)
          K[n*(col+1) + ii] -= H[(nRst+1)*col + row] * K[n*row + ii];
      } // for
      ierr = lac_L2Norm(K + n*(col+1), n, H + (nRst+1)*col + col + 1);
      if (ierr != lac_OK){
        delete[] r;
        delete[] K;
        delete[] H;
        delete[] e;
        delete[] F;
        delete[] y;
        delete[] temp;
        return ierr;
      } // if
      // normalize the Krylov subspace column
      for (int row = 0; row < n; ++row){
        K[n*(col+1) + row] /= H[(nRst+1)*col + col + 1];
      } // for

      // Perform the Givens rotation on the last column we added to H for Least
      // Squares
      ierr = lac_GivensRotation(col, H + (nRst+1)*col, e, F);
      if (ierr != lac_OK){
        delete[] r;
        delete[] K;
        delete[] H;
        delete[] e;
        delete[] F;
        delete[] y;
        delete[] temp;
        return ierr;
      } // if

      // check if we have already converged
      r_norm = fabs(e[col+1]);
      if (verbose){
        CONT("Outer: %6d Inner: %6d Residual: %16.10e\n", nOuter, col, r_norm);
      } // if
      if (r_norm <= tol){
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
    if (ierr != lac_OK){
      delete[] r;
      delete[] K;
      delete[] H;
      delete[] e;
      delete[] F;
      delete[] y;
      delete[] temp;
      return ierr;
    } // if
    if (precondition){
      ierr = obj->RightPrecondition(r, temp);
      if (ierr != lac_OK){
        delete[] r;
        delete[] K;
        delete[] H;
        delete[] e;
        delete[] F;
        delete[] y;
        delete[] temp;
        return ierr;
      } // if
      std::memcpy(r, temp, sizeof(double)*n);
    } // if
    for (int row = 0; row < n; ++row)
      x[row] += r[row];

    nOuter++;

  } // while not converged

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

#endif
