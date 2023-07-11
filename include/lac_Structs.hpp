/**
 * @file lac_Structs.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Structs used for the Linear Algebra Cookbook
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_STRUCTS_HPP_
#define _LAC_STRUCTS_HPP_

// Interface Parent Struct object for linear systems
struct lac_MatrixFreeLinearSystem{

  /**
   * @brief The matrix vector product class. This performs the linear operator A
   * [size n x n] to the vector x [size n]. 
   *
   * @param x the vector multiplied by A
   * @param Ax the result of Ax
   */
  virtual int MatVecProd(const double *x, double *Ax) = 0;

  /**
   * @brief Right Preconditioning operation on x. This performs the inverse of
   * the linear operator M [size n x n] on x [size n] -> (M^-1)x
   *
   * @param x the vector we are multiplying by M^-1
   * @param Mx the resulting vector
   */
  virtual int RightPrecondition(const double *x, double *Mx) = 0;

}; // lac_MatrixFreeLinearSystem
  
#endif

