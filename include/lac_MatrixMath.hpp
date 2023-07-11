/**
 * @file lac_MatrixMath.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Matrix math operations
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_MATRIX_MATH_HPP_
#define _LAC_MATRIX_MATH_HPP_

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


#endif
