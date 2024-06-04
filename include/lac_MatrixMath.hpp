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

#include "lac_BlockCRSMatrix.hpp"

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
 * @brief Computes the matrix vector product Ax -> b
 *
 * @param A an n x m row major matrix
 * @param x a length m vector
 * @param n the number of rows in A
 * @param m the number of columns in A
 * @param b a length n vector set to the product of Ax
 */
int lac_MatVecMultRow(const double *A, const double *x, const int n, 
                        const int m, double *b);

/**
 * @brief Computes the matrix vector product Ax -> b
 *
 * @param A an A->n * block_n x A->m * block_m block CRS matrix
 * @param block_n the number of rows in each A block 
 * @param block_m the number of columns in each A block
 * @param x a length m * block_m vector
 * @param b a length n * block_n vector set to the product of Ax
 */
int lac_MatVecMultBlockCRS(const lac_BlockCRSMatrix *A, const int block_n, 
                           const int block_m, const double *x, double *b);

/**
 * @brief Computes the Product C = AB
 *
 * @param A matrix that is n by k
 * @param B matrix that is k by m
 * @param n number of rows in A and C
 * @param m number of columns in B and C
 * @param k number of colunns in A and rows in B
 * @param C the matrix product
 */
int lac_MatRowMatRowMult(const double *A, const double *B, const int n, 
                         const int m, const int k, double *C);
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
 * @brief computes the dot product of a and b across all processors
 *
 * @param a first vector
 * @param b second vector
 * @param n the size of the vectors
 * @param dot the dot product of a and b
 */
int lac_DotProductAllReduce(const double *a, const double *b, const int n, double *dot);

#endif
