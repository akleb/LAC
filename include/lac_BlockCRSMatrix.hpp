/**
 * @file lac_BlockCRSMatrix.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Block CRS Matrix format
 * @version 0.1
 * @date 2023-12-05
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_BLOCK_CRS_MATRIX_HPP_
#define _LAC_BLOCK_CRS_MATRIX_HPP_

#include "lac_Default.hpp"
#include "lac_Error.hpp"

#include <string>

struct lac_BlockCRSMatrix{
  
  double *data;
  int *col_index;
  int *n_col;
  int n_nzero;
  int n_block;
  int n;
  int m;

}; // lac_BlockCRSMatrix

/**
 * @brief Initializes a Block CRS matrix without any data
 *
 * @param pp_Mat generates a pointer to the matrix at the address provided
 * @param n_block the size of a block (amount of data to store per block)
 * @param n the number of rows in the matrix
 * @param m the number of columns in the matrix
 * @param col_index an array of all of the columns in a row which are non-zero
 * each row comes after the previous
 * @param n_col a prefix sum of the number of columns in each row, good for
 * finding indices into col index
 */
void lac_BlockCRSInit(lac_BlockCRSMatrix **pp_Mat, const int n_block, const int n,
                      const int m, const int *col_index, const int *n_col);

/**
 * @brief frees the Block CRS matrix
 *
 * @param pp_Mat the matrx to free
 */
void lac_BlockCRSFree(lac_BlockCRSMatrix **pp_Mat);

/**
 * @brief Adds the data provided to the data in position row/col
 *
 * @param p_Mat the sparse matrix
 * @param row the row in the sparse matrix
 * @param col the col in the sparse matrix
 * @param data the data to add to the specified entry
 */
int lac_BlockCRSAddData(lac_BlockCRSMatrix *p_Mat, const int row, const int col, const double *data);

/**
 * @brief Sets the data provided to the data in position row/col
 *
 * @param p_Mat the sparse matrix
 * @param row the row in the sparse matrix
 * @param col the col in the sparse matrix
 * @param data the data to set to the specified entry
 */
int lac_BlockCRSSetData(lac_BlockCRSMatrix *p_Mat, const int row, const int col, const double *data);

/**
 * @brief Zeros all the data in the matrix
 *
 * @param p_Mat the matrix
 */
int lac_BlockCRSZeroData(lac_BlockCRSMatrix *p_Mat);

/**
 * @brief Sets the data provided to the data in position row/col
 *
 * @param p_Mat the sparse matrix
 * @param row the row in the sparse matrix
 * @param col the col in the sparse matrix
 * @param data the data retrieved from the sparse matrix
 */
int lac_BlockCRSGetData(lac_BlockCRSMatrix *p_Mat, const int row, const int col, double ** data);

/**
 * @brief Transposes the data in the block CRS matrix
 *
 * @param p_Mat the matrix will be transposed in place
 * @param block_n the number of rows in each block of the matrix before the
 * transpose
 * @param block_m the number of cols in each block of the matrix before the
 * transpose
 */
int lac_BlockCRSTranspose(lac_BlockCRSMatrix *p_Mat, const int block_n, const int block_m);

/**
 * @brief performs the ILU0 factorization on a BlockCRS Matrix
 *
 * @param p_A the block CRS matrix
 * @param block_n the size of each block in block_n by block_n
 * @param p_LU the resulting ILU0 factorization, needs to be iniitialized to the
 * same sparsity pattern as p_A
 */
int lac_BlockCRS_ILU0(const lac_BlockCRSMatrix *p_A, const int block_n, lac_BlockCRSMatrix *p_LU);

/**
 * @brief does a foward back sub to compute the solution to LUx = b
 * 
 * @param p_LU the LU factorization in BlockCRS format
 * @param block_n the size of the blocks in the blockCRS matrix
 * @param b the RHS, size block_n * p_LU->n
 * @param x the solution, size block_n * p_LU->n
 */
int lac_BlockCRS_LUForwardBackwardSub(lac_BlockCRSMatrix *p_LU, const int block_n, const double *b, double *x);

/**
 * @brief prints out the block CRS Matrix
 *
 * @param p_Mat the matrix to print
 * @param block_n the number of rows in each block
 * @param block_m the number of columns in each block
 */
int lac_BlockCRSPrint(lac_BlockCRSMatrix *p_Mat, const int block_n, const int block_m);

/**
 * @brief prints out the block CRS Matrix to file
 *
 * @param p_Mat the matrix to print
 * @param filename the name of the file to write the matrix to
 */
int lac_BlockCRSFileWrite(lac_BlockCRSMatrix *p_Mat, const std::string filename);

/**
 * @brief Initializes a block CRS matrix from file
 *
 * @param pp_Mat the matrix that is going to be intialized
 * @param filename the name of the file containing the matrix
 */
int lac_BlockCRSInitFromFile(lac_BlockCRSMatrix **pp_Mat, const std::string filename);

#endif

