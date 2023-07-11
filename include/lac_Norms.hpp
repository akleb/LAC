/**
 * @file lac_Norms.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief Computing various matrix/vector norms
 * @version 0.1
 * @date 2023-07-11
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_NORMS_HPP_
#define _LAC_NORMS_HPP_

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


#endif

