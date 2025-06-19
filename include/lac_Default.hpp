/**
 * @file lac_Default.hpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief The Default functions
 * @version 0.1
 * @date 2023-06-09
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef _LAC_DEFAULT_HPP_
#define _LAC_DEFAULT_HPP_

#include <cstdio>

#define LAC_ERR   printf("===> ERR: ");printf
#define LAC_NOTE  printf("        o ");printf
#define LAC_CONT  printf("          ");printf
#define LAC_WARN  printf("==> WARN: ");printf


/**
 * @brief Gathers an array split up across multiple processors on comm world
 * MPI_COMM WORLD onto rank 0, all other processors get *full_a = nullptr. They
 * get gathered in rank order. This should probably only be used for debugging
 * things in parallel to make them determinsitic.
 *
 * @param a the portion of the array on this processor
 * @param n the number of entries on this processor, not necesssarily equal to
 * all other processors
 * @param full_a the fully gather array in processor rank on rank zero all other
 * ranks set to nullptr
 * @param full_n set to the number of entries in the full_a on rank 0 all others,
 * nothing
 */
int _lac_GatherFullArray(const double *a, const int n, double **full_a, int *full_n);

#endif

