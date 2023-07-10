/**
 * @file test_GMRES.cpp
 * @author Alex Kleb (akleb@umich.edu)
 * @brief test the GMRES algorithm
 * @version 0.1
 * @date 2023-06-08
 *
 * @copyright Copyright (c) 2023
 *
 */

#include "lac_LinearAlgebra.hpp"
#include "unit_test_framework.h"

struct Identity{
private:
  int n;

public:
  Identity(int n_) : n(n_){};

  int MatVecProd(const double *x, double *b){
    std::memcpy(b, x, sizeof(double)*n);
    return lac_OK;
  } // MatVecProd

  int RightPrecondition(double *x, double *Mx){
    std::memcpy(Mx, x, sizeof(double)*n);
    return lac_OK;
  } // RightMInverse
};

#define N 10
TEST(test_identity){
  double b[N];
  srand(0);
  for (int i = 0; i < N; ++i)
    b[i] = (double)rand() / RAND_MAX;

  Identity identity(N);

  double x[N];
  std::memset(x, 0, sizeof(double)*N);
  lac_GMRES(&identity, b, N, 10, 1e-14, false, x);

  for (int i = 0; i < N; ++i){
    ASSERT_ALMOST_EQUAL(x[i], b[i], 1e-10);
  } // for

  return;

} // test_identity
#undef N

struct Diffusion{
private:
  int n;

public:
  Diffusion(int n_): n(n_) {if (n<=2) exit(1);};

  int MatVecProd(const double *x, double *b){
    b[0] = 2*x[0] - x[1] - x[n - 1];
    b[n-1] = 2*x[n-1] -x[n-2] + x[0];
    for (int i = 1; i < n-1; ++i)
      b[i] = 2*x[i] - x[i-1] - x[i+1];
    return lac_OK;
  } // MatVecProd

  // Diagonal preconditioning
  int RightPrecondition(double *x, double *Mx){
    
    Mx[0] = 0.5*x[0] + 0.25*(x[n-1] + x[1]);
    Mx[n-1] = 0.5*x[n-1] + 0.25*(x[n-2] + x[0]);
    for (int i = 1; i < n - 1; ++i)
      Mx[i] = 0.5*x[i] + 0.25*(x[i-1] + x[i+1]);
    return lac_OK;
  } // RightMInverse
};

#define N 10
TEST(test_diffusion){
  double b[N];
  std::memset(b, 0, sizeof(double)*N);

  Diffusion diff(N);

  double x[N];
  srand(0);
  for (int i = 0; i < N; ++i)
    x[i] = (double)rand() / RAND_MAX;

  lac_GMRES(&diff, b, N, 10, 1e-14, false, x);

  for (int i = 0; i < N; ++i){
    ASSERT_ALMOST_EQUAL(x[i], 0, 1e-10);
  } // for
    
  srand(0);
  for (int i = 0; i < N; ++i)
    x[i] = (double)rand() / RAND_MAX;
  lac_GMRES(&diff, b, N, 10, 1e-14, true, x);

  for (int i = 0; i < N; ++i){
    ASSERT_ALMOST_EQUAL(x[i], 0, 1e-10);
  } // for


  return;

} // test_identity
#undef N

TEST_MAIN()
