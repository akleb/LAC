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

#include "lac_GMRES.hpp"
#include "lac_Norms.hpp"
#include "lac_Structs.hpp"
#include "mpi.h"
#include "unit_test_framework.h"

struct Identity : public lac_MatrixFreeLinearSystem{
public:
  Identity(int n_) : n(n_) {};

  int MatVecProd(const double *x, double *Ax) override{
    std::memcpy(Ax, x, sizeof(double)*n);
    return lac_OK;
  } // MatVecProd

  int RightPrecondition(const double *x, double *Mx) override{
    std::memcpy(Mx, x, sizeof(double)*n);
    return lac_OK;
  } // RightMInverse
private:
  int n;
};

#define N 10
TEST(test_identity){
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double b[N];
  srand(rank);
  for (int i = 0; i < N; ++i)
    b[i] = (double)rand() / RAND_MAX;

  Identity identity(N);

  double x[N];
  double b_test[N];
  double norm, ierr;

  std::memset(x, 0, sizeof(double)*N);
  ierr = lac_GMRES(&identity, b, N, 10, 1e-14, 1e4, false, x, size, rank);

  ASSERT_EQUAL(ierr, lac_OK);

  identity.MatVecProd(x, b_test);
  for (int ii = 0; ii < N; ++ii)
    b_test[ii] -= b[ii];

  lac_L2Norm(b_test, N, &norm);
  ASSERT_TRUE(norm <= 1e-14);


  for (int i = 0; i < N; ++i){
    ASSERT_ALMOST_EQUAL(x[i], b[i], 1e-10);
  } // for

  return;

} // test_identity
#undef N

#define N 10
TEST(test_error){
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double b[N];
  srand(rank);
  for (int i = 0; i < N; ++i)
    b[i] = (double)rand() / RAND_MAX;

  Identity identity(N);

  double x[N];
  double b_test[N];
  double norm, ierr;

  std::memset(x, 0, sizeof(double)*N);
  ierr = lac_GMRES(&identity, b, N, 10, 1e-14, 0, false, x, size, rank);

  ASSERT_EQUAL(ierr, lac_MAX_ITER);

  return;

} // test_identity
#undef N

struct Diffusion : lac_MatrixFreeLinearSystem{
private:
  int n, rank, size;

  void GetHalos(const double *x, double *halo_L, double *halo_R){
    if (size == 1){
      *halo_L = x[n - 1];
      *halo_R = x[0];
    } // if
    else{
      MPI_Request req_in[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};
      MPI_Request req_out[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

      if (rank == 0){
        MPI_Irecv(halo_L, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, req_in + 0);
        MPI_Isend(x + 0, 1, MPI_DOUBLE, size - 1, 1, MPI_COMM_WORLD, req_out + 0);
      } // if
      else{
        MPI_Irecv(halo_L, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, req_in + 0);
        MPI_Isend(x + 0, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, req_out + 0);
      } // else
        
      MPI_Irecv(halo_R, 1, MPI_DOUBLE, (rank + 1) % size, 1, MPI_COMM_WORLD, req_in + 1);
      MPI_Isend(x + n - 1, 1, MPI_DOUBLE, (rank + 1) % size, 0, MPI_COMM_WORLD, req_out + 1);

      MPI_Waitall(2, req_in, MPI_STATUSES_IGNORE);
      MPI_Waitall(2, req_out, MPI_STATUSES_IGNORE);
    } // else
    
    return;

  } // GetHalos
    

public:
  Diffusion(int n_, int size_, int rank_): n(n_), size(size_), rank(rank_) {
    if (n<=2) exit(1);
  } // Constructor
    
  int MatVecProd(const double *x, double *b) override{
    double halo_R, halo_L;
    GetHalos(x, &halo_L, &halo_R);

    b[0] = 2*x[0] - x[1] - halo_L;
    if (rank == size - 1)
      b[n-1] = 2*x[n-1] -x[n-2] + halo_R;
    else 
      b[n - 1] = 2*x[n-1] - x[n-2] - halo_R;
    for (int i = 1; i < n-1; ++i)
      b[i] = 2*x[i] - x[i-1] - x[i+1];
    return lac_OK;
  } // MatVecProd

  // Diagonal preconditioning
  int RightPrecondition(const double *x, double *Mx) override{
    double halo_R, halo_L;
    GetHalos(x, &halo_L, &halo_R);
    
    Mx[0] = 0.5*x[0] + 0.25*(halo_L + x[1]);
    Mx[n-1] = 0.5*x[n-1] + 0.25*(x[n-2] + halo_R);
    for (int i = 1; i < n - 1; ++i)
      Mx[i] = 0.5*x[i] + 0.25*(x[i-1] + x[i+1]);
    return lac_OK;
  } // RightMInverse
};

#define N 10
TEST(test_diffusion){
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  double b[N];
  std::memset(b, 0, sizeof(double)*N);

  Diffusion diff(N, size, rank);

  double x[N];
  double b_test[N];
  double norm;
  int ierr;

  srand(rank);
  for (int i = 0; i < N; ++i)
    x[i] = (double)rand() / RAND_MAX;

  ierr = lac_GMRES(&diff, b, N, 10, 1e-14, 1e4, false, x, size, rank);
  ASSERT_EQUAL(ierr, lac_OK);

  diff.MatVecProd(x, b_test);
  for (int ii = 0; ii < N; ++ii)
    b_test[ii] -= b[ii];

  lac_L2Norm(b_test, N, &norm);
  ASSERT_TRUE(norm <= 1e-13);

  for (int i = 0; i < N; ++i){
    ASSERT_ALMOST_EQUAL(x[i], 0, 1e-10);
  } // for
    
  srand(rank);
  for (int i = 0; i < N; ++i)
    x[i] = (double)rand() / RAND_MAX;
  ierr = lac_GMRES(&diff, b, N, 10, 1e-14, 1e4, true, x, size, rank);
  ASSERT_EQUAL(ierr, lac_OK);

  diff.MatVecProd(x, b_test);
  for (int ii = 0; ii < N; ++ii)
    b_test[ii] -= b[ii];

  lac_L2Norm(b_test, N, &norm);
  ASSERT_TRUE(norm <= 1e-13);

  for (int i = 0; i < N; ++i){
    ASSERT_ALMOST_EQUAL(x[i], 0, 1e-10);
  } // for

  return;

} // test_diffusion
#undef N

TEST_PARALLEL_MAIN()

