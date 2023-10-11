add_unit_test(test_DirectSolve)
add_unit_test(test_DotProduct)
add_unit_test(test_InvertMatrix)
add_unit_test(test_MatVecMult)
add_unit_test(test_Norms)

# mpi tests
add_unit_test_mpi(test_GMRES "2 3 4")
add_unit_test_mpi(test_MPI "2 3 4")
