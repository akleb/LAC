add_unit_test(test_BlockCRSMatrix)
add_unit_test(test_DirectSolve)
add_unit_test(test_InvertMatrix)
add_unit_test(test_MatMatMult)
add_unit_test(test_MatVecMult)

# mpi tests
add_unit_test_mpi(test_DotProduct "2 3 4")
add_unit_test_mpi(test_GMRES "2 3 4")
add_unit_test_mpi(test_MPI "2 3 4")
add_unit_test_mpi(test_Norms "2 3 4")

# files needed for the tests
file(COPY test_matrix.dat DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
