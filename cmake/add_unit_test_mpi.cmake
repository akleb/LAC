macro(add_unit_test_mpi name nps)
	add_unit_test(${name})

	IF (NOT MEMCHECK)
		SET(procs ${nps})
		SEPARATE_ARGUMENTS(procs)
		FOREACH(np ${procs})
			IF (${MPIEXEC_MAX_NUMPROCS} LESS ${np})
				add_test(NAME ${name}_np${np}
						COMMAND ${MPIEXEC_EXECUTABLE} --oversubscribe ${MPIEXEC_NUMPROC_FLAG} ${np} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name}
						WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/input_files
						)
			ELSE()
				add_test(NAME ${name}_np${np}
						COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name}
						WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/input_files
						)
			ENDIF (${MPIEXEC_MAX_NUMPROCS} LESS ${np})
			
		ENDFOREACH()
	ENDIF (NOT MEMCHECK)

endmacro(add_unit_test_mpi)