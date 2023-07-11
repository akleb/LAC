macro(add_unit_test name)
    add_executable(${name}
                    ${name}.cpp)
    
    # add compilier options for making tests with debug flags 
    target_compile_options(${name} PRIVATE -g3
                                   PRIVATE -O0)

    target_include_directories(${name} PRIVATE ${LAC_INC_DIR} ${MPI_CXX_INCLUDE_DIRS})
    # link to the appropriate libraries
    target_link_libraries(${name} ${LAC_LIBS})
    
    add_test(NAME ${name}
            COMMAND ${name}
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endmacro(add_unit_test)
