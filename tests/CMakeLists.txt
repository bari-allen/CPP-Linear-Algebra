enable_testing()

add_executable(matrix_test matrix_test.cpp)
add_executable(vector_tests vector_tests.cpp)

target_include_directories(matrix_test PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)
target_include_directories(vector_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(matrix_test GTest::gtest_main imageMat_lib vector_exception_lib)
target_link_libraries(vector_tests GTest::gtest_main imageMat_lib vector_exception_lib)

target_include_directories(matrix_test PRIVATE ${gtest_SOURCE_DIR}/include)
target_include_directories(vector_tests PRIVATE ${gtest_SOURCE_DIR}/include)

set_target_properties(matrix_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set_target_properties(vector_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include(GoogleTest)
gtest_discover_tests(matrix_test)
gtest_discover_tests(vector_tests)