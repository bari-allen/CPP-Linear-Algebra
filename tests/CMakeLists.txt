enable_testing()

add_executable(matrix_test matrix_test.cpp)

target_link_libraries(matrix_test GTest::gtest_main)

target_include_directories(matrix_test PRIVATE ${gtest_SOURCE_DIR}/include)

set_target_properties(matrix_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

include(GoogleTest)
gtest_discover_tests(matrix_test)
