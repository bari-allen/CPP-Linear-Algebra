if(BUILD_TESTING)
    add_executable(matrix_test matrix_test.cpp)
    target_link_libraries(matrix_test PRIVATE GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(matrix_test)
endif()