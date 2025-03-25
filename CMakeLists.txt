cmake_minimum_required(VERSION 4.0)
project(PCA
    VERSION 0.0.1
    DESCRIPTION "Take the PCA of a set of images to compress them"
    LANGUAGES CXX)
include(CTest)
add_subdirectory(external/googletest)
add_subdirectory(tests)