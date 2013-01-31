
enable_testing()

find_package(Check REQUIRED)
include_directories(${CHECK_INCLUDE_DIRS})
include_directories(../src)

add_executable(lfft_tests lfft_tests.c)
target_link_libraries(lfft_tests ${CHECK_LIBRARIES} m)

add_test(test_lfft_tests ${PROJECT_BINARY_DIR}/bin/lfft_tests)

