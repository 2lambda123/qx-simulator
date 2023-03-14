cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

add_qx_test(BasisVectorTest BasisVectorTest.cpp .)
target_link_libraries(BasisVectorTest PRIVATE doctest)

add_qx_test(MatrixTest MatrixTest.cpp .)
target_link_libraries(MatrixTest PRIVATE doctest)

add_qx_test(MixedStateSimplifierTest MixedStateSimplifierTest.cpp .)
target_link_libraries(MixedStateSimplifierTest PRIVATE doctest)

add_qx_test(MixedStateTest MixedStateTest.cpp .)
target_link_libraries(MixedStateTest PRIVATE doctest)

add_qx_test(RandomTest RandomTest.cpp .)
target_link_libraries(RandomTest PRIVATE doctest)

add_qx_test(IntegrationTest IntegrationTest.cpp .)
target_link_libraries(IntegrationTest PRIVATE doctest)