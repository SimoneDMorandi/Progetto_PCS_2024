list(APPEND DFN_TEST_headers "${CMAKE_CURRENT_SOURCE_DIR}/test.hpp")


set(DFN_TEST_sources ${DFN_TEST_sources} PARENT_SCOPE)
set(DFN_TEST_headers ${DFN_TEST_headers} PARENT_SCOPE)
set(DFN_TEST_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)


#list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/test.hpp")


#set(src_test_sources ${src_sources} PARENT_SCOPE)
#set(src_test_headers ${src_headers} PARENT_SCOPE)
#set(src_test_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
