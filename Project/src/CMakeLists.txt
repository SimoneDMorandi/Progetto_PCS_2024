list(APPEND DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND DFN_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)
list(APPEND DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/Structures.hpp)

#list(APPEND DFN_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(DFN_sources ${DFN_sources} PARENT_SCOPE)
set(DFN_headers ${DFN_headers} PARENT_SCOPE)
set(DFN_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

#list(APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp")
#list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp")
#list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/Structures.hpp")


#set(src_sources ${src_sources} PARENT_SCOPE)
#set(src_headers ${src_headers} PARENT_SCOPE)
#set(src_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)