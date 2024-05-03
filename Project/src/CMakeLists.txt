list(APPEND DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)
list(APPEND DFN_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)

list(APPEND DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/Traccia.hpp)

list(APPEND DFN_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(DFN_sources ${DFN_sources} PARENT_SCOPE)
set(DFN_headers ${DFN_headers} PARENT_SCOPE)
set(DFN_includes ${DFN_includes} PARENT_SCOPE)
