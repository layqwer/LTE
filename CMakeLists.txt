cmake_minimum_required(VERSION 2.8.4)
include_directories(${CMAKE_SOURCE_DIR}/include)
LIST (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

project(CellSearch)

# 添加源文件
set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/searcher.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lte_lib.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/itpp_ext.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/constants.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/common.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/macros.cpp
)

# 添加你的源文件
add_executable(CellSearch ${SOURCES})

# 设置C++标准
set_property(TARGET CellSearch PROPERTY CXX_STANDARD 11)

# 查找IT++库
FIND_PACKAGE( Boost COMPONENTS system REQUIRED )
INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}")
FIND_PACKAGE( ITPP REQUIRED )
INCLUDE_DIRECTORIES("${ITPP_INCLUDE_DIR}")
message("${ITPP_INCLUDE_DIR}")
FIND_PACKAGE( FFTW REQUIRED )
INCLUDE_DIRECTORIES("${FFTW_INCLUDE_DIR}")
FIND_PACKAGE( Curses REQUIRED )
INCLUDE_DIRECTORIES("${CURSES_INCLUDE_DIR}")

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")
INCLUDE_DIRECTORIES("${PROJECT_BINARY_DIR}")

target_link_libraries(CellSearch ${Boost_LIBRARIES})
target_link_libraries(CellSearch ${ITPP_LIBRARIES})
target_link_libraries(CellSearch ${FFTW_LIBRARIES})
target_link_libraries(CellSearch ${CURSES_LIBRARIES})