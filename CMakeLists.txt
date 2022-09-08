cmake_minimum_required(VERSION 3.20)
project(PanGenie)


set(CMAKE_CXX_STANDARD 17)
set (CMAKE_VERBOSE_MAKEFILE on)
enable_testing()

find_package(PkgConfig REQUIRED)
pkg_check_modules(JELLYFISH REQUIRED jellyfish-2.0)

include_directories(include)

link_directories(${JELLYFISH_LIBRARY_DIRS})

include_directories(${JELLYFISH_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/src)
add_compile_options(${JELLYFISH_CFLAGS_OTHER})

set(CMAKE_CXX_FLAGS "-O3 -march=native -W -Wall -pedantic -fiopenmp -fma -funroll-loops")

add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "${dir}")
endforeach()
