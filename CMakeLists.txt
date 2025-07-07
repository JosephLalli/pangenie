cmake_minimum_required(VERSION 3.5)
project(PanGenie)

set (CMAKE_CXX_STANDARD 20)

# Set install prefix if not specified
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/home/lalli/usr/local" CACHE PATH "Default install prefix" FORCE)
endif()

# Option to enable KFF file format support
option(ENABLE_KFF_SUPPORT "Enable KFF file format support" OFF)

enable_testing()

find_package(PkgConfig REQUIRED)
find_package(ZLIB REQUIRED)

# Find cereal using modern CMake approach
find_path(CEREAL_INCLUDE_DIR 
    NAMES cereal/cereal.hpp
    PATHS ${CMAKE_PREFIX_PATH}/include
          /usr/local/include
          /usr/include
)

if(CEREAL_INCLUDE_DIR)
    message(STATUS "Found cereal headers: ${CEREAL_INCLUDE_DIR}")
else()
    message(FATAL_ERROR "cereal headers not found. Set CEREAL_ROOT environment variable or CMAKE_PREFIX_PATH")
endif()

pkg_check_modules(JELLYFISH REQUIRED jellyfish-2.0)

# Handle KFF support option
if(ENABLE_KFF_SUPPORT)
    message(STATUS "KFF support requested by user")
    
    # Find KFF C++ API using modern CMake approach
    set(KFF_ROOT "$ENV{KFF_ROOT}" CACHE PATH "KFF root directory")
    find_path(KFF_INCLUDE_DIR 
        NAMES kff_io.hpp
        HINTS $ENV{KFF_ROOT}/include
              ${KFF_ROOT}/include
        PATHS ${CMAKE_PREFIX_PATH}/include
              /usr/local/include
              /usr/include
        PATH_SUFFIXES kff)
    
    find_library(KFF_LIBRARY 
        NAMES kff
        HINTS $ENV{KFF_ROOT}/lib
              ${KFF_ROOT}/lib
        PATHS ${CMAKE_PREFIX_PATH}/lib
              /usr/local/lib
              /usr/lib)

    if(KFF_INCLUDE_DIR AND KFF_LIBRARY)
        set(KFF_FOUND TRUE)
        message(STATUS "Found KFF C++ API: ${KFF_LIBRARY}")
    else()
        set(KFF_FOUND FALSE)
        message(FATAL_ERROR "KFF C++ API not found but ENABLE_KFF_SUPPORT=ON. Please install KFF C++ API or set ENABLE_KFF_SUPPORT=OFF")
    endif()
else()
    set(KFF_FOUND FALSE)
    message(STATUS "KFF support disabled (use -DENABLE_KFF_SUPPORT=ON to enable)")
endif()

link_directories(${JELLYFISH_LIBRARY_DIRS} ${ZLIB_LIBRARY_DIRS} ${cereal_LIBRARY_DIRS})

if(KFF_FOUND)
    include_directories(${KFF_INCLUDE_DIR})
    link_directories(${KFF_LIBRARY_DIRS})
    add_definitions(-DKFF_SUPPORT)
endif()

include_directories(${JELLYFISH_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS} ${CEREAL_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/src)
add_compile_options(${JELLYFISH_CFLAGS_OTHER})

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-long-long -Wno-unknown-pragmas -O3 -I${ZLIB_INCLUDE_DIRS}")

add_subdirectory(src)
add_subdirectory(tests)

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "${dir}")
endforeach()
