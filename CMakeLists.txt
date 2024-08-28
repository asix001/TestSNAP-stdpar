cmake_minimum_required (VERSION 3.12)
project(matvec
  LANGUAGES C CXX
  )


# don't allow in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(STATUS "Warning! Building from the source directory is not allow")
  message(STATUS "Remove 'CMakeCache.txt' and 'CMakeFiles' and build from a separate directory")
  message(ERROR "In-source build")
endif()

SET(MyTarget test_snap)
FILE(GLOB sources
    ${CMAKE_SOURCE_DIR}/src/*.cpp
    ${CMAKE_SOURCE_DIR}/src/*.h
  )

ADD_EXECUTABLE(
  ${MyTarget}
  ${sources}
  )

target_compile_definitions(${MyTarget} PRIVATE REFDATA_TWOJ=${ref_data})

### NVIDIA GPU ###
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "NVHPC")
  find_package(CUDA REQUIRED) 
  add_subdirectory(NVTX/c)
  set(CMAKE_CXX_FLAGS "-O3 -Wall -gpu=cc80 -stdpar=gpu -gpu=mem:unified -std=c++20") #-g -gpu=debug") # -Minfo")
  target_link_libraries(${MyTarget} PRIVATE nvtx3-cpp) # use nvtx tool
endif()

### AMD GPU ###
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "HIP")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -std=c++17 --hipstdpar --hipstdpar-path=/home/gp.sc.cc.tohoku.ac.jp/joanna/research/himeno/stdpar/roc-stdpar/include --offload-arch=gfx90a:xnack+")
endif()
