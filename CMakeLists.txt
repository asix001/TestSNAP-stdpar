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
target_link_libraries(${MyTarget} PRIVATE nvtx3-cpp) # use nvtx tool

find_package(CUDA REQUIRED)
# find_package(ROCM REQUIRED)

### NVIDIA GPU ###
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "NVHPC")
    set(CMAKE_CXX_FLAGS "-O3 -Wall -gpu=cc80 -stdpar=gpu -std=c++20")
endif()

### AMD GPU ###
if(${CMAKE_CXX_COMPILER_ID} STREQUAL "HIP")
  set(CMAKE_CXX_FLAGS "-O3 -Wall -std=c++17 -v --hipstdpar --hipstdpar-path=/home/gp.sc.cc.tohoku.ac.jp/joanna/research/himeno/stdpar/roc-stdpar/include --offload-arch=gfx90a:xnack+ -B /home/gp.sc.cc.tohoku.ac.jp/joanna/llvm-project/lld/build/bin")

endif()

## roc-stdpar path
# /home/gp.sc.cc.tohoku.ac.jp/joanna/research/himeno/stdpar/roc-stdpar/include
