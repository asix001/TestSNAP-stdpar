# TestSNAP-stdpar

Original work: https://github.com/FitSNAP/TestSNAP/tree/master

Implementation based on: https://github.com/FitSNAP/TestSNAP/tree/OpenMP4.5

## Compiler versions
AMD clang version 18.0.0 (ROCM 6.2.0)
nvc++ 22.11-0  

## Build system
### CMake
#### NVIDIA GPU
```
mkdir build_cuda && cd build_cuda
cmake -S .. -B . -Dref_data=14 -DCMAKE_CXX_COMPILER=nvc++  -DCMAKE_CXX_COMPILER_ID=NVHPC
cmake --build .
./test_snap
```    
#### AMD GPU
```
mkdir build_hip && cd build_hip
cmake -S .. -B . -Dref_data=14 -DCMAKE_CXX_COMPILER=amdclang++  -DCMAKE_CXX_COMPILER_ID=HIP
cmake --build .
./test_snap
```    
