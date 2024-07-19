# TestSNAP-stdpar

Original work: https://github.com/FitSNAP/TestSNAP/tree/master

Implementation based on: https://github.com/FitSNAP/TestSNAP/tree/OpenMP4.5

## Build system
### CMake
Download NVTX tool into code directory or comment out NVTX commands (in `CMakeLists.txt`, `sna.cpp`, and `test_snap.cpp`)

https://github.com/NVIDIA/NVTX?tab=readme-ov-file

#### NVIDIA GPU
```
mkdir build_cuda && cd build_cuda
cmake -S .. -B . -Dref_data=14 -DCMAKE_CXX_COMPILER=nvc++  -DCMAKE_CXX_COMPILER_ID="NVHPC" 
cmake --build .
./test_snap
```    

