<!----------------BEGIN-HEADER------------------------------------>
# TestSNAP-stdpar

Original work: https://github.com/FitSNAP/TestSNAP/tree/master
Implementation based on: https://github.com/FitSNAP/TestSNAP/tree/OpenMP4.5

## Build system
### CMake
#### NVIDIA GPU
```
mkdir build_cuda && cd build_cuda
cmake -S .. -B . -Dref_data=14 -DCMAKE_CXX_COMPILER=nvc++  -DCMAKE_CXX_COMPILER_ID="NVHPC" 
cmake --build .
./test_snap
```    
<!-----------------END-HEADER------------------------------------->

