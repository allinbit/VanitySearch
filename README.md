# VanitySearch
A version support custom range scanning and multi address scanning.

# Build
- Edit the makefile and set up the appropriate CUDA SDK and compiler paths for nvcc.
    ```
    ccap=86
    
    ...
    
    CXX        = g++-9
    CUDA       = /usr/local/cuda-11.7
    CXXCUDA    = /usr/bin/g++-9
    ```

 - Build:
    ```
    $ make all
    ```

# Usage
- Example for bitcoin puzzle 66
    ```
    ./vanitysearch -t 0 -gpu -gpuId 0 -i in.txt -o out.txt --keyspace 2xxxxxx0000000000:+FFFFFFFFFF
    ```
