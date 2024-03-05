# VanitySearch
A version support custom range scanning and multi address scanning.

This is a modified version of VanitySearch by [JeanLucPons](https://github.com/JeanLucPons/VanitySearch/).

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
    ./vanitysearch -t 0 -gpu -gpuId 0 -i in.txt -o out.txt --keyspace 20000000000000000:+FFFFFFFFFF
    ```

    ```
    in.txt
    1KZegRZVRkVPUcR6iCZWukdeeP6xuEwpo9 /privatekey:20000000000000001
    1JWHHMMjU4YjVd4gU6jVw3zPMsGKMJUiQw /privatekey:20000000100000000
    13zb1hQbWVsc2S7ZTZnP2G4undNNpdh5so /targetaddress
    ```
