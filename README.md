How to build and install the libraries
======================================


Prerequisites
-------------

- C++ compiler with C++11 support (e.g. g++ >= 4.8.2)
- ROOT (>= 5.34.30), http://root.cern.ch


Build with cmake
----------------

Checkout the code using the following command:

    git clone https://github.com/star-bnl/star-vertex.git

Compile and build libraries:

    cd star-vertex/
    mkdir build && cd build
    cmake -D CMAKE_INSTALL_PREFIX=./ ../
    make install
