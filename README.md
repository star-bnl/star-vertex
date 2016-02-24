How to build and install the libraries
======================================


Prerequisites
-------------

- C++ compiler with C++11 support (e.g. g++ >= 4.8.2)
- ROOT (>= 5.34.30), http://root.cern.ch


Build with cmake
----------------

Checkout the code using the following command:

    git clone --recursive https://github.com/star-bnl/star-vertex.git

Compile and build libraries:

    cd star-vertex/
    mkdir build && cd build
    cmake -D CMAKE_INSTALL_PREFIX=./ ../
    make install


How to run and get results
==========================

    root4star -b -q -l 'bfc.C(1, 1000, "fzin tpcRS y2015 AgML pxlFastSim istFastSim usexgeom FieldOn
       MakeEvent Sti NoSsdIt NoSvtIt StiHftC TpcHitMover VFPPVnoCTB beamline TpxClu Idst BAna l0 Tree
       logger genvtx tpcDB bbcSim btofsim mtdSim tags emcY2 EEfs evout geantout -dstout IdTruth big
       MiniMcMk clearmem", "/scratch/smirnovd/public/star-vertex-eval/w_sim_nopileup_fzd/pythia8.starsim.ppW.100.fzd")'
