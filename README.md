How to build and install the libraries
======================================


Prerequisites
-------------

- C++ compiler with C++11 support (e.g. g++ >= 4.8.2)
- ROOT (>= 5.34.30), http://root.cern.ch
- cmake (>= 2.8)


Build with cmake
----------------

Checkout the code using the following command:

    $ git clone --recursive https://github.com/star-bnl/star-vertex.git

Compile and build libraries:

    $ cd star-vertex/
    $ mkdir build && cd build
    $ cmake -D CMAKE_INSTALL_PREFIX=./ ../
    $ make install


How to run and get results
==========================

Starting v2.0, one can use the `beamline` and `beamline3D` options (or no beam
line option) in addition to the preferred vertex finder such as `MinuitVF` or
`VFPPVnoCTB`. For example, within the STAR environment the command may look
as:

    root4star -b -q -l 'bfc.C(1, 1000, "fzin tpcRS y2015 AgML usexgeom FieldOn
       MakeEvent Sti NoSsdIt NoSvtIt TpcHitMover VFPPVnoCTB beamline3D TpxClu
       tpcDB btofsim IdTruth clearmem", "/path/to/fzd/pythia8.starsim.ppW.100.fzd")'


Release History
===============

Highlighted features of the past and future releases.

__v3.x__

* Add primary track association when reconstruct from `muDst`
  * Support for `MinuitVF`
* New seed finder based on `TSpectrum`
* Get rid of debugging histograms in production code


__v3.0-rc__

* Vertex reconstruction using `muDst` trees and `PPV` vertex finder
* Common approach to fitting in `PPV` and `MinuitVF`
* Consolidate track selection in `PPV` and `MinuitVF`
* Unified selection of seeding and fitting algorithms:
  * `SeedFinder: MinuitVF, PPVLikelihood, TSpectrum`
  * `VertexFit: NoBeamline, Beamline1D, Beamline3D`


__v2.0__   `2016-12-12`

* Support 3D vertex fits in `MinuitVF` and `PPV`
* Optional vertex position constraint by the beam line
  * Routine to calculate `chi^2` for the vertex and the beam line
  * Proper calculation of uncertainties


__v1.0__   `2016-03-11`

* Original code inherited from previous developers
* Build with `cmake`, continuous integration
