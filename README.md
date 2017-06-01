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

Notable features of the past and future releases.

__v3.2-rc__

* `muDst` re-vertexing for `MinuitVF`
* Seed finder based on `TSpectrum`
* Consolidate track selection in `PPV` and `MinuitVF`
* Remove dependence on StMessMgr


__v3.1-rc__

* Primary track association when re-vertexing from `muDst`
* For `muDst` additional cut on absolute track-to-vertex distance (&lt; 3cm)
to match `Sti` primary track selection
* Full covariance matrix for vertex position
* Track weight adjustment based on BTOF matching (used for vertex ranking)
* New type `BeamLine` enhancing `vertexSeed_st`


__v3.0.2__

* Vetted by LFS/UPC group
* Option to save other than default number of low-rank vertices (five) in `PPV`
* Option to exclusively select tracks pointing to BTOF detector


__v3.0.1__

* Vetted by LFS/UPC group
* In `PPV` exclude vertex seeds which cannot be fit, i.e. require at least two
tracks or one track + beam line
* For `PPV` added "real" 1D fit with beam line (identical to `MinuitVF`). Not to
be confused with the original `PPV` vertexing based on peak search in 1D
weighted track DCA histogram
* (Partially) Got rid of debug histograms in production code


__v3.0__   `2017-03-16`

* Vertex reconstruction using `muDst` trees and `PPV` vertex finder
  * Primary track are within three sigma from the vertex in transverse and
  longitudinal directions. The sigma is defined as the total uncorrelated
  uncertainty of the track and vertex
* Common fitting routines in `PPV` and `MinuitVF`
* Universal support with extensibility in mind for user-selected seeding and
fitting algorithms:
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
