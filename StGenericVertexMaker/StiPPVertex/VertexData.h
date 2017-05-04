#ifndef VertexData_h
#define VertexData_h
/**
 * full description of found vertex
 */

#include <array>

#include <TVector3.h>

#include "StGenericVertexMaker/BeamLine.h"


class VertexData {
 public:
  int id; // vertex ID assigned by PPV
  bool isTriggered; ///< Indicates whether the vertex potentially belongs to triggered event
  short mIdTruth;
  TVector3 r; // vertex position

  /// Covariance matrix is a symmetric square matrix given by six elements below
  /// the diagonal.
  ///
  /// xx
  /// xy yy
  /// xz yz zz
  ///
  std::array<double, 6>  mCovMatrix;

  int nUsedTrack; // # of tracks used to identify the vertex
  float Lmax; // maximum of the likelhood function.
  float gPtSum; // total tranverse momentum of used tracks.
  int nBtof,nCtb,nBemc,nEemc,nTpc,nAnyMatch; // number of matched tracks
  int nBtofV,nCtbV,nBemcV,nEemcV,nTpcV,nAnyVeto; // number of vetoed tracks
  
  VertexData(int vertexId, const TVector3& position);
  VertexData(int vertexId = 0) : VertexData(vertexId, TVector3(999, 999, 999)) { }
  VertexData(const TVector3& position) : VertexData(0, position) { }


  const std::array<double, 6>&  errMatrix() const { return mCovMatrix; }

  double errPerp2() const { return mCovMatrix[0] + mCovMatrix[2]; }

  /// Modifiers for vertex position and its error matrix. The beam line can be
  /// used to constrain the position and uncertainties when appropriate.
  ///@{
  void setXY(const star_vertex::BeamLine& bl);
  void setXYZ(const star_vertex::BeamLine& bl, double z, double cov_z);
  void setXYZ(std::array<double, 3> xyz, std::array<double, 6> cov_xyz);
  ///@}

  void print(ostream& os) const;
};

#endif
