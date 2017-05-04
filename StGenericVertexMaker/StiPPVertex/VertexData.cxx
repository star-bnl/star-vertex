#include <stdio.h>
#include <cmath>

#include "StGenericVertexMaker/StiPPVertex/VertexData.h"


VertexData::VertexData(int vertexId, const TVector3& position) :
  id(vertexId),
  isTriggered(false),
  mIdTruth(0),
  r(position),
  mCovMatrix{{}},
  nUsedTrack(0), Lmax(0), gPtSum(0),
  nBtof(0),  nCtb(0),  nBemc(0),  nEemc(0),  nTpc(0),  nAnyMatch(0),
  nBtofV(0), nCtbV(0), nBemcV(0), nEemcV(0), nTpcV(0), nAnyVeto(0)
{
}


void VertexData::setXY(const star_vertex::BeamLine& bl)
{
   setXYZ(bl, r.Z(), mCovMatrix[5]);
}


void VertexData::setXYZ(const star_vertex::BeamLine& bl, double z, double cov_z)
{
   r.SetXYZ( bl.X(z), bl.Y(z), z );

   // FIXME: The errors can be calculated at z exactly for the simple beam
   // line model but for now, as a first order approximation, we use the beam
   // line errors at z = 0
   mCovMatrix = std::array<double, 6>{
      bl.err_x0*bl.err_x0,
                        0, bl.err_y0*bl.err_y0,
                        0,                   0, cov_z
   };
}


void VertexData::setXYZ(std::array<double, 3> xyz, std::array<double, 6> cov_xyz)
{
   r.SetXYZ( xyz[0], xyz[1], xyz[2] );
   mCovMatrix = cov_xyz;
}


void VertexData::print(ostream& os) const {
  os << " Vertex ID="<<id<< " isTriggered: " << isTriggered << " nUsedTrack="<<nUsedTrack<<" gPtSum="<< gPtSum<<" Lmax="<< Lmax << " idTruth: " << mIdTruth
     << " match: any="<<nAnyMatch<<"-"<<nAnyVeto<<" CTB="<<nCtb<<"-"<<nCtbV<<" BEMC="<<nBemc<<"-"<<nBemcV<<" EEMC="<<nEemc<<"-"<<nEemcV<<" TPC="<<nTpc<<"-"<<nTpcV << "\n"
     << Form(" xyz: (%5.3f, %5.3f, %5.3f) +/- (%5.3f, %5.3f, %5.3f)\n", r.x(), r.y(), r.z(), std::sqrt(mCovMatrix[0]), std::sqrt(mCovMatrix[2]), std::sqrt(mCovMatrix[5]) );
}
