#include <stdio.h>
#include <cmath>

#include "TCernLib.h"

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


/**
 * The basis is the same as in StDcaGeometry, i.e. (imp, z, psi).
 *
 * imp: impact parameter such that
 *      x =  -impact*sin(psi) = imp * cos(phi)
 *      y =   impact*cos(psi) = imp * sin(phi)
 * z  : z-coordinate of the track
 * psi: the angle of the track, e.g. for DCA psi = phi + pi/2
 */
std::array<double, 3> VertexData::positionAsDcaGeometry() const
{
   // phi: [-pi, +pi]
   // psi = r.Phi() + M_PI_2
   // if (psi > M_PI) psi -= 2*M_PI;
   double psi = r.Phi() > M_PI_2 ? r.Phi() - 3*M_PI_4 : r.Phi() + M_PI_2;

   return std::array<double, 3>{r.Perp(), r.Z(), psi};
}


/**
 * The jacobian matrix for the transformation from (imp, z, psi) to (x, y, z)
 *
 * note: imp is the impact paramer
 * note: psi = phi + pi/2
 *
 *     [  sin(psi),  0,  imp*cos(psi)]
 * J = [ -cos(psi),  0,  imp*sin(psi)]
 *     [         0,  1,           0]
 */
std::array<double, 6> VertexData::errMatrixAsDcaGeometry() const
{
   double imp = r.Perp();
   double psi = r.Phi() > M_PI_2 ? r.Phi() - 3*M_PI_4 : r.Phi() + M_PI_2;

   // The inverse of the jacobian matrix J
   double jacobian_inverse[9] =
   {
          std::sin(psi),   -std::cos(psi),   0,
                      0,                0,   1,
      std::cos(psi)/imp,     sin(psi)/imp,   0
   };

   std::array<double, 6>  newCovMatrix{0};

   TCL::trasat(jacobian_inverse, mCovMatrix.data(), newCovMatrix.data(), 3, 3);

   return newCovMatrix;
}


void VertexData::print(ostream& os) const {
  os << " Vertex ID="<<id<< " isTriggered: " << isTriggered << " nUsedTrack="<<nUsedTrack<<" gPtSum="<< gPtSum<<" Lmax="<< Lmax << " idTruth: " << mIdTruth
     << " match: any="<<nAnyMatch<<"-"<<nAnyVeto<<" CTB="<<nCtb<<"-"<<nCtbV<<" BEMC="<<nBemc<<"-"<<nBemcV<<" EEMC="<<nEemc<<"-"<<nEemcV<<" TPC="<<nTpc<<"-"<<nTpcV << "\n"
     << Form(" xyz: (%5.3f, %5.3f, %5.3f) +/- (%5.3f, %5.3f, %5.3f)\n",
             r.x(), r.y(), r.z(),
             std::sqrt(mCovMatrix[0]), std::sqrt(mCovMatrix[2]), std::sqrt(mCovMatrix[5]) );
}
