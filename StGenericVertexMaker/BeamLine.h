#ifndef BeamLine_h
#define BeamLine_h

#include "tables/vertexSeed.h"


namespace star_vertex {


struct BeamLine : public vertexSeed_st
{
   /// Returns x on the beamline corresponding to the passed value of z
   double X(double z) const { return x0 + dxdz*z; }

   /// Returns y on the beamline corresponding to the passed value of z
   double Y(double z) const { return y0 + dydz*z; }

   template<typename XYZ>
   XYZ GetXYZAt(double z) { return XYZ(X(z), Y(z), z); }
};


}

#endif
