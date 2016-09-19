#ifndef DecayVertexFinder_h
#define DecayVertexFinder_h

#include <vector>

#include "StarClassLibrary/StThreeVectorD.hh"
#include "StVertexRootIO/TDecayVertex.h"

class StMuDst;
class StMuTrack;


/**
 * This class tries to reconstruct secondary vertices from Lambda0 and K0
 * particles.
 */
class DecayVertexFinder
{
public:

   DecayVertexFinder() {}

   void Find(const StMuDst &muDst, std::vector<TDecayVertex>& vtxCont);

private:

   // Particle Standard Mass from PDG
   static const double M_proton;   // Proton
   static const double M_pion;     // Pion+/-

   bool Accept(StMuTrack *muTrack);

   void AddSecondaryVertex(const StMuTrack &ptra, const StThreeVectorD &ptra_p,
                           const StMuTrack &ntra, const StThreeVectorD &ntra_p,
      const StThreeVectorD &primVtx, const StThreeVectorD &secdVtx,
      std::vector<TDecayVertex>& vtxCont);
};

#endif
