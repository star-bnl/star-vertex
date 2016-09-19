#ifndef TPrimaryVertex_h
#define TPrimaryVertex_h

#include "TObject.h"

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"


class TPrimaryVertex : public TObject
{
public:

   TPrimaryVertex() : TObject() { Init(); }

   int rank, mult, refMult;
   bool maxmult;
   float primX, primY, primZ, zVpd;
   StThreeVectorF positionError;
   float McX, McY, McZ;
   float chi2;
   int beam, postx, prompt, cross, tof, notof, EEMC, noEEMC, BEMC, noBEMC;


   void Init() {
      rank = mult = refMult = -999;
      maxmult = false;
      primX = primY = primZ = zVpd = 999.f;
      positionError.set(999.f, 999.f, 999.f);
      McX = McY = McZ = chi2 = 999.f;
   }


   void Set(const StMuPrimaryVertex &recoVertex, const StMuMcVertex* mcVertex=nullptr, bool isMaxMult=false, float z_Vpd=999.)
   {
      rank    = recoVertex.ranking();
      mult    = recoVertex.noTracks();
      refMult = recoVertex.refMult();
      maxmult = isMaxMult;
      primX   = recoVertex.position().x();
      primY   = recoVertex.position().y();
      primZ   = recoVertex.position().z();
      zVpd    = z_Vpd;
      positionError = recoVertex.posError();
      McX     = mcVertex ? mcVertex->XyzV().x() : 999.f;
      McY     = mcVertex ? mcVertex->XyzV().y() : 999.f;
      McZ     = mcVertex ? mcVertex->XyzV().z() : 999.f;
      chi2    = recoVertex.chiSquared();
      beam    = recoVertex.isBeamConstrained() ? 1 : 0;
      postx   = recoVertex.nPostXtracks();
      prompt  = recoVertex.nPromptTracks();
      cross   = recoVertex.nCrossCentralMembrane();
      tof     = recoVertex.nCTBMatch()    + recoVertex.nBTOFMatch();
      notof   = recoVertex.nCTBNotMatch() + recoVertex.nBTOFNotMatch();
      EEMC    = recoVertex.nEEMCMatch();
      noEEMC  = recoVertex.nEEMCNotMatch();
      BEMC    = recoVertex.nBEMCMatch();
      noBEMC  = recoVertex.nBEMCNotMatch();
   }

   ClassDef(TPrimaryVertex, 1)
};

#endif
