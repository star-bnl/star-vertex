#ifndef TDecayVertex_h
#define TDecayVertex_h

#include "TObject.h"


#ifdef __CINT__
enum DecayParticle_t {Undefined, Lambda, AntiLambda, Kaon};
#else
enum class DecayParticle_t {Undefined, Lambda, AntiLambda, Kaon};
#endif


class TDecayVertex : public TObject
{
public:

   TDecayVertex() : TObject() {}

   DecayParticle_t parent;
   double dEdx_d1;
   double nSigP_d1;
   double nSigPi_d1;
   double pt_d1;
   double eta_d1;
   double phi_d1;
   double dca_d1;
   int qaT_d1;
   int idT_d1;
   int tof_d1;
   double beta_d1;

   double dEdx_d2;
   double nSigP_d2;
   double nSigPi_d2;
   double pt_d2;
   double eta_d2;
   double phi_d2;
   double dca_d2;
   int qaT_d2;
   int idT_d2;
   int tof_d2;
   double beta_d2;

   double dl_p;
   double V0x_p;
   double V0y_p;
   double V0z_p;
   double im_p;
   double pt_p;
   double eta_p;
   double phi_p;

   ClassDef(TDecayVertex, 1)
};

#endif
