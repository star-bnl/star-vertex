#include <cmath>

#include "StSecondaryVertexMaker/DecayVertexFinder.h"

#include "StarClassLibrary/SystemOfUnits.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"


const double DecayVertexFinder::M_proton = 0.938272;
const double DecayVertexFinder::M_pion   = 0.139570;


void DecayVertexFinder::Find(const StMuDst &muDst, std::vector<TDecayVertex>& vtxCont)
{
   vtxCont.clear();

   // Get a pointer to the class holding event-wise information
   StMuEvent &muEvent = *muDst.event();

   if ( muDst.numberOfGlobalTracks() < 2 )
   {
      std::cout << "DecayVertexFinder: not enough tracks" << std::endl;
      return;
   }

   int nTracks = muDst.numberOfGlobalTracks();

   // Loop over all tracks in the event
   for ( int p = 0; p < nTracks; ++p )
   {
      StMuTrack *ptra = muDst.globalTracks(p);

      // Select positive tracks
      if ( ptra->charge() < 0 ) continue;
      if ( !Accept( ptra ) ) continue;

      double ptra_nSigp  = ptra->nSigmaProtonFit();
      double ptra_nSigpi = ptra->nSigmaPionFit();

      // Do not consider this track ptra if it is not consistent with either proton or pion
      if ( std::fabs(ptra_nSigp) > 3.0 && std::fabs(ptra_nSigpi) > 3.0 ) continue;

      // Another loop over all tracks in the event to build secondary vertex
      // from track pairs
      for ( int n = 0; n < nTracks; ++n )
      {
         if ( n == p ) continue;

         StMuTrack *ntra = muDst.globalTracks(n);

         // Select negative tracks
         if ( ntra->charge() > 0 ) continue;
         if ( !Accept( ntra ) ) continue;

         double ntra_nSigp  = ntra->nSigmaProtonFit();
         double ntra_nSigpi = ntra->nSigmaPionFit();

         // Do not consider this track ntra if it is not consistent with either proton or pion
         if (std::fabs(ntra_nSigp)  > 3 && std::fabs(ntra_nSigpi) > 3) continue;

         // Do not consider this pair of tracks if they are not consistent with both being pions
         if (std::fabs(ptra_nSigpi) > 3 && std::fabs(ntra_nSigpi) > 3) continue;

         StPhysicalHelixD ptra_helix = ptra->helix();
         StPhysicalHelixD ntra_helix = ntra->helix();

         // Calculate the distance between the two tracks/helices
         pair<double, double> ss = ptra_helix.pathLengths(ntra_helix);
         double s1 = ss.first;
         double s2 = ss.second;

         StThreeVectorD tdca = ptra_helix.at(s1) - ntra_helix.at(s2);
         double tdca2 = tdca.mag();  // closest distence between two daughters particle

         if ( tdca2 > 2.0 ) continue;

         StThreeVectorD secondary_vtx = ( ptra_helix.at(s1) + ntra_helix.at(s2) ) / 2.0;
         StThreeVectorD primary_vtx   = muEvent.primaryVertexPosition();
         StThreeVectorD decay_length  = secondary_vtx - primary_vtx;  // V0_distance to Primary Vetex

         if ( decay_length.mag() < 0.2 ) continue;

         StThreeVectorD ptra_p = ptra_helix.momentumAt(s1, muEvent.magneticField()*kilogauss);
         StThreeVectorD ntra_p = ntra_helix.momentumAt(s2, muEvent.magneticField()*kilogauss);
         StThreeVectorD tV0_p = ptra_p + ntra_p;

         if ( tV0_p.mag() < 1E-5 ) continue;

         double dot = tV0_p.dot(decay_length);

         // Calculate the DCA of decaying particle to the secondary vertex
         double tdcaV0 = sqrt(decay_length.mag2() - dot*dot / tV0_p.mag2());

         if ( tdcaV0 > 1.5 ) continue;

         // Cosine?
         double tcrp = dot / (tV0_p.mag() * decay_length.mag());

         if ( tcrp < 0.9 ) continue;

         AddSecondaryVertex(*ptra, ptra_p, *ntra, ntra_p, primary_vtx, secondary_vtx, vtxCont);
      }
   }
}


bool DecayVertexFinder::Accept(StMuTrack *stMuTrack)
{
   if ( stMuTrack->nHitsFit()  < 15 ) return false;
   if ( stMuTrack->nHitsPoss() < 30 ) return false;

   if ( double(stMuTrack->nHitsFit()) / double(stMuTrack->nHitsPoss()) < 0.51 ) return false;

   if ( stMuTrack->flag() < 0 || stMuTrack->flag() > 1000 ) return false;

   return true;
}


void DecayVertexFinder::AddSecondaryVertex(const StMuTrack &ptra, const StThreeVectorD &ptra_p,
      const StMuTrack &ntra, const StThreeVectorD &ntra_p,
      const StThreeVectorD &primary_vtx, const StThreeVectorD &secondary_vtx,
      std::vector<TDecayVertex>& vtxCont)
{
   double ptra_nSigp  = ptra.nSigmaProtonFit();
   double ptra_nSigpi = ptra.nSigmaPionFit();

   double ntra_nSigp  = ntra.nSigmaProtonFit();
   double ntra_nSigpi = ntra.nSigmaPionFit();

   StThreeVectorD decay_length = secondary_vtx - primary_vtx;

   double Mp = 0;
   double Mn = 0;

   if      ( std::fabs(ptra_nSigp)  < std::fabs(ptra_nSigpi) && std::fabs(ptra_nSigp)  < 3) Mp = M_proton;
   else if ( std::fabs(ptra_nSigpi) < std::fabs(ptra_nSigp)  && std::fabs(ptra_nSigpi) < 3) Mp = M_pion;

   if      ( std::fabs(ntra_nSigp)  < std::fabs(ntra_nSigpi) && std::fabs(ntra_nSigp)  < 3) Mn = M_proton;
   else if ( std::fabs(ntra_nSigpi) < std::fabs(ntra_nSigp)  && std::fabs(ntra_nSigpi) < 3) Mn = M_pion;

   // Do not add vertex if either of the daughters cannot be identified as a proton or a pion
   // OR when both of the daughers are protons
   if (!Mp || !Mn || (Mp == M_proton && Mn == M_proton) ) {
      return;
   }

   double Ep = sqrt( Mp*Mp + ptra_p.mag2() );
   double En = sqrt( Mn*Mn + ntra_p.mag2() );

   StLorentzVectorD tp_lv(Ep, ptra_p);
   StLorentzVectorD tn_lv(En, ntra_p);
   StLorentzVectorD tV0_lv = tp_lv + tn_lv;

   DecayParticle_t parent = DecayParticle_t::Undefined;

   if ( tV0_lv.m() >= 0.42 && tV0_lv.m() < 0.58 && Mp == M_pion && Mn == M_pion )
      parent = DecayParticle_t::Kaon;
   else if ( tV0_lv.m() >= 1.08 && tV0_lv.m() < 1.16 &&
            ( (Mp == M_proton && Mn == M_pion) || (Mp == M_pion && Mn == M_proton) )
           )
      parent = DecayParticle_t::Lambda;
   else {
      // Cannot identify decaying particle
      return;
   }

   TDecayVertex vtx;

   vtx.parent    = parent;

   vtx.dEdx_d1   = ptra.dEdx()*1E+6;
   vtx.nSigP_d1  = ptra_nSigp;
   vtx.nSigPi_d1 = ptra_nSigpi;
   vtx.pt_d1     = tp_lv.perp();
   vtx.eta_d1    = tp_lv.pseudoRapidity();
   vtx.phi_d1    = tp_lv.phi();
   vtx.dca_d1    = ptra.dca().mag();
   vtx.qaT_d1    = ptra.qaTruth();
   vtx.idT_d1    = ptra.idTruth();
   vtx.tof_d1    = ptra.btofPidTraits().matchFlag();
   vtx.beta_d1   = ptra.btofPidTraits().beta();

   vtx.dEdx_d2   = ntra.dEdx()*1E+6;
   vtx.nSigP_d2  = ntra_nSigp;
   vtx.nSigPi_d2 = ntra_nSigpi;
   vtx.pt_d2     = tn_lv.perp();
   vtx.eta_d2    = tn_lv.pseudoRapidity();
   vtx.phi_d2    = tn_lv.phi();
   vtx.dca_d2    = ntra.dca().mag();
   vtx.qaT_d2    = ntra.qaTruth();
   vtx.idT_d2    = ntra.idTruth();
   vtx.tof_d2    = ntra.btofPidTraits().matchFlag();
   vtx.beta_d2   = ntra.btofPidTraits().beta();

   vtx.dl_p      = decay_length.mag();
   vtx.V0x_p     = secondary_vtx.x();
   vtx.V0y_p     = secondary_vtx.y();
   vtx.V0z_p     = secondary_vtx.z();
   vtx.im_p      = tV0_lv.m();
   vtx.pt_p      = tV0_lv.perp();
   vtx.eta_p     = tV0_lv.pseudoRapidity();
   vtx.phi_p     = tV0_lv.phi();

   vtxCont.push_back(vtx);
}
