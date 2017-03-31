#include <stdio.h>
#include <cmath>
#include <St_base/StMessMgr.h>

#include "StEvent/StDcaGeometry.h"
#include "Sti/StiKalmanTrackNode.h"

#include "StGenericVertexMaker/StiPPVertex/TrackData.h"
#include "StGenericVertexMaker/StiPPVertex/VertexData.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"



template<class OriginalTrack_t>
TrackData<OriginalTrack_t>::TrackData(const OriginalTrack_t& motherTrack, const StDcaGeometry* motherDca) :
  vertexID(0),
  mother(motherTrack),
  dca(motherDca),
  mIdTruth(0),
  mQuality(0),
  mIdParentVx(0),
  dcaTrack(), zDca(0), ezDca(0), rxyDca(0),
  gPt(0),
  mBtof(0), mCtb(0), mBemc(0), mEemc(0), mTpc(0),
  anyMatch(false), anyVeto(false),
  weight(1),
  btofBin(-1), ctbBin(-1), bemcBin(-1), eemcBin(-1)
{ }


/** Specialized constructor to create tracks from StMuTrack-s */
template<>
TrackData<StMuTrack>::TrackData(const StMuTrack &motherTrack, const StDcaGeometry* trackDca) :
  vertexID(0),
  mother(motherTrack),
  dca(trackDca),
  mIdTruth(0),
  mQuality(0),
  mIdParentVx(0),
  dcaTrack(), zDca(0), ezDca(0), rxyDca(0),
  gPt(0),
  mBtof(0), mCtb(0), mBemc(0), mEemc(0), mTpc(0),
  anyMatch(false), anyVeto(false),
  weight(1),
  btofBin(-1), ctbBin(-1), bemcBin(-1), eemcBin(-1)
{
  mIdTruth    = motherTrack.idTruth();
  mQuality    = motherTrack.qaTruth();
  mIdParentVx = motherTrack.idParentVx();

  if (trackDca) {
    zDca   = trackDca->z();
    ezDca  = std::sqrt(trackDca->errMatrix()[2]);
    rxyDca = trackDca->impact();
    gPt    = trackDca->pt();
  }
}


template<>
TrackData<StiKalmanTrack>::TrackData(const StiKalmanTrack &motherTrack) :
  TrackData(motherTrack, new StDcaGeometry())
{
  // This code is adopted from StiStEventFiller::fillDca()
  StiKalmanTrackNode *tNode = const_cast<StiKalmanTrack&>(motherTrack).extrapolateToBeam();

  if (tNode)
  {
    const StiNodePars &pars = tNode->fitPars();
    const StiNodeErrs &errs = tNode->fitErrs();
    float alfa = tNode->getAlpha();
    float setp[7] = {(float)pars.y(),    (float)pars.z(),    (float)pars.phi(),
                     (float)pars.ptin(), (float)pars.tanl(), (float)pars.curv(), (float)pars.hz()};
    setp[2] += alfa;
    float sete[15];

    for (int i=1, li=1, jj=0; i<kNPars; li += ++i) {
      for (int j=1;j<=i;j++) {
         sete[jj++] = errs.G()[li+j];
      }
    }

    const_cast<StDcaGeometry*>(dca)->set(setp, sete);
  }
  else
  {
    delete dca;
    dca = nullptr;
  }

  // The following code is moved from the original StPPVertexFinder::examinTrackDca()
  StiKalmanTrackNode* bmNode = motherTrack.getInnerMostNode();

  if ( bmNode && bmNode->isDca() )
  {
    zDca   = bmNode->getZ(); // FIXME: Why not bmNode->z_g() instead?
    ezDca  = std::sqrt(bmNode->getCzz());
    rxyDca = std::sqrt(bmNode->x_g()*bmNode->x_g() + bmNode->y_g()*bmNode->y_g());
    gPt    = bmNode->getPt();
  }
}


template<>
TrackData<StiKalmanTrack>::~TrackData()
{
  delete dca;
  dca = nullptr;
}


//==========================================================
//==========================================================
template<class OriginalTrack_t>
bool TrackData<OriginalTrack_t>::matchVertex(VertexData &V, float dzMax)
{
  float dz = zDca - V.r.z();
  bool ret = fabs(dz) < dzMax + ezDca;

  if (ret)
     LOG_DEBUG<< Form("PPV::matchTr2Ver VerID=%d  weight=%.2f anyM=%d anyV=%d  m: ctb=%d  bemc=%d eemc=%d tpc=%d dz=%.2f +/- %.2f\n",V.id,weight,anyMatch,anyVeto,mCtb,mBemc,mEemc,mTpc,dz,ezDca)<<endm;

  return ret;
}


template<class OriginalTrack_t>
double TrackData<OriginalTrack_t>::calcChi2DCA(const VertexData &V) const
{
   double err2;
   double vxyz[3];
   V.r.GetXYZ(vxyz);

   double dist = dca->thelix().Dca(vxyz, &err2);
   double chi2 = dist*dist/err2;

   return chi2;
}


//==========================================================
//==========================================================
template<class OriginalTrack_t>
void TrackData<OriginalTrack_t>::scanNodes(std::vector<int> &hit, int jz0)
{
  /* INPUT: vector of hits for active nodes
     i=[0,jz0-1] is on one side of z-Axis
     i=[jz0,mx-1] on another side
  */

  // params
  const int minCenter=4 , minMiss=6; // criteria for Match & Veto
  const int mxDev=2; // max # of deviations from expected pattern

  // printf("patt size=%d, jz0=%d\n",hit.size(),jz0);

  int nPatt[2]{};
  int i;
  bool vetoL=false, vetoR=false, matchL=false, matchR=false;

  // Scan Left end for missing hits
  for(i=0;i<(int)hit.size();i++) {
    // printf("i=%d hit=%d\n",i,hit[i]);
    nPatt[hit[i]]++;
    if(nPatt[1]>mxDev) break;
    if(nPatt[0]<minMiss) continue;
    vetoL=true;
    break;
  }

  // printf("vetoL=%d   nUp=%d nDwn=%d\n\n",vetoL, nPatt[1],nPatt[0]);

  // Scan Right end for missing hits
  memset(nPatt,0,sizeof(nPatt));
  for(i=hit.size()-1; i>=0;i--) {
    // printf("i=%d hit=%d\n",i,hit[i]);
    nPatt[hit[i]]++;
    if(nPatt[1]>mxDev) break;
    if(nPatt[0]<minMiss) continue;
    vetoR=true;
    break;
  }

  // printf("vetoR=%d nUp=%d nDwn=%d\n\n",vetoR, nPatt[1],nPatt[0]);

  if(jz0>minCenter && jz0<(int)hit.size()-minCenter) { // examin membrane
    // Scan Left half at membrane
    memset(nPatt,0,sizeof(nPatt));
    for(i=jz0-1; i>=0;i--) {
      //  printf("i=%d hit=%d\n",i,hit[i]);
      nPatt[hit[i]]++;
      if(nPatt[0]>mxDev) break;
      if(nPatt[1]<minCenter) continue;
      matchL=true;
      break;
    }
    
    // printf("matchL=%d   nUp=%d nDwn=%d\n\n",matchL, nPatt[1],nPatt[0]);
    
    
    // Scan Right half at membrane
    memset(nPatt,0,sizeof(nPatt));
    for(i=jz0;i<(int)hit.size();i++) {
      // printf("i=%d hit=%d\n",i,hit[i]);
      nPatt[hit[i]]++;
      if(nPatt[0]>mxDev) break;
      if(nPatt[1]<minCenter) continue;
      matchR=true;
      break;
    }
    
    // printf("matchR=%d   nUp=%d nDwn=%d\n\n",matchR, nPatt[1],nPatt[0]);
    
  } // end of membrane test

  bool match=  matchL && matchR;
  bool veto=(vetoL || vetoR) && !match;

  //  printf(" TPC Conclusion: tpcMatch=%d  tpcVeto=%d \n\n",match,veto);
  updateAnyMatch(match,veto,mTpc);
  weight*=getTpcWeight();
}

//==========================================================
//==========================================================
template<class OriginalTrack_t>
void TrackData<OriginalTrack_t>::updateAnyMatch(bool match, bool veto, int & mXXX)
{
  if(match) {
    anyMatch=true;
    anyVeto=false;
    mXXX=1;
    //  } else if(veto && (!anyMatch) ) {
  } else if(veto && (!match) ) {
    anyVeto=true;
    mXXX=-1;
  } else {
    mXXX=0;
  }
}

//==========================================================
//==========================================================
template<class OriginalTrack_t>
float TrackData<OriginalTrack_t>::getTpcWeight()
{
  const float Wdunno=1, Wmatch=5, Wveto=0.2;
  if(mTpc>0) return Wmatch;  
  if(mTpc<0) return Wveto;
  return Wdunno;
}


template<class OriginalTrack_t>
void TrackData<OriginalTrack_t>::print(ostream& os) const
{
   os << Form("vertID=%d track@z0=%.2f +/- %.2f gPt=%.3f, rxyDca: %.3f, idTruth: %d, qaTruth: %d, idParentVx: %d",
	      vertexID, zDca, ezDca, gPt, rxyDca, mIdTruth, mQuality, mIdParentVx);

   if (dca) {
      os << Form(" dca: (%5.3f, %5.3f, %5.3f) +/- (imp: %5.3f, %5.3f)",
            dca->origin().x(), dca->origin().y(), dca->origin().z(),
            std::sqrt( dca->errMatrix()[0] ), std::sqrt( dca->errMatrix()[2] ) );
   }

   os << endl;
}


template class TrackData<StiKalmanTrack>;
template class TrackData<StMuTrack>;
