/***************************************************************************
 *
 * $Id: StMinuitVertexFinder.cxx,v 1.8 2006/05/31 04:09:52 fisyak Exp $
 *
 * Author: Thomas Ullrich, Feb 2002
 ***************************************************************************
 *
 * Description: 
 *
 ***************************************************************************
 *
 * $Log: StMinuitVertexFinder.cxx,v $
 * Revision 1.8  2006/05/31 04:09:52  fisyak
 * Use dca track parameters for primary vertex fit
 *
 * Revision 1.7  2006/05/09 17:51:05  mvl
 * Added protection against event->emcCollection()==0
 *
 * Revision 1.6  2006/05/04 20:01:31  jeromel
 * Switched to logger
 *
 * Revision 1.5  2006/04/26 15:37:04  jeromel
 * mVertexOrderMethod (To be tested)
 *
 * Revision 1.4  2006/04/25 13:06:44  mvl
 * Seed-finding range extended to -200<vtx_z<200
 *
 * Revision 1.3  2006/04/08 23:21:15  fisyak
 * Add protection for  bemcDet==0
 *
 * Revision 1.2  2006/04/08 19:06:29  mvl
 * Update for multiple vertex finding and rank calculation for identifying the
 * triggered vertex. Ranks are based on mean dip angle of tracks, BEMC matches
 * and tracks crossing the central membrane and optimised for Cu+Cu.
 * The track cuts are now bit tighter (dca<2 in transverse direction and
 * nfitpoints > 15) to suppress 'fake' vertices.
 * In addition, a lower multiplicity cut of 5 tracks is implemented.
 *
 * Revision 1.18  2005/12/08 16:54:13  fisyak
 * Fix mRequireCTB, one more non initialized variable
 *
 * Revision 1.17  2005/09/29 22:17:30  fisyak
 * more strict cut for failed vertex
 *
 * Revision 1.16  2005/07/19 21:52:45  perev
 * MultiVertex
 *
 * Revision 1.15  2005/06/21 02:16:36  balewski
 * multiple prim vertices are stored in StEvent
 *
 * Revision 1.14  2004/12/13 20:39:58  fisyak
 * Add initaition of StGenericVertexFinder variables, replace mDumMaker by StMaker::GetChain() method
 *
 * Revision 1.13  2004/08/17 20:41:18  perev
 * LeakOff
 *
 * Revision 1.12  2004/08/04 21:57:56  balewski
 * toward smarter ppLMV5
 *
 * Revision 1.11  2004/07/23 01:28:55  jeromel
 * Typo corrected
 *
 * Revision 1.10  2004/07/23 00:59:10  jeromel
 * Removed methods (moved in base class). Changed setFlagBase().
 *
 * Revision 1.9  2004/04/06 02:43:43  lbarnby
 * Fixed identification of bad seeds (no z~0 problem now). Better flagging. Message manager used.
 *
 * Revision 1.8  2004/04/04 23:20:13  jeromel
 * isfinite() -> finite()
 *
 * Revision 1.7  2004/03/23 16:15:04  lbarnby
 * Extra protection for non-finite track length. User function to not use ITTF tracks
 *
 * Revision 1.6  2003/10/09 16:40:12  perev
 * delete helix object added
 *
 * Revision 1.5  2003/10/06 04:37:58  perev
 * delete helix
 *
 * Revision 1.4  2003/09/02 17:58:19  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.3  2003/05/12 21:10:06  lbarnby
 * Made destructor virtual
 *
 * Revision 1.2  2003/05/09 22:20:00  lbarnby
 * Now also calculates and reports error on vertex. Corrected filter to use ITTF tracks. Some temporary protections against inf/Nan. Skip delete of TMinuit class since causing seg. fault.
 *
 * Revision 1.1  2002/12/05 23:42:46  hardtke
 * Initial Version for development and integration
 *
 **************************************************************************/
#include "StMinuitVertexFinder.h"
#include "StEventTypes.h"
#include "StEnumerations.h"
#include "TMinuit.h"
#include "StGlobals.hh"
#include "SystemOfUnits.h"
#include "StCtbMatcher.h"
#include "StMessMgr.h"
#include <math.h>
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StDcaGeometry.h"
vector<StDcaGeometry*>   StMinuitVertexFinder::mDCAs;
vector<StPhysicalHelixD> StMinuitVertexFinder::mHelices;
vector<unsigned short>   StMinuitVertexFinder::mHelixFlags;
vector<double>           StMinuitVertexFinder::mSigma;
vector<double>           StMinuitVertexFinder::mZImpact;
double                   StMinuitVertexFinder::mWidthScale = 1;
double                   StMinuitVertexFinder::mX0;
double                   StMinuitVertexFinder::mY0;
double                   StMinuitVertexFinder::mdxdz;
double                   StMinuitVertexFinder::mdydz;
bool                     StMinuitVertexFinder::requireCTB;
int                      StMinuitVertexFinder::nCTBHits;
bool                     StMinuitVertexFinder::mUseDCA;
//==========================================================
//==========================================================
void StMinuitVertexFinder::Clear(){
  StGenericVertexFinder::Clear();
  mStatusMin    = 0;
  mNSeed = 0;
}

void
StMinuitVertexFinder::setExternalSeed(const StThreeVectorD& s)
{
    mExternalSeedPresent = true;
    mExternalSeed = s;
}


StMinuitVertexFinder::StMinuitVertexFinder() {
  LOG_INFO << "StMinuitVertexFinder::StMinuitVertexFinder is in use." << endm;
  mBeamHelix =0;

  mMinNumberOfFitPointsOnTrack = 15; 
  mDcaZMax = 3;     // Note: best to use integer numbers
  mMinTrack = 5;
  mRImpactMax = 2;

  mMinuit = new TMinuit(3);         
  mMinuit->SetFCN(&StMinuitVertexFinder::fcn);
  mMinuit->SetPrintLevel(-1);
  mMinuit->SetMaxIterations(1000);
  mExternalSeedPresent = false;
  mVertexConstrain = false;
  mRequireCTB = false;
  requireCTB = false;
  mUseITTF   = false;
  mUseDCA    = false;
  mVertexOrderMethod = orderByRanking; // change ordering by ranking
}
 

 StMinuitVertexFinder::~StMinuitVertexFinder()
 {
   delete mBeamHelix; mBeamHelix=0;
   LOG_WARN << "Skipping delete Minuit in StMinuitVertexFinder::~StMinuitVertexFinder()" << endm;
   //delete mMinuit;
   mDCAs.clear();
   mHelices.clear();
   mHelixFlags.clear();
   mZImpact.clear();
   mSigma.clear();
 }


void
StMinuitVertexFinder::setFlagBase(){
  if(mUseITTF){
    mFlagBase = 8000;
  } else {
    mFlagBase = 1000;
  }
}

int StMinuitVertexFinder::findSeeds() {
  mNSeed = 0;

  int zImpactArr[400]; // simple array to 'histogram' zImpacts
  for (int i=0; i < 400; i++)
    zImpactArr[i]=0;

  Int_t nTrk = mZImpact.size();
  for (int iTrk=0; iTrk < nTrk; iTrk++) {
    if (fabs(mZImpact[iTrk]) < 200)
      zImpactArr[int(mZImpact[iTrk]+200)]++;
  }

  // Search for maxima using sliding 3-bin window
  Int_t nOldBin = 0;
  int slope = 0;
  int nBinZ = 3;
  for (int iBin=0; iBin < 400 - nBinZ; iBin++) {
    int nTrkBin = 0;
    for (int iBin2=0; iBin2 < nBinZ; iBin2++) {
      nTrkBin += zImpactArr[iBin + iBin2];
    }
    if (nTrkBin > nOldBin)
      slope = 1;
    else if (nTrkBin < nOldBin) {
      if (slope == 1) {
	if (mNSeed < maxSeed) {
	  float seed_z = -200 + iBin + (Float_t)nBinZ / 2 - 1;
	  Double_t meanZ = 0;
	  int nTrkZ = 0;
	  for (int iTrk = 0; iTrk < nTrk; iTrk ++ ) {
	    if ( fabs(mZImpact[iTrk] - seed_z ) < mDcaZMax ) {
	      meanZ += mZImpact[iTrk];
	      nTrkZ++;
	    }
	  }
	  if (nTrkZ > mMinTrack) {
	    if (mDebugLevel) 
	      cout << "Seed " << mNSeed << ", z " << seed_z << " nTrk " << nTrkZ << " meanZ/nTrkZ " << meanZ/nTrkZ << endl;
	    seed_z = meanZ/nTrkZ;
	    mSeedZ[mNSeed] = seed_z;
	    mNSeed ++;
	  }
	}
	else {
	  LOG_WARN << "Too many vertex seeds, limit=" << maxSeed << endm;
	}	
      }
      slope = -1;
    }
    if (mDebugLevel > 1) 
      cout << "iBin " << iBin << " nTrkBin " << nTrkBin << " nOldBin " << nOldBin << ", slope " << slope << " mNSeed " << mNSeed << endl; 
    nOldBin = nTrkBin;
  }

  LOG_INFO << "Found " << mNSeed << " seeds" << endm;
  return mNSeed;
}

void StMinuitVertexFinder::fillBemcHits(StEvent *event){
  static int nMod = 120;
  static int nEta = 20;
  static int nSub = 2;
  static float mEmcThres = 0.15;
  for (int m=0; m < nMod; m++) 
    for (int e=0; e < nEta; e++) 
      for (int s=0; s < nSub; s++) 
	mBemcHit[m][e][s]=0;

  int n_emc_hit=0;
  if (event->emcCollection() && event->emcCollection()->detector(kBarrelEmcTowerId)) {
    StEmcDetector* bemcDet = event->emcCollection()->detector(kBarrelEmcTowerId);
    for (int iMod=0; iMod < nMod; iMod++) {
      if (!bemcDet->module(iMod)) 
        continue;
      StEmcModule *mod = bemcDet->module(iMod);
      StSPtrVecEmcRawHit&  hits = mod->hits();
      for (StEmcRawHitIterator hitIter = hits.begin(); hitIter != hits.end(); hitIter++) {
        StEmcRawHit *hit = *hitIter;
        if (hit->energy() > mEmcThres) {
          mBemcHit[hit->module()-1][hit->eta()-1][hit->sub()-1]=1;
          n_emc_hit++; 
        }
      }
    }
  }
  if (mDebugLevel) 
    cout << "Found " << n_emc_hit << " emc hits" << endl;
}

int  
StMinuitVertexFinder::matchTrack2BEMC(const StTrack *track){
  static const double rBemc = 242; // middle of tower
  static StEmcGeom *bemcGeom = StEmcGeom::getEmcGeom("bemc");
  //static double rBemc = bemcGeom->Radius(); // front face??
  //cout << "rBemc: " << rBemc << endl;

  if (track->outerGeometry()==0) {
    if (mDebugLevel) // Happens only rarely
      cout << "No outer track geom" << endl;
    return 0;
  }

  StPhysicalHelixD helix = track->outerGeometry()->helix();

  if (!helix.valid()) {
    if (mDebugLevel) // Happens only rarely
      cout << "Invalid helix" << endl;
    return 0;
  }
     
  pairD  d2;
  d2 = helix.pathLength(rBemc);
  if (d2.first > 99999999 && d2.second > 99999999)
    return 0;
  double path=d2.second;
  if (d2.first >= 0 && d2.first < d2.second)
    path = d2.first;
  else if(d2.first>=0 || d2.second<=0) {
    LOG_WARN << "Unexpected solution for track crossing Cyl:" 
			<< " track mom = " 
			<< track->geometry()->momentum().mag() 
			<< ", d2.first=" << d2.first 
			<< ", second=" << d2.second <<", trying first" << endm;
    path=d2.first;
  }
  StThreeVectorD posCyl = helix.at(path);

  float phi=posCyl.phi();
  float eta=posCyl.pseudoRapidity();

  if (fabs(eta) < 1) {
    int mod, ieta, sub;
    if (bemcGeom->getBin(phi, eta, mod, ieta, sub) == 0) {
      // There is some edge effect leading to sub=-1. 
      // Strange, but leave in for now// HOW CAN IT be: sub == -1????
      if (sub > 0 && mBemcHit[mod-1][ieta-1][sub-1]) {
	return 1;
      }
    }
  }
  return 0;
}

int StMinuitVertexFinder::checkCrossMembrane(const StTrack *track) {
  int n_pnt_neg=0, n_pnt_pos=0;
  
  StPtrVecHit hits = track->detectorInfo()->hits(kTpcId);
  for (StHitIterator hitIter = hits.begin(); hitIter != hits.end(); hitIter++) {
    if ((*hitIter)->position().z() < 0)
      n_pnt_neg++;
    if ((*hitIter)->position().z() > 0)
      n_pnt_pos++;
  }
  return (n_pnt_pos > 5 && n_pnt_neg > 5);
}

void StMinuitVertexFinder::calculateRanks() {    
  
  // Calculation of veretx ranks to select 'best' (i.e. triggered)
  // vertex
  // Three ranks are used, based on average dip, number of BEMC matches 
  // and tarcks crossing central membrane
  // Each rank is normalised to be 0 for 'average' valid vertices and
  // has a sigma of 1. 
  //
  // Values are limited to [-5,1] for each rank
  //
  // Note that the average dip angle ranking is naturally peaked to 1,
  // while the others are peaked at 1 (intentionally) due to rounding.

  // A fancier way would be to calculate something more like 
  // a likelihood based on the expected distributions. 
  // That's left as an excercise to the reader.


  int nBemcMatchTot = 0;
  int nVtxTrackTot = 0;
  for (int iVertex=0; iVertex < size(); iVertex++) {
    StPrimaryVertex *primV = getVertex(iVertex);
    nVtxTrackTot += primV->numTracksUsedInFinder();
    nBemcMatchTot += primV->numMatchesWithBEMC();
  }

  mBestRank = -999;
  mBestVtx  = 0;
  for (int iVertex=0; iVertex < size(); iVertex++) {
    StPrimaryVertex *primV = getVertex(iVertex);
    // expected values based on Cu+Cu data
    float avg_dip_expected = -0.0033*primV->position().z();
    float n_bemc_expected = 0;
    if (nVtxTrackTot) 
      n_bemc_expected = (1-0.25*(1-(float)primV->numTracksUsedInFinder()/nVtxTrackTot))*nBemcMatchTot; 

    float n_cross_expected = abs(primV->position().z())*0.0020*primV->numTracksUsedInFinder(); // old coeff 0.0016 with dca 3 and 10 points on track

    if (mDebugLevel)
      cout << "vertex z " << primV->position().z() << " dip expected " << avg_dip_expected << " bemc " << n_bemc_expected << " cross " << n_cross_expected << endl;
    float rank_avg_dip = 1 - fabs(primV->meanDip() - avg_dip_expected)*sqrt((float)primV->numTracksUsedInFinder())/0.67;  // Sigma was 0.8 for old cuts
    if (rank_avg_dip < -5)
      rank_avg_dip = -5;

    float rank_bemc = 0;
    if (n_bemc_expected >= 1) { 
      //float sigma = 0.12*n_bemc_match_tot;
      float sigma = 0.5*sqrt(n_bemc_expected);
      if ( sigma < 0.75 ) { // limit sigma to avoid large weights 
	// at small multiplicity
	sigma = 0.75;
      }
      rank_bemc = (primV->numMatchesWithBEMC() - n_bemc_expected)/sigma+0.5; // distribution is asymmetric; add 0.5 
    }
    if (rank_bemc < -5)
      rank_bemc = -5;
    if (rank_bemc > 1)
      rank_bemc = 1;
    
    float rank_cross = 0;
    if ( n_cross_expected >= 1 ) {
      float sigma=1.1*sqrt(n_cross_expected);
      rank_cross = (primV->numTracksCrossingCentralMembrane() - n_cross_expected)/sigma;
    }
    if (rank_cross < -5)
      rank_cross = -5;
    if (rank_cross > 1)
      rank_cross = 1;
    
    if (mDebugLevel)
      cout << "rankings: " << rank_avg_dip << " " << rank_bemc << " " << rank_cross << endl;
    primV->setRanking(rank_cross+rank_bemc+rank_avg_dip);
    if (primV->ranking() > mBestRank) {
      mBestRank = primV->ranking();
      mBestVtx = primV;
    }
  }
}

int
StMinuitVertexFinder::fit(StEvent* event)
{
    double arglist[4];

    setFlagBase();

    // get CTB info
    StCtbTriggerDetector* ctbDet = 0;
    vector<ctbHit> ctbHits;

    StTriggerDetectorCollection* trigCol = event->triggerDetectorCollection();
    if(trigCol){
      ctbDet = &(trigCol->ctb());

      float ctbSum = 0;
	    	    
      for (UInt_t slat = 0; slat < ctbDet->numberOfSlats(); slat++) {
	for (UInt_t tray = 0; tray < ctbDet->numberOfTrays(); tray++) {
	  ctbHit curHit;
	  curHit.adc = ctbDet->mips(tray,slat,0);
	  if(curHit.adc > 0){
	    ctbSum += curHit.adc;
	    ctb_get_slat_from_data(slat, tray, curHit.phi, curHit.eta);
	    ctbHits.push_back(curHit);
	  }
	}
      }
    }

    //
    //  Loop all global tracks (TPC) and store the
    //  refering helices and their estimated DCA
    //  resolution in vectors.
    //  Quality cuts are applied (see accept()).
    //  The helices and the sigma are used in
    //  fcn to calculate the fit potential which
    //  gets minimized by Minuit.
    //
    mDCAs.clear();
    mHelices.clear();
    mHelixFlags.clear();
    mSigma.clear();
    mZImpact.clear();

    fillBemcHits(event);

    double sigma;
    bool ctb_match;

    int n_ctb_match_tot = 0;
    int n_bemc_match_tot = 0;
    int n_cross_tot = 0;

    StSPtrVecTrackNode& nodes = event->trackNodes();
    mUseDCA = kFALSE;
    for (unsigned int k=0; k<nodes.size(); k++) {
      StGlobalTrack* g = ( StGlobalTrack*) nodes[k]->track(global);
      if (g && g->dcaGeometry()) {mUseDCA = kTRUE; break;}
    }
    LOG_QA << "QAInfo: StMinuitVertexFinder::fit use DCA track parameters" <<  endm;
    for (unsigned int k=0; k<nodes.size(); k++) {
      StGlobalTrack* g = ( StGlobalTrack*) nodes[k]->track(global);
      if (accept(g)&&
	  ((mUseITTF&&g->fittingMethod()==kITKalmanFitId)||
	   (!mUseITTF&&g->fittingMethod()!=kITKalmanFitId))) {	
	///LSB This should not be necessary and could be removed in future
#ifndef __alpha__
	if (!finite(g->geometry()->helix().curvature()) ){
	  LOG_WARN << "NON-FINITE curvature in track !!" << endm;
	  continue;
	}
#endif
	if (! mUseDCA ) {
	  mWidthScale = 1;
	  const StPhysicalHelixD helix = g->geometry()->helix(); 
	  Double_t r_c = sqrt(helix.xcenter()*helix.xcenter()+helix.ycenter()*helix.ycenter());
	  if (fabs(r_c - 1./helix.curvature()) > mRImpactMax)
	    continue;
	  
	  mHelices.push_back(g->geometry()->helix());
	  mHelixFlags.push_back(1);
	  
	  Double_t path=(TMath::ATan2(-helix.ycenter(),-helix.xcenter())-helix.phase())/2/TMath::Pi();
	  path *= helix.h()*helix.period();
	  StThreeVectorD tmp_pos = helix.at(path);
	  
	  Double_t z_lin=tmp_pos.z();
	  if (r_c - 1./helix.curvature() > 0) 
	    z_lin -= tmp_pos.perp()*tan(helix.dipAngle());
	  else
	    z_lin += tmp_pos.perp()*tan(helix.dipAngle());
	  
	  mZImpact.push_back(z_lin);

	  // sigma = 0.45+0.0093*::sqrt(g->length())/abs(g->geometry()->momentum()); HIJING + TRS
	  sigma = 0.6+0.0086*::sqrt(g->length())/abs(g->geometry()->momentum());
	  mSigma.push_back(sigma);
	} else { // use DCA parameters
	  mWidthScale = 0.1;// 1./TMath::Sqrt(5.);
	  StDcaGeometry* gDCA = g->dcaGeometry();
	  if (! gDCA) continue;
	  if (gDCA->impact() >  mRImpactMax) continue;
	  mDCAs.push_back(gDCA);
// 	  StPhysicalHelixD helix = gDCA->helix(); 
// 	  mHelices.push_back(helix);
	  mHelices.push_back(g->geometry()->helix());
	  mHelixFlags.push_back(1);
	  Double_t z_lin = gDCA->z();
	  mZImpact.push_back(z_lin);
	  mSigma.push_back(-1);
	  
	}
	bool shouldHitCTB = false;
	double etaInCTBFrame = -999;
	ctb_match =  EtaAndPhiToOrriginAtCTB(g,&ctbHits,shouldHitCTB,etaInCTBFrame);
	if (ctb_match) {
	  mHelixFlags[mHelixFlags.size()-1] |= kFlagCTBMatch;
	  n_ctb_match_tot++;
	}

	if (matchTrack2BEMC(g)) {
	  mHelixFlags[mHelixFlags.size()-1] |= kFlagBEMCMatch;
	  n_bemc_match_tot++;
	}

	if (checkCrossMembrane(g)) {
	  mHelixFlags[mHelixFlags.size()-1] |= kFlagCrossMembrane;
	  n_cross_tot++;
	}
      }
    }
    if (mDebugLevel)
      cout << "Found " << n_ctb_match_tot << " ctb matches, " << n_bemc_match_tot << " bemc matches, " << n_cross_tot << " tracks crossing central membrane" << endl; 
    //
    //  In case there are no tracks left we better quit
    //
    if (mHelices.empty()) {
	LOG_WARN << "StMinuitVertexFinder::fit: no tracks to fit." << endm;
	mStatusMin = -1;
	return 0;
    }
    LOG_INFO << "StMinuitVertexFinder::fit size of helix vector: " << mHelices.size() << endm;

    // Set some global pars
    if (mRequireCTB) requireCTB = true;
    
    //
    //  Reset and clear Minuit parameters
    // mStatusMin
    mMinuit->mnexcm("CLEar", 0, 0, mStatusMin);
    
    //
    //  Set parameters and start values. We do
    //  constrain the parameters since it harms
    //  the fit quality (see Minuit documentation).
    //
    // Initialize the seed with a z value which is not one of the discrete 
    // values which it can tend to, implies zero not allowed.
    // Also need different initialization when vertex constraint.

    static double step[3] = {0.03, 0.03, 0.03};

    //
    //  Scan z to find best seed for the actual fit.
    //  Skip this step if an external seed is given.
    //
    if (!mExternalSeedPresent) {
      if (!mVertexConstrain){ 
	arglist[0] = 3;
      }
      else {
        arglist[0]=1;
      }
      findSeeds();
    }
    else {
      mNSeed = 1;
    }

    Float_t old_vtx_z = -999;
    Double_t seed_z = -999;
    Double_t chisquare = 0;
    for (Int_t iSeed = 0; iSeed < mNSeed; iSeed++) {
      //
      //  Reset and clear Minuit parameters
      //  mStatusMin
      mMinuit->mnexcm("CLEar", 0, 0, mStatusMin);

      seed_z= mSeedZ[iSeed]; 

      if (mExternalSeedPresent)
	seed_z = mExternalSeed.z();
      if (mDebugLevel)
	LOG_INFO << "Vertex seed = " << seed_z << endm;
      
      if (!mVertexConstrain){ 
	mMinuit->mnparm(0, "x", 0, step[0], 0, 0, mStatusMin);
	mMinuit->mnparm(1, "y", 0, step[1], 0, 0, mStatusMin);
	mMinuit->mnparm(2, "z", seed_z, step[2], 0, 0, mStatusMin);
      }
      else {
	mMinuit->mnparm(0, "z", seed_z, step[2], 0, 0, mStatusMin);
      }

      if (!mVertexConstrain){ 
	arglist[0] = 3;
      }
      else {
	arglist[0] = 1;
      }

      int done = 0;
      int iter = 0;

      int n_trk_vtx = 0;
      Int_t n_helix = mHelices.size();
      do {  
	// For most vertices one pass is fine, but multiple passes 
	// can be done
	n_trk_vtx = 0;
	for (int i=0; i < n_helix; i++) {
	  if (fabs(mZImpact[i]-seed_z) < mDcaZMax) {
	    mHelixFlags[i] |= kFlagDcaz;
	    n_trk_vtx++;
	  }
	  else
	    mHelixFlags[i] &= ~kFlagDcaz;
	}
      	
	if (mDebugLevel) 
	  cout << n_trk_vtx << " tracks within dcaZ cut (iter " << iter <<" )" << endl;
	if (n_trk_vtx < mMinTrack) {
	  if (mDebugLevel) 
	    cout << "Less than mMinTrack (=" << mMinTrack << ") tracks, skipping vtx" << endl;
	  continue;
	}
	mMinuit->mnexcm("MINImize", 0, 0, mStatusMin);
	done = 1;

	//
	//  Check fit result
	//

	if (mStatusMin) {
	  LOG_WARN << "StMinuitVertexFinder::fit: error in Minuit::mnexcm(), check status flag. ( iter=" << iter << ")" << endm;
	  done = 0; // refit
	}

	Double_t fedm, errdef;
	Int_t npari, nparx;

	mMinuit->mnstat(chisquare, fedm, errdef, npari, nparx, mStatusMin);

	if (mStatusMin != 3) {
	  cout << "Warning: Minuit Status: " << mStatusMin << ", func val " << chisquare<< endl;
	  done = 0;  // refit
	}
	mMinuit->mnhess();

	double new_z, zerr;
	if (!mVertexConstrain) {
	  mMinuit->GetParameter(2, new_z, zerr); 
	}
	else {
	  mMinuit->GetParameter(0, new_z, zerr); 
	}

	if (fabs(new_z - seed_z) > 1) // refit if vertex shifted
	  done = 0;

	int n_trk = 0;
	for (int i=0; i < n_helix; i++) {
	  if ( fabs(mZImpact[i] - new_z) < mDcaZMax ) {
	    n_trk++;
	  }
	}
	if ( 10 * abs(n_trk - n_trk_vtx) >= n_trk_vtx ) // refit if number of track changed by more than 10%
	  done = 0;

	iter++;
	seed_z = new_z; // seed for next iteration
      } while (!done && iter < 5 && n_trk_vtx >= mMinTrack);

      if (n_trk_vtx < mMinTrack)
	continue;

      if (!done) { 
	LOG_WARN << "Vertex unstable: no convergence after " << iter << " iterations. Skipping vertex " << endm;
	continue;
      }

      if (!mExternalSeedPresent && fabs(seed_z-mSeedZ[iSeed]) > mDcaZMax)
	LOG_WARN << "Vertex walks during fits: seed was " << mSeedZ[iSeed] << ", vertex at " << seed_z << endl;

      if (fabs(seed_z - old_vtx_z) < 0.1) {
	if (mDebugLevel) 
	  cout << "Vertices too close (<0.1). Skipping" << endl;
	continue;
      }

      // Store vertex
      StThreeVectorD XVertex;
      Float_t cov[6];
      memset(cov,0,sizeof(cov));
    
      Double_t val, verr;
#if 1
      if (!mVertexConstrain) {
	XVertex = StThreeVectorD(mMinuit->fU[0],mMinuit->fU[1],mMinuit->fU[2]);
	Double_t emat[9];
	/* 0 1 2
	   3 4 5
           6 7 8 */
	mMinuit->mnemat(emat,3);
	cov[0] = emat[0];
	cov[1] = emat[3];
	cov[2] = emat[4];
	cov[3] = emat[6];
	cov[4] = emat[7];
	cov[5] = emat[8];
      }
      else {
	mMinuit->GetParameter(0, val, verr); 
	XVertex.setZ(val);  cov[5]=verr*verr;
	
	// LSB Really error in x and y should come from error on constraint
	// At least this way it is clear that those were fixed paramters
	XVertex.setX(beamX(val));  cov[0]=0.1; // non-zero error values needed for Sti
	XVertex.setY(beamY(val));  cov[2]=0.1;
      }
#else
      if (!mVertexConstrain) {
	mMinuit->GetParameter(0, val, verr); 
	XVertex.setX(val);  cov[0]=verr*verr;
	
	mMinuit->GetParameter(1, val, verr); 
	XVertex.setY(val);  cov[2]=verr*verr;

	mMinuit->GetParameter(2, val, verr); 
	XVertex.setZ(val);  cov[5]=verr*verr;
      
      }
      else {
	mMinuit->GetParameter(0, val, verr); 
	XVertex.setZ(val);  cov[5]=verr*verr;
      
	// LSB Really error in x and y should come from error on constraint
	// At least this way it is clear that those were fixed paramters
	XVertex.setX(beamX(val));  cov[0]=0.1; // non-zero error values needed for Sti
	XVertex.setY(beamY(val));  cov[2]=0.1;
      }
#endif

      StPrimaryVertex primV;
      primV.setPosition(XVertex);
      primV.setChiSquared(chisquare);  // this is not really a chisquare, but anyways
      primV.setCovariantMatrix(cov); 
      primV.setVertexFinderId(minuitVertexFinder);
      primV.setFlag(1); // was not set earlier by this vertex finder ?? Jan
      primV.setRanking(333);
      primV.setNumTracksUsedInFinder(n_trk_vtx);

      int n_ctb_match = 0;
      int n_bemc_match = 0;
      int n_cross = 0;
      n_trk_vtx = 0;

      double mean_dip = 0;
      double sum_pt = 0;
      for (unsigned int i=0; i < mHelixFlags.size(); i++) {
	if (!(mHelixFlags[i] & kFlagDcaz))
	  continue;
	n_trk_vtx++;
	if (mHelixFlags[i] & kFlagCTBMatch)
	  n_ctb_match++;
	if (mHelixFlags[i] & kFlagBEMCMatch)
	  n_bemc_match++;
	if (mHelixFlags[i] & kFlagCrossMembrane)
	  n_cross++;
	mean_dip += mHelices[i].dipAngle();
	sum_pt += mHelices[i].momentum(0).perp();
      }
      
      mean_dip /= n_trk_vtx;

      if (mDebugLevel) {
	cout << "check n_trk_vtx " << n_trk_vtx << ", found " << n_ctb_match << " ctb matches, " << n_bemc_match << " bemc matches, " << n_cross << " tracks crossing central membrane" << endl; 
	cout << "mean dip " << mean_dip << endl;
      }
      primV.setNumMatchesWithCTB(n_ctb_match);      
      primV.setNumMatchesWithBEMC(n_bemc_match);
      primV.setNumTracksCrossingCentralMembrane(n_cross);
      primV.setMeanDip(mean_dip);
      primV.setSumOfTrackPt(sum_pt);

      //..... add vertex to the list
      addVertex(&primV);

      old_vtx_z = XVertex.z();
    }

    calculateRanks();

    //  Reset the flag which tells us about external
    //  seeds. This needs to be provided for every event.
    mExternalSeedPresent = false;
    
    requireCTB = false;

    return 1;
} 
//________________________________________________________________________________
Double_t StMinuitVertexFinder::Chi2atVertex(StThreeVectorD &vtx) {
  Double_t f = 0;
  Double_t e, s;
  nCTBHits = 0;
  if (! mUseDCA) {
    for (unsigned int i=0; i<mHelices.size(); i++) {
      if ((mHelixFlags[i] & kFlagDcaz) && (!requireCTB||(mHelixFlags[i] & kFlagCTBMatch))) {
	s = mWidthScale*mSigma[i];
	e = mHelices[i].distance(vtx, false);  // false: don't do multiple loops
	f -= exp(-(e*e)/(2*s*s));  // robust potential
	if((mHelixFlags[i] & kFlagCTBMatch) && e<3.0) nCTBHits++;
      }
    }
  } else {// UseDCA
    for (unsigned int i=0; i<mDCAs.size(); i++) {
      if ((mHelixFlags[i] & kFlagDcaz) && (!requireCTB||(mHelixFlags[i] & kFlagCTBMatch))) {
	const StDcaGeometry* gDCA = mDCAs[i];
	if (! gDCA) continue;
	const StPhysicalHelixD helix = gDCA->helix();
	e = helix.distance(vtx, false);  // false: don't do multiple loops
#if 0
	StPhysicalHelixD helixA = mHelices[i];
        Double_t eA = helixA.distance(vtx, false);
	// linear approximation, ignore error in angles
	StThreeVectorD dcaP = helix.at(e) - vtx;
	StThreeVectorD dcaPA = helixA.at(eA) - vtx;
	Double_t impact = gDCA->impact();
	if (impact);
	Double_t x = dcaP.x();
	Double_t y = dcaP.y();
	StThreeVectorF mom = gDCA->momentum();
	Double_t px = mom.x();
	Double_t py = mom.y();
	Double_t Imp = TMath::Sqrt(x*x + y*y);
	Double_t sgn = x*py-y*px;
	Double_t sgnA = helixA.curvatureSignedDistance(dcaPA);
        if (sgnA);
	if (sgn < 0) Imp = - Imp;
	Double_t Z  = dcaP.z();
#endif
	const Float_t *errMatrix = gDCA->errMatrix();
#if 0
	Double_t d = errMatrix[0]*errMatrix[2] - errMatrix[1]*errMatrix[1];
	if (TMath::Abs(d) < 1e-7) continue;
	Double_t a[3] = {errMatrix[2]/d, -errMatrix[1]/d, errMatrix[0]/d};
// 	Double_t chi2     = a[0]*Imp*Imp + 2*a[1]*Imp*Z + a[2]*Z*Z; 
#endif
	Double_t chi2     = e*e/(errMatrix[0] + errMatrix[2]);
	Double_t scale = 1./(mWidthScale*mWidthScale);
	f += scale*(1. - TMath::Exp(-chi2/scale)); // robust potential
	//	f -= scale*TMath::Exp(-chi2/scale); // robust potential
	if((mHelixFlags[i] & kFlagCTBMatch) && e<3.0) nCTBHits++;
      }
    }
  }
  return f;
}
//________________________________________________________________________________

void StMinuitVertexFinder::fcn1D(int& npar, double* gin, double& f, double* par, int iflag)
{
    double z = par[0];
    double x = beamX(z);
    double y = beamY(z);
    StThreeVectorD vtx(x,y,z);
    f = Chi2atVertex(vtx);
}
void StMinuitVertexFinder::fcn(int& npar, double* gin, double& f, double* par, int iflag)
{
  StThreeVectorD vtx(par);
  f = Chi2atVertex(vtx);
}


bool
StMinuitVertexFinder::accept(StTrack* track) const
{
    //
    //   Accept only tracks which fulfill certain
    //   quality criteria.
    //

    return (track &&
	    track->flag() >= 0 &&
	    track->fitTraits().numberOfFitPoints() >= mMinNumberOfFitPointsOnTrack &&
	    !track->topologyMap().trackFtpc() &&
	    finite(track->length()) &&  //LSB another temporary check
	    track->geometry()->helix().valid());
}


/// Use mMinuit print level 
void
StMinuitVertexFinder::setPrintLevel(int level) 
{
  mMinuit->SetPrintLevel(level);
}

void
StMinuitVertexFinder::printInfo(ostream& os) const
{

    os << "StMinuitVertexFinder - Statistics:" << endl;
    os << "Number of vertices found ......." << size() << endl;
    os << "Rank of best vertex ............" << mBestRank << endl;
    if(mBestVtx) {
      os << "Properties of best vertex:" << endl;
      os << "position ..................... " << mBestVtx->position() << endl;
      os << "position errors .............. " << mBestVtx->positionError()<< endl;
      os << "# of used tracks ............. " << mBestVtx->numTracksUsedInFinder() << endl;
      os << "Chisquare .................... " << mBestVtx->chiSquared() << endl;
    }
    os << "min # of fit points for tracks . " << mMinNumberOfFitPointsOnTrack << endl;
    os << "final potential width scale .... " << mWidthScale << endl;
}

void StMinuitVertexFinder::UseVertexConstraint(double x0, double y0, double dxdz, double dydz, double weight) {
  mVertexConstrain = true;
  mX0 = x0;
  mY0 = y0;
  mdxdz = dxdz;
  mdydz = dydz;
  mWeight = weight;
  LOG_INFO << "StMinuitVertexFinder::Using Constrained Vertex" << endm;
  LOG_INFO << "x origin = " << mX0 << endm;
  LOG_INFO << "y origin = " << mY0 << endm;
  LOG_INFO << "slope dxdz = " << mdxdz << endm;
  LOG_INFO << "slope dydz = " << mdydz << endm;
  LOG_INFO << "weight in fit = " << weight <<  endm;
  StThreeVectorD origin(mX0,mY0,0.0);
  double pt  = 88889999;   
  double nxy=::sqrt(mdxdz*mdxdz +  mdydz*mdydz);
    if(nxy<1.e-5){ // beam line _MUST_ be tilted
      LOG_WARN << "StMinuitVertexFinder:: Beam line must be tilted!" << endm;
      nxy=mdxdz=1.e-5; 
    }
    double p0=pt/nxy;  
    double px   = p0*mdxdz;
    double py   = p0*mdydz;
    double pz   = p0; // approximation: nx,ny<<0
    StThreeVectorD MomFstPt(px*GeV, py*GeV, pz*GeV);
    delete mBeamHelix;
    mBeamHelix = new StPhysicalHelixD(MomFstPt,origin,0.5*tesla,1.);

    //re-initilize minuit for 1D fitting
    mMinuit = new TMinuit(1);         
    mMinuit->SetFCN(&StMinuitVertexFinder::fcn1D);
    mMinuit->SetPrintLevel(1);
    mMinuit->SetMaxIterations(1000);
    mExternalSeedPresent = false;


}


double StMinuitVertexFinder::beamX(double z) {
  float x = mX0 + mdxdz*z;
  return x;
}

double StMinuitVertexFinder::beamY(double z) {
  float y = mY0 + mdydz*z;
  return y;
}

int  StMinuitVertexFinder::NCtbMatches() { 
  return nCTBHits;
}
int  StMinuitVertexFinder::NCtbSlats() { 
  return -777; // dum result, perhaps not needed at all,JB
}


//void StMinuitVertexFinder::SetFitPointsCut(int fitpoints) {mMinNumberOfFitPointsOnTrack = fitpoints;return;}


