/*!
 * \class  StKinkMaker
 * \brief  Class to find kink secondary vertices
 * \author Camelia Mironov, KSU
 * \date   Jan,2004
 *
 * This is the C++ Kink maker replacement for the table-based
 * older code.
 *
 */


#ifndef STAR_StKinkMaker
#define STAR_StKinkMaker

#include "StMaker.h"
#include "StThreeVectorD.hh"
#include "TMath.h"
#include "StTrackGeometry.h"
#include "StV0FinderMaker.h"

class StKinkLocalTrack;
class St_tkf_tkfpar;
class tkf_tkfpar_st;
class StPhysicalHelixD;
class StEvent;
class StKinkVertex;
class StTrack;
class pairD;

//enum TrackerUsage{ See StV0FinderMaker.h
//  kTrackerUseTPT  = 0,
//  kTrackerUseITTF = 1,
//  kTrackerUseBOTH = 2
//};

class StKinkMaker : public StMaker {
public: 
  //  StKinkMaker(const char* name);
  StKinkMaker(const char* name="KinkMaker", TrackerUsage use=kTrackerUseTPT);
  virtual  ~StKinkMaker(); 
  virtual  Int_t  Init();
  virtual  Int_t  Make();
  virtual void SetTrackerUsage(Int_t opt=0){UseTracker=opt;}
  //  virtual Int_t GetTrackerUsage(){return UseTracker;}

  virtual const char *GetCVS() const
  {static const char cvs[]="Tag $Name:  $ $Id: StKinkMaker.h,v 1.1 2004/03/03 01:44:04 cmironov Exp $ built "__DATE__" "__TIME__ ; return cvs;}


 private:
  // Bool_t            m_kinkEvalOn;   //switch for the evaluation
  
  StKinkLocalTrack *mTrack1;   //!       
  StKinkLocalTrack *mTrack2;   //!    
  StTrack *mParentTrackCandidate;              //! 
  StTrack *mDaughterTrackCandidate;              //!
         
  StThreeVectorD   mEventVertex;       //position of primary vertex
  StThreeVectorD   mParentMoment, mDaughterMoment;
  StThreeVectorD   mKinkVertex;
  Float_t          mParentImpact, mDaughterImpact;
  Float_t          mDca, mDecayAngle;

  double mBfield;//mg field 
  protected:
  
  void FillEvent(StTrackGeometry *myDaughterGeometry1, 
                 StTrackGeometry *myParentGeometry11,
                 StTrack *mDaughterTrackCandidate,
                 const StThreeVectorD &mParentMomentum,
                 const StThreeVectorD &mDaughterMomentum);
  
  bool acceptTrack(StTrack *tr, int mITTFflag, Int_t mUseTracker);
 

  St_tkf_tkfpar* m_tkfpar;
  StEvent* event;
  StKinkVertex* kinkVertex;

  int UseTracker;
  int ITTFflag;
  ClassDef(StKinkMaker, 1)  
};
    
#endif
    


   
  