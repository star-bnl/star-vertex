/*!
 * \class StGenericVertexFinder
 *
 * \author Lee Barnby, April 2003
 *
 * (pseudo) Base class for vertex finders
 *
 *
 * We have the option of requiring that at least one track matches
 * the CTB during the first two scans for the vertex (these are coarse scans
 * to locate the probable z vertex
 *
 * myvertex.CTBforSeed();
 * myvertex.NoCTBforSeed();
 *
 * During the final fit (once the z position of the vertex has been constrained)
 * there is no CTB requirement.  To get the number of tracks with CTB match:
 *
 *  myvertex.NCtbMatches()
 *
 * Because NCtbMatches() may be handled multiple ways, the implementation
 * is enforced.
 *
 *
 * $Id$
 */

#ifndef STAR_StGenericVertexFinder
#define STAR_StGenericVertexFinder

#include "StEventTypes.h"


class StEvent;

class StGenericVertexFinder {
 public:
  // virtual and '=0' ; those MUST be implemented
  virtual ~StGenericVertexFinder(){};                         // virtual destructor
  virtual bool           fit(StEvent*)=0;                     // fit the vertex
  virtual int            NCtbMatches()=0;                     // returns the number of CTB match

  // General (default)
  virtual StThreeVectorD result() const {return mFitError;};  // result of fit
  virtual StThreeVectorD error()  const {return mFitResult;}; // error on fit result
  virtual int            status() const {return mStatus;};    // status flag

  void                   FillStEvent(StEvent*) const;
  void                   CTBforSeed(){   mRequireCTB = true;}
  void                   NoCTBforSeed(){ mRequireCTB = false;}

  void                   DoUseITTF(){    mUseITTF=kTRUE; };
  void                   DoNotUseITTF(){ mUseITTF=kFALSE;};
  void                   setFlagBase();

  void                   setExternalSeed(const StThreeVectorD&);
  void                   NoVertexConstraint();
  void                   SetFitPointsCut(int fitpoints) {mMinNumberOfFitPointsOnTrack = fitpoints;};



 protected:
  bool                   mUseITTF;          // Use only tracks with ITTF encoded method
  UInt_t                 mFlagBase;         // ITTF track flag

  StPhysicalHelixD*      mBeamHelix;        // Beam Line helix
  bool                   mVertexConstrain;  // Use vertex constraint from db

  double                 mWeight ;          // Weight in fit for vertex contraint
  bool                   mRequireCTB;       // Set maker to use CTB

  StThreeVectorD         mExternalSeed;
  bool                   mExternalSeedPresent;


  StThreeVectorD         mFitResult;        // fit result
  StThreeVectorD         mFitError;         // fit errors
  int                    mStatus;           // status flag 

  unsigned int           mMinNumberOfFitPointsOnTrack;


};



// $Log$
// Revision 1.2  2004/04/06 02:43:43  lbarnby
// Fixed identification of bad seeds (no z~0 problem now). Better flagging. Message manager used.
//
// Revision 1.1  2003/05/09 22:22:36  lbarnby
// Initial revision: a base class for STAR (StEvent-based) vertex finders
//
#endif
