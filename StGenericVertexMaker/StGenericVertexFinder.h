/*!
 * \class StGenericVertexFinder
 *
 * \author Lee Barnby, April 2003
 *
 * (pseudo) Base class for vertex finders
 *
 *
 * $Id: StGenericVertexFinder.h,v 1.18 2009/07/09 00:16:12 genevb Exp $
 */

#ifndef STAR_StGenericVertexFinder
#define STAR_StGenericVertexFinder

#include "StEventTypes.h"
#include "StPrimaryVertex.h"

class StEvent;

class StGenericVertexFinder {
 public:
  // virtual and '=0' ; those MUST be implemented
  virtual ~StGenericVertexFinder();                           // virtual destructor
  virtual int            fit(StEvent*)=0;                     // fit the vertex

  StPrimaryVertex*       getVertex(int idx) const;
  void                   addVertex(StPrimaryVertex*);
  int                    size() const;
  virtual void           UseVertexConstraint(double, double, double, double, double)=0;
          void           NoVertexConstraint();
          int            IsVertexConstraint() const {return mVertexConstrain;}
  virtual void           UsePCT(bool usePCT = true);
  virtual void           CalibBeamLine(){ /* noop */;} // overload if useful

  virtual void           printInfo(ostream& = cout) const=0;

  // General (default)
  virtual void           SetMode(Int_t mode=0 ) {mMode = mode;}
  virtual int            GetMode() const 	{return mMode;}
          void           SetDebugLevel(Int_t level) {mDebugLevel=level;}
  virtual void           Init(){ /* noop */;}
  virtual void           Finish(){ /* noop */;}
  virtual void           InitRun  (int runumber){ /* noop */;}
  virtual void           Clear();
  const std::vector<StPrimaryVertex> *result() {return &mVertexList;} 
 
  void                   FillStEvent(StEvent*) const;

 protected: //................................

  StGenericVertexFinder();
 private:
  vector<StPrimaryVertex> mVertexList;      // Holds all found prim veritcess

 protected: //................................
  StPrimaryVertexOrder   mVertexOrderMethod; // will default to 0 i.e. orderByNumberOfDaughters
  bool                   mVertexConstrain;   // Use vertex constraint from db
  int                    mMode;              // used for any Finder behavior change
  int                    mDebugLevel;

};



// $Log: StGenericVertexFinder.h,v $
// Revision 1.18  2009/07/09 00:16:12  genevb
// Create a calib mode for StGenericVertex when using VtxSeedCalG
//
// Revision 1.17  2008/10/23 20:37:31  genevb
// Add switches for turning on/off use of Post-Crossing Tracks [default:off]
//
// Revision 1.16  2006/04/26 15:37:04  jeromel
// mVertexOrderMethod (To be tested)
//
// Revision 1.15  2006/04/08 00:18:10  mvl
// Added member for debuglevel
//
// Revision 1.14  2005/07/19 21:45:53  perev
// MultiVertex
//
// Revision 1.13  2005/06/21 02:16:36  balewski
// multiple prim vertices are stored in StEvent
//
// Revision 1.12  2005/03/11 22:23:53  balewski
// towards PPV
//
// Revision 1.11  2005/03/09 19:24:18  balewski
// preparation for PPV vertex finder
//
// Revision 1.10  2004/12/13 20:39:58  fisyak
// Add initaition of StGenericVertexFinder variables, replace mDumMaker by StMaker::GetChain() method
//
// Revision 1.9  2004/09/13 15:41:30  balewski
// fix bug in ppLMV4/5 switch
//
// Revision 1.8  2004/09/03 00:09:08  jeromel
// Modified code to Implement Init() and SetMode() and allow passing a switch
// to chose the vertex finder from within the same code implementation. Was
// needed for ppLMV (one implementation, two algorithm)
//
// Revision 1.7  2004/08/04 21:57:56  balewski
// toward smarter ppLMV5
//
// Revision 1.6  2004/07/24 19:40:38  balewski
// fix swap of vert & errVert
//
// Revision 1.5  2004/07/24 02:57:40  balewski
// clean up of ppLMV, CTB-util separated
//
// Revision 1.4  2004/07/23 02:24:38  jeromel
// Oops ... Worng swithc (had twice Minuit). Now corrected.
//
// Revision 1.3  2004/07/23 00:58:19  jeromel
// Base class method+data member (was duplicated in implementation)
//
// Revision 1.2  2004/04/06 02:43:43  lbarnby
// Fixed identification of bad seeds (no z~0 problem now). Better flagging. Message manager used.
//
// Revision 1.1  2003/05/09 22:22:36  lbarnby
// Initial revision: a base class for STAR (StEvent-based) vertex finders
//
#endif
