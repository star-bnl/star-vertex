/***************************************************************************
 * $Id: StGenericVertexFinder.cxx,v 1.6 2004/12/13 20:39:58 fisyak Exp $
 *
 * Author: Lee Barnby, April 2003
 *
 ***************************************************************************
 * Description: Base class for vertex finders
 *
 ***************************************************************************/
#include "StGenericVertexFinder.h"
#include "StMessMgr.h"
#include "StMaker.h"


StGenericVertexFinder::StGenericVertexFinder() : 
  mUseITTF(false), mFlagBase(0), mBeamHelix(0), mVertexConstrain(false), mWeight(0),
  mRequireCTB(false),  mExternalSeedPresent(false), mStatus(0), mMode(0), 
  mMinNumberOfFitPointsOnTrack(0) {}

/*!
  Adds the vertex to StEvent (currently as a primary)
  Here we invent our own flag and other data to put in
  In real life we have to get it from somewhere (as done for position)
*/
void StGenericVertexFinder::FillStEvent(StEvent* event) const{

  Float_t ex,ey,ez; // Position errors 
  ex = this->error().x();ey = this->error().y();ez = this->error().z();
  Float_t cov[6] = {ex*ex,0.0,ey*ey,0.0,0.0,ez*ez};

  Float_t xSq = 5.43;
  Float_t probXSq = 0.2468;

  StPrimaryVertex* primV = new StPrimaryVertex();
  primV->setPosition(this->result());             //requires StThreeVectorF
//   primV->setFlag(mFlagBase+this->status());       //requires unsigned int
  primV->setFlag(1);                              // Should conform to dst_vertex.idl flag definition.
  primV->setCovariantMatrix(cov);                 //requires float[6]
  primV->setChiSquared(xSq);                      //requires float
  primV->setProbChiSquared(probXSq);              //requires float

  //primV->setParent();  //requires StTrack* but we won't use this, also
  //addDaughter(StTrack*) and removeDaughter(StTrack*) not used here
  //addDaughter would be used when filling primary tracks in later maker

  event->addPrimaryVertex(primV);
  gMessMgr->Debug() << "StGenericVertexFinder::FillStEvent: Added new primary vertex" << endm;
}



void
StGenericVertexFinder::setExternalSeed(const StThreeVectorD& s)
{
    mExternalSeedPresent = true;
    mExternalSeed = s;
}


void StGenericVertexFinder::NoVertexConstraint() 
{
  mVertexConstrain = false; 
  gMessMgr->Info() << "StGenericVertexFinder::No Vertex Constraint" << endm;
}

void StGenericVertexFinder::setFlagBase()
{
  if(mUseITTF){
    mFlagBase = 8000;
  } else {
    mFlagBase = 1000;
  }
}


// $Log: StGenericVertexFinder.cxx,v $
// Revision 1.6  2004/12/13 20:39:58  fisyak
// Add initaition of StGenericVertexFinder variables, replace mDumMaker by StMaker::GetChain() method
//
// Revision 1.5  2004/07/30 22:59:00  calderon
// Setting the primary vertex flag to 1 for the moment, as per
// dst_vertex.idl.  This was causing the FTPC code to reject the
// primary vertex used as their seed.
//
// Revision 1.4  2004/07/24 02:57:40  balewski
// clean up of ppLMV, CTB-util separated
//
// Revision 1.3  2004/07/23 00:57:43  jeromel
// Base class method implementation
//
// Revision 1.2  2004/04/06 02:43:43  lbarnby
// Fixed identification of bad seeds (no z~0 problem now). Better flagging. Message manager used.
//
// Revision 1.1  2003/05/09 22:22:46  lbarnby
// Initial revision: a base class for STAR (StEvent-based) vertex finders
//
