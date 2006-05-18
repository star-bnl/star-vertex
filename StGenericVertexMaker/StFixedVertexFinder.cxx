/*
 *  StFixedVertexFinder.cxx
 *  $Id$
 *
 *  Author Lee Barnby (University of Birmingham) May 2006.
 *
 */

#include <StEventTypes.h>
#include <StEnumerations.h>
#include <StMessMgr.h>


#include "StFixedVertexFinder.h"
#include "StMcEvent.hh"
#include "StMcVertex.hh"
#include "StMaker.h"


StFixedVertexFinder::StFixedVertexFinder(){
    mFixedX=mFixedY=mFixedZ=0.0;
    mMode=0;
}

int StFixedVertexFinder::fit(StEvent* event)
{
  StPrimaryVertex  prim;
  StVertexFinderId VFId;
  Float_t cov[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; // All errors are zero

  if (mMode == 0){
    // Really the default which takes the SetPos() TBI
    LOG_DEBUG << "StFixedVertexFinder::fit() fixing a vertex" << endm;
    StThreeVectorD pos(mFixedX,mFixedY,mFixedZ);
    prim.setPosition(pos);
    VFId = undefinedVertexFinder;

  } else {
    
    LOG_DEBUG << "StFixedVertexFinder::fit() reading a vertex from StMcEvent" << endm;
    StMcEvent* mcEvent= (StMcEvent*) StMaker::GetChain()->GetDataSet("StMcEvent");
    
    if ( ! mcEvent){
      // Not sure what to do here, will leave and do nothing (so no vertex)
      // and it will be super obvious that something went terribly wrong
      LOG_ERROR << "StFixedVertexFinder::fit() no StMcEvent" << endm;
      return 0;

    } else {
      StMcVertex *pv = mcEvent->primaryVertex();
      StThreeVectorD pos(pv->position().x(),
    			 pv->position().y(),
    			 pv->position().z());
      //cout << "DEBUG :: saving McVertex " << pos << endl;
      prim.setPosition(pos);
    }
    VFId = mcEventVertexFFinder;
  }

  prim.setCovariantMatrix(cov);
  prim.setFlag(1);                    // So that we know its the primary vertex
  prim.setRanking(-5);                // Have to have something
  prim.setVertexFinderId(VFId);       // Id depends on MC or fixed position used
  addVertex(&prim);

  return size();
}

void StFixedVertexFinder::printInfo(ostream& os)const{
    os << "StFixedVertexFinder - fixed vertex" << endl;
    os << "Fixed position: x=" << mFixedX << " y=" << mFixedY << " z=" << mFixedZ << endl;
}

void StFixedVertexFinder::UseVertexConstraint(double x0, double y0, double dxdz, double dydz, double weight){
    LOG_WARN << "StFixedVertexFinder::UseVertexConstraint() - vertex beam constraint NOT implemented in context of fixed vertex finder" << endm;

}

void StFixedVertexFinder::SetVertexPosition(double x, double y, double z){
    mFixedX=x;
    mFixedY=y;
    mFixedZ=z;
}

/*
 * $Log$
 */
