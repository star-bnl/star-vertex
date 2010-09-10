#ifndef VertexData_h
#define VertexData_h
/*********************************************************************
 * $Id$
 *********************************************************************
 * full description of found vertex
 */

#include <TVector3.h>

class VertexData {
 public:
  int id; // vertex ID assigned by PPV
  TVector3 r,er; // vertex position and its error
  int nUsedTrack; // # of tracks used to identify the vertex
  float Lmax; // maximum of the likelhood function.
  float gPtSum; // total tranverse momentum of used tracks.
  int nBtof,nCtb,nBemc,nEemc,nTpc,nAnyMatch; // number of matched tracks - dongx
  int nBtofV,nCtbV,nBemcV,nEemcV,nTpcV,nAnyVeto; // number of vetoed tracks - dongx
  
  // methods
  VertexData();
  void print(ostream& os) const;
};
#endif


/*
 * $Log$
 * Revision 1.2  2005/08/30 22:08:43  balewski
 * drop '*' from declaration of   mTrackData &  mVertexData
 *
 * Revision 1.1  2005/07/11 20:38:13  balewski
 * PPV added for real
 *

 *
 *
 *********************************************************************/
