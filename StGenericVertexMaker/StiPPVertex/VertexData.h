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
  int nCtb,nBemc,nEemc,nTpc,nAnyMatch; // number of matched tracks 
  int  nCtbV,nBemcV,nEemcV,nTpcV,nAnyVeto; // number of vetoed tracks
  
  // methods
  VertexData();
  void print(ostream& os);
};
#endif


/*
 * $Log$

 *
 *
 *********************************************************************/
