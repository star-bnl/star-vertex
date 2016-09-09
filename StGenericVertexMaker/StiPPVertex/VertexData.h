#ifndef VertexData_h
#define VertexData_h
/*********************************************************************
 * $Id: VertexData.h,v 1.3 2010/09/10 21:08:35 rjreed Exp $
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
  int nBtof,nCtb,nBemc,nEemc,nTpc,nAnyMatch; // number of matched tracks
  int nBtofV,nCtbV,nBemcV,nEemcV,nTpcV,nAnyVeto; // number of vetoed tracks
  
  // methods
  VertexData();
  VertexData(const TVector3& position);
  void print(ostream& os) const;
};
#endif


/*
 * $Log: VertexData.h,v $
 * Revision 1.3  2010/09/10 21:08:35  rjreed
 * Added function UseBOTF and bool mUseBtof to switch the use of the TOF on and off in vertex finding.  Default value is off (false).
 * Added functions, and variables necessary to use the TOF in PPV for vertex finding.  Includes matching tracks to the TOF and changing the track weight based on its matched status with the TOF.
 *
 * Revision 1.2  2005/08/30 22:08:43  balewski
 * drop '*' from declaration of   mTrackData &  mVertexData
 *
 * Revision 1.1  2005/07/11 20:38:13  balewski
 * PPV added for real
 *

 *
 *
 *********************************************************************/
