/***************************************************************************
 * $Id$
 *
 * Author: Lee Barnby, April 2003
 *
 ***************************************************************************
 * Description: Base class for vertex finders
 *
 ***************************************************************************/

#ifndef STAR_StGenericVertexFinder
#define STAR_StGenericVertexFinder

#include "StEventTypes.h"


class StEvent;

class StGenericVertexFinder {
 public:
  virtual bool fit(StEvent*)=0;
  virtual int status() const =0;
  virtual StThreeVectorD result() const=0;
  virtual StThreeVectorD error() const=0;
  //virtual ~StGenericVertexFinder();
  void FillStEvent(StEvent*) const;
};

// $Log$
#endif
