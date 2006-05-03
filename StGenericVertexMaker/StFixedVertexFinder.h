/*!
 *  \class StFixedVertexFinder
 *  \author Lee Barnby (University of Birmingham) May 2006.
 *  \brief StGenericVertexFinder implementation for fixing vertex.
 *
 *  $Id$
 *
 */

#ifndef STAR_StFixedVertexFinder
#define STAR_StFixedVertexFinder


#include "StGenericVertexFinder.h"

class StEvent;

class StFixedVertexFinder: public StGenericVertexFinder{
public:
    StFixedVertexFinder();
    // should we have destructor too?
    
    // mandatory implementations
    int fit(StEvent*);
    void printInfo(ostream& = cout)const;
    void UseVertexConstraint(double x0, double y0, double dxdz, double dydz, double weight);
    
private:
    Double_t mFixedX;
    Double_t mFixedY;
    Double_t mFixedZ;

};

/***************************************************************************
*
* $Log$
*
**************************************************************************/
#endif
