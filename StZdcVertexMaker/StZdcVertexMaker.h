/***************************************************************************
 *
 * $Id: StZdcVertexMaker.h,v 1.1 2001/08/30 16:17:36 macross Exp $
 *
 * Author:  Johan E. Gonzalez, August 2001
 ***************************************************************************
 *
 * Description: This is the .h file for StZdcVertexMaker.cxx
 *
 ***************************************************************************
 *
 * $Log: StZdcVertexMaker.h,v $
 * Revision 1.1  2001/08/30 16:17:36  macross
 * Initial Revition.
 *
 **************************************************************************/
#ifndef StZdcVertexMaker_hh
#define StZdcVertexMaker_hh

#include "StMaker.h"

class St_ZdcCalPars;

class StZdcVertexMaker : public StMaker {
public:
    StMaker* currentChain;
    StZdcVertexMaker(const char* name = "StZdcVertexMaker",
                     const char* title = "event/StZdcVertexMaker");
    ~StZdcVertexMaker();
    void  Clear(const char* opt="");
    Int_t Init();
    Int_t Make();
    Int_t Finish();
    
private:    
    float mEAP0;
    float mEAP1;
    float mEAP2;
    float mEAP3;
    float mWAP0;
    float mWAP1;
    float mWAP2;
    float mWAP3;
    float mVPAR;
    float mOFF;
    
    // the following is a ROOT macro  that is needed in all ROOT accessible code
    ClassDef(StZdcVertexMaker, 1)
};
#endif



