/***************************************************************************
 *
 * $Id$
 *
 * Author:  Johan E. Gonzalez, August 2001
 ***************************************************************************
 *
 * Description: This is the .h file for StZdcVertexMaker.cxx
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.3  2003/09/10 19:47:42  perev
 * ansi corrs
 *
 * Revision 1.2  2001/08/31 19:07:36  macross
 * Modified code to retrieve ADC and TDC pulses from TrgDet table
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
    virtual const char *GetCVS() const
    {static const char cvs[]="Tag $Name$ $Id$ built "__DATE__" "__TIME__ ; return cvs;}
    
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
    ClassDef(StZdcVertexMaker,0)
};
#endif



