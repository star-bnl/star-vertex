/*!
 * \class StGenericVertexMaker
 * \author  David Hardtke, based on Jan Balewskis template
 *
 * Maker for minuit based vertex finder
 * Lee Barnby - modification, becomes StGenericVertexMaker
 *
 * $Id: StGenericVertexMaker.h,v 1.3 2004/07/23 00:56:58 jeromel Exp $
 *
 */


#ifndef STAR_StGenericVertexMaker
#define STAR_StGenericVertexMaker

#ifndef StMaker_H
#include "StMaker.h"
#endif


class StEvent;
class StPrimaryVertex;
class StMinuitVertexFinder;
class StGenericVertexFinder;
class TNtuple;

class StGenericVertexMaker : public StMaker 
{
 private: 
  // control and cuts
  Bool_t  use_ITTF;
  Bool_t  usebeamline;
  Bool_t  useCTB;
  Bool_t  eval;
  Bool_t  externalFindUse; /// Finder will by called externally (by StiMaker)
  Float_t EtaCut;          /// Tracks with larger eta not considered

  TNtuple *mEvalNtuple;    /// Ntuple for evaluation purposes

  StEvent *mEvent;
  StPrimaryVertex* primV;
  StMinuitVertexFinder *theFinder;

  Bool_t DoFit(); ///Find and fit the primary vertex
  void const FillStEvent();
  void MakeEvalNtuple();

  Int_t nEvTotal,nEvGood;

 protected:

 public: 
  StGenericVertexMaker(const char *name="GenericVertex");
  virtual       ~StGenericVertexMaker();
  virtual Int_t Init();
  virtual Int_t InitRun  (int runumber);
  virtual Int_t Finish();
  virtual Int_t  Make();

  inline StGenericVertexFinder* GetGenericFinder(){return (StGenericVertexFinder*)theFinder;};

  inline void UseBeamLine(){usebeamline = kTRUE;};
  inline void DoNotUseBeamLine(){usebeamline = kFALSE;};
  inline void UseCTB(){useCTB = kTRUE;};
  inline void DoNotUseCTB(){useCTB = kFALSE;};
  inline void DoEval(){eval= kTRUE;};
  inline void SetInternalFind(){externalFindUse=kFALSE;};
  inline void SetUseITTF(){use_ITTF=kTRUE;};
  inline void SetDoNotUseITTF(){use_ITTF=kFALSE;};

  virtual const char *GetCVS() const
    {static const char cvs[]="Tag $Name:  $ $Id: StGenericVertexMaker.h,v 1.3 2004/07/23 00:56:58 jeromel Exp $ built "__DATE__" "__TIME__ ; return cvs;}
  
  ClassDef(StGenericVertexMaker, 0)   //StAF chain virtual base class for Makers
};
    
#endif
    



