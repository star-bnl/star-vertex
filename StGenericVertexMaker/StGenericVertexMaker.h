/*!
 * \class StGenericVertexMaker
 * \author David Hardtke, based on Jan Balewskis template
 *
 * Maker for minuit based vertex finder
 * Lee Barnby - modification, becomes StGenericVertexMaker
 *
 * $Id: StGenericVertexMaker.h,v 1.16 2014/08/06 11:43:19 jeromel Exp $
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
  Bool_t  useITTF;
  Bool_t  useBeamline;
  Bool_t  calibBeamline;
  Bool_t  useCTB;
  Bool_t  usePCT;
  Bool_t  useBTOF;
  Bool_t  eval;
  Bool_t  externalFindUse; /// Finder will by called externally (by StiMaker)
  Int_t   minTracks;

  TNtuple *mEvalNtuple;    /// Ntuple for evaluation purposes

  StEvent *mEvent;
  StPrimaryVertex* primV;
  StGenericVertexFinder *theFinder;

  Bool_t DoFit(); ///Find and fit the primary vertex
  void const FillStEvent();
  void MakeEvalNtuple();

  Int_t nEvTotal,nEvGood;

 public: 
  StGenericVertexMaker(const char *name="GenericVertex");
  virtual       ~StGenericVertexMaker();
  virtual Int_t Init();
  virtual Int_t InitRun (Int_t runumber);
  virtual void  Clear(const char* opt="");
  virtual Int_t Finish();
  virtual Int_t Make();
  StGenericVertexFinder* GetGenericFinder(){return (StGenericVertexFinder*)theFinder;};

  void UseBeamLine()            {SetAttr("BeamLine"       , kTRUE );}
  void DoNotUseBeamLine()       {SetAttr("BeamLine"       , kFALSE);}
  void CalibBeamLine()          {SetAttr("calibBeamline"  , kTRUE );}
  void UseCTB()                 {SetAttr("CTB"            , kTRUE );}
  void DoNotUseCTB()            {SetAttr("CTB"            , kFALSE);}
  void DoEval()                 {SetAttr("eval"           , kTRUE );}
  void SetInternalFind()        {SetAttr("externalFindUse", kFALSE);}
  void SetUseITTF()             {SetAttr("ITTF"           , kTRUE );}
  void SetDoNotUseITTF()        {SetAttr("ITTF"           , kFALSE);}
  void SetMinimumTracks(Int_t n){SetAttr("minTracks"      , n     );}
  void UsePCT()                 {SetAttr("PCT"            , kTRUE );}
  void DoNotUsePCT()            {SetAttr("PCT"            , kFALSE);}

  virtual const char *GetCVS() const
    {static const char cvs[]="Tag $Name:  $ $Id: StGenericVertexMaker.h,v 1.16 2014/08/06 11:43:19 jeromel Exp $ built " __DATE__ " " __TIME__ ; return cvs;}
  
  ClassDef(StGenericVertexMaker, 0)   //StAF chain virtual base class for Makers
};
    
#endif
    



