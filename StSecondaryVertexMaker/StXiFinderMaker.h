/**
  \class StXiFinderMaker
  
  StXiFinderMaker finds Xi secondary vertices

*/

#ifndef StXiFinderMaker_hh
#define StXiFinderMaker_hh

#include "StV0FinderMaker.h"

class St_exi_exipar;
class exi_exipar_st;
class StXiVertex;


class StXiFinderMaker : public StV0FinderMaker {

 public:
  StXiFinderMaker(const char* name="XiFinderMaker");
  virtual ~StXiFinderMaker();

  virtual Int_t Init();
  virtual Int_t Make();

  virtual Bool_t UseV0();

  virtual const char *GetCVS() const
  {static const char cvs[]="Tag $Name$ $Id$ built "__DATE__" "__TIME__ ; return cvs;}

 protected:
  St_exi_exipar* exipar;           //!
  exi_exipar_st* parsXi;           //!
  StXiVertex* xiVertex;            //!
  int det_id_xi;

  ClassDef(StXiFinderMaker,0)

};

#endif

//_____________________________________________________________________________
// $Id$
// $Log$
//
