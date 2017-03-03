#ifndef EemcHitList_h
#define EemcHitList_h
#include <sys/types.h>
#include "StGenericVertexMaker/StiPPVertex/ScintHitList.h"
#include "StEEmcUtil/EEfeeRaw/EEdims.h"
#include "StEEmcUtil/database/cstructs/eemcConstDB.hh"
class StEmcDetector;
class StEEmcDb;
class EEmcGeomSimple;

class EemcHitList : public ScintHitList {
 private:

  StEEmcDb* eeDb; 
  EEmcGeomSimple *geomE;
  int name2bin[MaxSectors][MaxSubSec][MaxEtaBins]; // map --> my bin
  const Float_t *etaHL; // limits of eta bins

  //params
  unsigned int killStatEEmc;
 
 public:
 EemcHitList(StEEmcDb* x=nullptr, unsigned int y=EEMCSTAT_ONLPED|EEMCSTAT_STKBT|EEMCSTAT_HOTHT|EEMCSTAT_HOTJP|EEMCSTAT_JUMPED, EEmcGeomSimple *z=nullptr);
  virtual  ~EemcHitList();
  void clear();
  void initRun();
  void build( StEmcDetector*det, float adcMin);
  virtual  int etaBin(float eta);
  virtual float bin2EtaLeft(int iEta);

};

#endif
