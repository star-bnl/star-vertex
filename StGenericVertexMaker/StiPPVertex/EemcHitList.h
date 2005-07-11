#ifndef EemcHitList_h
#define EemcHitList_h

#include "ScintHitList.h"
#include "StEEmcUtil/EEfeeRaw/EEdims.h"
class StEmcDetector;
class StEEmcDbMaker;
class EEmcGeomSimple;

class EemcHitList : public ScintHitList {
 private:

  StEEmcDbMaker* eeDb; 
  EEmcGeomSimple *geomE;
  int name2bin[MaxSectors][MaxSubSec][MaxEtaBins]; // map --> my bin
  const float *etaHL; // limits of eta bins

  //params
  uint killStatEEmc;
 
 public:
 EemcHitList(StEEmcDbMaker* x, uint y, EEmcGeomSimple *z);
  virtual  ~EemcHitList();
  void clear();
  void initRun();
  void build( StEmcDetector*det, float adcMin);
  virtual  int etaBin(float eta);
  virtual float bin2EtaLeft(int iEta);

};

#endif