#ifndef BemcHitList_h
#define BemcHitList_h

#include "ScintHitList.h"
class StEmcDetector;
class StBemcTables;
class StEmcGeom ;

class BemcHitList : public ScintHitList {
 private:
  enum {mxm=120,mxe=20,mxs=2};
  int mes2bin[mxm][mxe][mxs]; // map {m,e,s}--> my bin
  StBemcTables *myTable;
  StEmcGeom *geomB;
  // params:
  float kSigPed;

 public:
 BemcHitList();
  virtual  ~BemcHitList();
  void clear();
  void initRun(bool isMC);
  void build( StEmcDetector*det, float adcMin);
  virtual   int etaBin(float eta);
  virtual float bin2EtaLeft(int iEta);
};

#endif
