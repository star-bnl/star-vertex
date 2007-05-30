/// *-- Author : David Hardtke, based on Jan Balewski
// Revision 1.1.1.1  2001/01/31 14:00:07  balewski
// First release
//
//
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//   Maker to run minuit based vertex finder                            //
//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <strings.h>
#include <math.h>

#include "StGenericVertexMaker.h"
#include "St_DataSetIter.h"
#include "StEventTypes.h"
#include "TH2.h"
#include "TNtuple.h"
#include "StMessMgr.h"

#include "StGenericVertexFinder.h"
#include "StppLMVVertexFinder.h"
#include "StFixedVertexFinder.h"

#include "StTreeMaker/StTreeMaker.h"

// Vertex finder implemtations
#include "Minuit/StMinuitVertexFinder.h"
#include "StiPPVertex/StPPVertexFinder.h"


#include "tables/St_g2t_vertex_Table.h" // tmp for Dz(vertex)
#include "tables/St_vertexSeed_Table.h" //

// for Helix model
#include "StarCallf77.h"
extern "C" {void type_of_call F77_NAME(gufld,GUFLD)(float *x, float *b);}
#define gufld F77_NAME(gufld,GUFLD)
#include "SystemOfUnits.h"
#ifndef ST_NO_NAMESPACES
using namespace units;
#endif


ClassImp(StGenericVertexMaker)
//___________________________________________________________
StGenericVertexMaker::StGenericVertexMaker(const char *name):StMaker(name)
{
  usebeamline = kFALSE;
  useCTB = kFALSE;
  eval = kFALSE;
  nEvTotal=nEvGood=0;
  externalFindUse=kTRUE; ///Default means that no finding actually done
  use_ITTF=kTRUE;
  m_Mode2=0;
}
//_____________________________________________________________________________
StGenericVertexMaker::~StGenericVertexMaker()
{

  if(theFinder) delete theFinder;

}

/*!
  The Init() method instantiates the VertexFinder() method.
  Since this is  a Maker, the switch between the diverse methods
  will be made as part of the m_Mode mechanism. m_Mode will be
  a bit set to later allow multiple vertex finder running in the
  same pass. The structure does not allow this to happen for now
  (will need to have stack-like of VertexFinders and loop over
  them adding vertices in the collection)

  m_Mode = 0x1     Minuit
  m_Mode = 0x2     ppLMV4  This will not be able to run in parrallele of ppLMV5
  m_Mode = 0x4     ppLMV5  This will not be able to run in parrallele of ppLMV4
  m_Mode = 0x8     PPV with CTB matching
  m_Mode = 0x10    PPV without CTB matching
  m_Mode = 0x20    Fixed vertex finder
  m_Mode = 0x40    Fixed vertex finder, read from MC event
 
  Default          Minuit  (to preserver backward compatibility)

  All VertexFinder-s need to have the same methods (like DoUseITTF()
  NCtbMatches() etc ...) described in the GenericVertexFinder() class).
  Currentely, methods are not part of the base class and need
  cleanup.

*/
Int_t StGenericVertexMaker::Init()
{
  // setup params

  LOG_INFO << "StGenericVertexMaker::Init: m_Mode=" <<  m_Mode <<" m_Mode2=" <<  m_Mode2 <<  endm;
  bool isMinuit=false;

  if ( m_Mode & 0x1){
    theFinder= new StMinuitVertexFinder();
    isMinuit=true;

  } else if ( m_Mode & 0x2){
    theFinder= new StppLMVVertexFinder();
    theFinder->SetMode(0);                 // this mode is an internal to ppLMV option switch

  } else if ( m_Mode & 0x4){
    theFinder= new StppLMVVertexFinder();
    theFinder->SetMode(1);                 // this mode is an internal to ppLMV option switch

  } else if ( m_Mode & 0x8 ||  m_Mode & 0x10){ // 2 version of PPV w/ & w/o CTB
    LOG_INFO << "StGenericVertexMaker::Init: uses PPVertex finder"<<  endm;
    theFinder= new StPPVertexFinder();
    if ( m_Mode & 0x10) ((StPPVertexFinder*) theFinder)->useCTB(false);	
    if(GetMaker("emcY2")) {//very dirty, but detects if it is M-C or real data
      ((StPPVertexFinder*) theFinder)->setMC(true);
    }

  } else if ( m_Mode & 0x20 || m_Mode & 0x40) {
      theFinder = new StFixedVertexFinder();
      if (m_Mode & 0x40){
	LOG_INFO << "StGenericVertexMaker::Init: fixed vertex using MC vertex" << endm;
	theFinder->SetMode(1);
      } else {
	LOG_INFO << "StGenericVertexMaker::Init: fixed vertex 'finder' selected" << endm;
      }

  } else {
    // Later, this would NEVER make multiple possible vertex
    // finder unlike for option 0x1 .
    theFinder= new StMinuitVertexFinder();
    isMinuit=true;

  }
  
  if(isMinuit) { // this is ugly, one should abort at 'else' above, Jan
    if (use_ITTF)  ((StMinuitVertexFinder*)theFinder)->DoUseITTF();
    if (useCTB) ((StMinuitVertexFinder*)theFinder)->CTBforSeed();
  } else {
    assert(!eval); // current implementation support only Minuit Vertex finder, JB 
  }
 
  //    theFinder->UseVertexConstraint(-0.265,0.4088,-0.00135,0.0004333,0.0001);
  //theFinder->UseVertexConstraint(0.0,0.0,0.0,0.0,0.0001);
  if (eval) mEvalNtuple = new TNtuple("results","results","thX:thY:thZ:thStat:goodGlob:evX:evY:evZ:evStat:nPrim:nCTB:geantX:geantY:geantZ");

  theFinder->Init();
  return StMaker::Init();
}

//_____________________________________________________________________________
void StGenericVertexMaker::Clear(const char* opt){
  LOG_INFO <<" StGenericVertexMaker::Clear()"<<endm; 
  theFinder->Clear();
}


//_____________________________________________________________________________
Int_t StGenericVertexMaker::InitRun(int runnumber){
  theFinder->InitRun(runnumber);
  if (usebeamline) {
     double x0 = 0.;
     double y0 = 0.;
     double dxdz = 0.;
     double dydz = 0.;

     // Get Current Beam Line Constraint from database
     TDataSet* dbDataSet = this->GetDataBase("Calibrations/rhic");

     if (dbDataSet) {
       vertexSeed_st* vSeed = ((St_vertexSeed*) (dbDataSet->FindObject("vertexSeed")))->GetTable();

     x0 = vSeed->x0;
     y0 = vSeed->y0;
     dxdz = vSeed->dxdz;
     dydz = vSeed->dydz;
     }
     else {
       LOG_INFO << "StGenericVertexMaker -- No Database for beamline" << endm;
     }
     LOG_INFO << "BeamLine Constraint: " << endm;
     LOG_INFO << "x(z) = " << x0 << " + " << dxdz << " * z" << endm;
     LOG_INFO << "y(z) = " << y0 << " + " << dydz << " * z" << endm;
     theFinder->UseVertexConstraint(x0,y0,dxdz,dydz,0.0001);
  }
  return StMaker::InitRun(runnumber);
}


//_____________________________________________________________________________
Int_t StGenericVertexMaker::Finish()
{

  LOG_INFO << "StGenericVertexMaker::Finish " <<GetName() <<endm;
  LOG_INFO << " Total events: " << nEvTotal << endm;
  LOG_INFO << " Good events:  " << nEvGood  << endm;


  //LSB TODO Leave this for now. Should really be using STAR/ROOT I/O scheme?
  if (eval) {
   TFile out("MinuitVertexEval.root","RECREATE");
   mEvalNtuple->Write();
   out.Close();
  }
  
  if(theFinder) theFinder->Finish();  
  return  kStOK;
}

//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________
Bool_t StGenericVertexMaker::DoFit(){
  StThreeVectorD myvertex;

  StEvent *event = (StEvent *) GetInputDS("StEvent");
  assert(event);

  if (theFinder->fit(event)) {
    theFinder->printInfo();
  }  else {
    LOG_ERROR << "StGenericVertexMaker::DoFit: vertex fit failed, no vertex." << endm;
    return kFALSE;
  }

  return kTRUE;

}
//_____________________________________________________________________________
//_____________________________________________________________________________
//_____________________________________________________________________________


//_____________________________________________________________________________
Int_t StGenericVertexMaker::Make()
{
  nEvTotal++;
  primV  = NULL;
  mEvent = NULL;
  mEvent = (StEvent *)GetInputDS("StEvent");
  LOG_DEBUG << "StGenericVertexMaker::Make: StEvent pointer " << mEvent << endm;
  LOG_DEBUG << "StGenericVertexMaker::Make: external find use " << externalFindUse << endm;

  if(!externalFindUse){
    DoFit();
  }

//   //For testing purposes
//   theFinder->DoNotUseITTF();
//   if (theFinder->fit(mEvent)) theFinder->printInfo();
//   if (theFinder->status()!=-1) theFinder->FillStEvent(mEvent);

//   theFinder->DoUseITTF();
//   if (theFinder->fit(mEvent)) theFinder->printInfo();
//   if (theFinder->status()!=-1) {
//     theFinder->FillStEvent(mEvent);
//     nEvGood++;
//   }


  if (eval)MakeEvalNtuple();

  if(!externalFindUse){
    ///Only fill StEvent when successful
    if (theFinder->size()>0){
      theFinder->FillStEvent(mEvent);
      nEvGood++;
    }
  }
  return kStOK;
}

//-----------------------------------------------------------------------------

void StGenericVertexMaker::MakeEvalNtuple(){ // only for Minuit vertex finder

  // get geant vertex
  St_DataSet *gds=GetDataSet("geant");
  St_g2t_vertex  *g2t_ver=0;
  g2t_vertex_st *gver=0;
  if(gds)  g2t_ver=( St_g2t_vertex *)gds->Find("g2t_vertex");
  if(g2t_ver)gver=g2t_ver->GetTable();

  double gx = -999.;
  double gy = -999.;
  double gz = -999.;

  if(gver) {
    gx=gver->ge_x[0];
    gy=gver->ge_x[1];
    gz=gver->ge_x[2];
  }


  //  G E T     P R I M     V E R T E X
  primV=mEvent->primaryVertex();
  int nCtb= ((StMinuitVertexFinder*)theFinder)->NCtbMatches(); 
  int stat= ((StMinuitVertexFinder*)theFinder)->statusMin();
  
  if(!primV) {
    LOG_INFO <<"primaryVertex()=NULL"<<endm;
    // why would one access x,y,z of the vertex if it is not found, Jan ???
    float x=999,y=999,z=999;
    mEvalNtuple->Fill(x,y,z,stat,mEvent->summary()->numberOfGoodTracks(),-999.,-999.,-999.,-999.,-999.,nCtb,gx,gy,gz);
    }
  else {
     LOG_INFO << Form("primaryVertex()= %f, %f %f, nTracks=%d\n",primV->position().x(),primV->position().y(),primV->position().z(),primV->numberOfDaughters())<<endm;
     mEvalNtuple->Fill(primV->position().x(),primV->position().y(),primV->position().z(),stat               ,mEvent->summary()->numberOfGoodTracks(),primV->position().x(),primV->position().y(),primV->position().z(),primV->flag(),primV->numberOfDaughters(), nCtb,gx,gy,gz);
  }
}

//____________________________________________________________________________
