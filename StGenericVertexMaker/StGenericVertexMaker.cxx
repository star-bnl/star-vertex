//*-- Author : David Hardtke, based on Jan Balewski
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
#include "StChain.h"
#include "St_DataSetIter.h"
#include "StEventTypes.h"
#include "TH2.h"
#include "TNtuple.h"
#include "StMessMgr.h"

#include "StGenericVertexFinder.h"  
#include "StMinuitVertexFinder.h"  

#include "StTreeMaker/StTreeMaker.h"

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
//_____________________________________________________________________________
StGenericVertexMaker::StGenericVertexMaker(const char *name):StMaker(name)
{
  usebeamline = kFALSE;
  useCTB = kFALSE;
  eval = kFALSE;
  nEvTotal=nEvGood=0;
  externalFindUse=kTRUE; ///Default means that no finding actually done
}
//_____________________________________________________________________________
StGenericVertexMaker::~StGenericVertexMaker()
{
}

//_____________________________________________________________________________
Int_t StGenericVertexMaker::Init()
{
  // setup params
  EtaCut=1.4; // Sensible default cut
  theFinder=new StMinuitVertexFinder();
  if (use_ITTF) theFinder->DoUseITTF();
  //    theFinder->CTBforSeed();
  //    theFinder->UseVertexConstraint(-0.265,0.4088,-0.00135,0.0004333,0.0001);
  //theFinder->UseVertexConstraint(0.0,0.0,0.0,0.0,0.0001);
  if (eval) mEvalNtuple = new TNtuple("results","results","thX:thY:thZ:thStat:goodGlob:evX:evY:evZ:evStat:nPrim:nCTB:geantX:geantY:geantZ");
  return StMaker::Init();
}

//_____________________________________________________________________________
Int_t StGenericVertexMaker::InitRun(int runnumber){
  if (useCTB) theFinder->CTBforSeed();
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
       cout << "StGenericVertexMaker -- No Database for beamline" << endl;
     }   
     cout << "BeamLine Constraint: " << endl;
     cout << "x(z) = " << x0 << " + " << dxdz << " * z" << endl;
     cout << "y(z) = " << y0 << " + " << dydz << " * z" << endl << endl;
     theFinder->UseVertexConstraint(x0,y0,dxdz,dydz,0.0001);
  }
  return StMaker::InitRun(runnumber);
}


//_____________________________________________________________________________
Int_t StGenericVertexMaker::Finish()
{
  //LSB TODO change over to using message manager
  cout <<" Finish::"<<GetName() <<endl;
  printf(" Total events: %d \n", nEvTotal);
  printf(" Good events: %d \n", nEvGood);

  //LSB TODO Leave this for now. Should really be using STAR/ROOT I/O scheme?
  if (eval) {
   TFile out("MinuitVertexEval.root","RECREATE");
   mEvalNtuple->Write();
   out.Close();
  }

  //LSB TODO check whether this is correct usage
  if(theFinder) delete theFinder;

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
    myvertex = theFinder->result();
    theFinder->printInfo();
  }  else {
    cout << "Error: vertex fit failed, no vertex." << endl;
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
  primV=NULL;
  mEvent=NULL;
  mEvent = (StEvent *)GetInputDS("StEvent"); 

  if(!externalFindUse){
    DoFit();
  } 

  if (eval)MakeEvalNtuple();

  if(!externalFindUse){
    ///Only fill StEvent when successful
    if (theFinder->status()!=-1){
      theFinder->FillStEvent(mEvent); 
      nEvGood++;
    }
  }
  return kStOK;
}

//-----------------------------------------------------------------------------

void StGenericVertexMaker::MakeEvalNtuple(){

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
  if(!primV) {
    printf("primaryVertex()=NULL\n");
    mEvalNtuple->Fill(theFinder->result().x(),theFinder->result().y(),theFinder->result().z(),theFinder->status(),mEvent->summary()->numberOfGoodTracks(),-999.,-999.,-999.,-999.,-999.,theFinder->NCtbMatches(),gx,gy,gz);
    }
  else {
    printf("primaryVertex()= %f, %f %f, nTracks=%d\n",primV->position().x(),primV->position().y(),primV->position().z(),primV->numberOfDaughters());  
  mEvalNtuple->Fill(theFinder->result().x(),theFinder->result().y(),theFinder->result().z(),theFinder->status(),mEvent->summary()->numberOfGoodTracks(),primV->position().x(),primV->position().y(),primV->position().z(),primV->flag(),primV->numberOfDaughters(),theFinder->NCtbMatches(),gx,gy,gz);
  }
}

//____________________________________________________________________________
// LSB Commented out since moved to finder
// void const StGenericVertexMaker::FillStEvent(){
//   //Adds the vertex to StEvent (currently as a primary)
//   // Here we invent our own flag and other data to put in
//   // In real life we have to get it from somewhere (as done for position)
//   UInt_t minuitFlag=1000;
//   Float_t cov[6] = {0.1,0.2,0.3,0.4,0.5,0.6};
//   Float_t xSq = 5.43;
//   Float_t probXSq = 0.2468;

//   StPrimaryVertex* primV = new StPrimaryVertex();
//   primV->setPosition(theFinder->result());    //requires StThreeVectorF
//   primV->setFlag(minuitFlag+theFinder->status());       //requires unsigned int
//   primV->setCovariantMatrix(cov);      //requires float[6]
//   primV->setChiSquared(xSq);           //requires float
//   primV->setProbChiSquared(probXSq);       //requires float
//   //primV->setParent();  //requires StTrack* but we won't use this, also
//   //addDaughter(StTrack*) and removeDaughter(StTrack*) not used here
//   //addDaughter would be used when filling primary tracks in later maker

//   mEvent->addPrimaryVertex(primV);
//   gMessMgr->Debug()
//     << "StGenericVertexMaker::FillStEvent: Added new primary vertex" << endm;

// }

//------------  N O T   U S E D   -------------------------------


