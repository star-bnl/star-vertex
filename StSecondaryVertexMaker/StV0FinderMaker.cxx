//////////////////////////////////////////////////////////////////////////
///                                                                    ///
/// StV0FinderMaker class (finds V0 secondary vertices)                ///
///                                                                    ///
//////////////////////////////////////////////////////////////////////////
//
//  Cuts can be found in the code by comments beginning with "Cut:"
//
//


#include "StV0FinderMaker.h"
#include "StMessMgr.h"
#include "StEvent/StEventTypes.h"
#include "TMath.h"
#include "TVector2.h"
#include "tables/St_ev0_ev0par2_Table.h"

///Begin Betty
///#include "StEstMaker/StEstTracker.h"
#include "StEvent/StTrack.h"
///End Betty

#include "math_constants.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"

//Stuff to get correct Bfield
#include "StarCallf77.h"   
 extern "C" {void type_of_call F77_NAME(gufld,GUFLD)(float *x, float *b);}   
 #define gufld F77_NAME(gufld,GUFLD) 


#define MAXTRACKS 10000
#define MAXSTRACK  6500




static StTrack* STATtrk[MAXTRACKS];
static unsigned short STATptrk[MAXSTRACK];
static unsigned short STATntrk[MAXSTRACK];
static short STAThits[MAXTRACKS];
static short STATdetId[MAXTRACKS];
static double STATpt[MAXTRACKS];
static double STATptot[MAXTRACKS];
static StPhysicalHelixD STATheli[MAXTRACKS];
static unsigned short STATtrkID[MAXTRACKS];



ClassImp(StV0FinderMaker)
 
//_____________________________________________________________________________
  StV0FinderMaker::StV0FinderMaker(const char *name):StMaker(name),
         ev0par2(0),pars(0),pars2(0),event(0),v0Vertex(0),
         prepared(kFALSE),useExistingV0s(kFALSE),dontZapV0s(kFALSE),
         useTracker(kTrackerUseTPT),useSVT(kNoSVT),useV0Language(kV0LanguageUseCpp),
         useXiLanguage(kXiLanguageUseCppOnCppV0),useLanguage(kLanguageUseRun),
         useLikesign(kLikesignUseStandard),useRotating(kRotatingUseStandard)
{
  // Assign pointers for static arrays
  maxtracks = MAXTRACKS;
  trk = STATtrk;
  ptrk = STATptrk;
  ntrk = STATntrk;
  hits = STAThits;
  detId = STATdetId;
  pt = STATpt;
  ptot = STATptot;
  heli = STATheli;
  trkID = STATtrkID;

  // Check for multiple instances
  if (hits[maxtracks-1] == -2)
    gMessMgr->Warning() << "StV0FinderMaker(" << name <<
      ") : MORE THAN ONE INSTANCE!" << endm;
  else hits[maxtracks-1] = -2;

  ptV0sq = 3.5*3.5;
  Bfield = 1.e-10; //Random value for initialisation.
}










//_____________________________________________________________________________
StV0FinderMaker::~StV0FinderMaker() {
}












//_____________________________________________________________________________
void StV0FinderMaker::GetPars()
{
  TDataSet* dbDataSet = GetDataBase("global/vertices");
  if (!dbDataSet) {
    gMessMgr->Error(
      "StV0FinderMaker::Init() : could not find appropriate database.");
    return; 
  }
  ev0par2 = (St_ev0_ev0par2*) (dbDataSet->FindObject("ev0par2"));
  if (!ev0par2) {
    gMessMgr->Error(
      "StV0FinderMaker::Init() : could not find ev0par2 in database.");
    return;
  }
  ///AddRunCont(ev0par2);
}













//_____________________________________________________________________________
Int_t StV0FinderMaker::Init()
{bool a,b,c;
 
 if ((useTracker!=kTrackerUseTPT) && (useTracker!=kTrackerUseITTF) && (useTracker!=kTrackerUseBOTH))
    {gMessMgr->Error("StV0FinderMaker::Init() : wrong TrackerUsage parameter set.");
     return kStErr;
     }
 if ((useSVT!=kNoSVT) && (useSVT!=kUseSVT))
    {gMessMgr->Error("StV0FinderMaker::Init() : wrong SVTUsage parameter set.");
     return kStErr;
     }
 if ((useLikesign!=kLikesignUseStandard) && (useLikesign!=kLikesignUseLikesign))
    {gMessMgr->Error("StV0FinderMaker::Init() : wrong LikesignUsage parameter set.");
     return kStErr;
     }
 if ((useRotating!=kRotatingUseStandard) && (useRotating!=kRotatingUseRotating) && (useRotating!=kRotatingUseSymmetry) && (useRotating!=kRotatingUseRotatingAndSymmetry))
    {gMessMgr->Error("StV0FinderMaker::Init() : wrong RotatingUsage parameter set.");
     return kStErr;
     }
 
 if (useTracker == kTrackerUseTPT) gMessMgr->Info()<<"StV0FinderMaker : use TPT tracks."<<endm;
 if (useTracker == kTrackerUseITTF) gMessMgr->Info()<<"StV0FinderMaker : use ITTF tracks."<<endm;
 if (useTracker == kTrackerUseBOTH) gMessMgr->Info()<<"StV0FinderMaker : use TPT *and* ITTF tracks."<<endm;
 if (useSVT == kUseSVT) gMessMgr->Info()<<"StV0FinderMaker : use SVT points if possible."<<endm;///Betty
 if (useLikesign == kLikesignUseLikesign) gMessMgr->Info()<<"StV0FinderMaker : does like-sign finding."<<endm;
 if (useRotating == kRotatingUseRotating) gMessMgr->Info()<<"StV0FinderMaker : does rotating finding."<<endm;
 if (useRotating == kRotatingUseSymmetry) gMessMgr->Info()<<"StV0FinderMaker : does symmetry finding."<<endm;
 if (useRotating == kRotatingUseRotatingAndSymmetry) gMessMgr->Info()<<"StV0FinderMaker : does rotating + symmetry finding."<<endm;

 if (useLanguage != kLanguageUseSpecial)
    {a=(bool)(1&(useLanguage>>2));
     b=(bool)(1&(useLanguage>>1));
     c=(bool)(1&useLanguage);
     useV0Language=2*(!(a^c))+(a|c);
     useXiLanguage=4*(b&(!(a^c)))+2*(a&b&(!c))+(a|c);
     }
 switch (useLanguage)
    {case kLanguageUseOldRun : gMessMgr->Info()<<"StV0FinderMaker : Fortran run."<<endm;
                               break;
     case kLanguageUseRun : gMessMgr->Info()<<"StV0FinderMaker : C++ run."<<endm;
                            gMessMgr->Info()<<"StV0FinderMaker : BE CAREFUL : you are NOT running the XiFinder !"<<endm;
                            break;
     case kLanguageUseTestV0Finder : gMessMgr->Info()<<"StV0FinderMaker : Test V0Finder."<<endm;
                                     break;
     case kLanguageUseTestXiFinder : gMessMgr->Info()<<"StV0FinderMaker : Test XiFinder."<<endm;
                                     gMessMgr->Info()<<"StV0FinderMaker : BE CAREFUL : you are NOT running the XiFinder !"<<endm;
                                     break;
     case kLanguageUseTestBothFinders : gMessMgr->Info()<<"StV0FinderMaker : Test V0Finder and XiFinder."<<endm;
                                        gMessMgr->Info()<<"StV0FinderMaker : BE CAREFUL : you are NOT running the XiFinder !"<<endm;
                                        break;
     case kLanguageUseSpecial : break;
     default : gMessMgr->Error("StV0FinderMaker::Init() : wrong LanguageUsage parameter set.");
               return kStErr;
     }
 if ((useV0Language!=kV0LanguageUseFortran) && (useV0Language!=kV0LanguageUseCpp) && (useV0Language!=kV0LanguageUseBoth))
    {gMessMgr->Error("StV0FinderMaker::Init() : wrong V0LanguageUsage parameter set.");
     return kStErr;
     }
 if (1&useV0Language) gMessMgr->Info()<<"StV0FinderMaker :    Will store Fortran V0."<<endm;
 if (2&useV0Language) gMessMgr->Info()<<"StV0FinderMaker :    Will store C++ V0."<<endm;
 if (1&useXiLanguage) gMessMgr->Info()<<"StV0FinderMaker :    Will store Fortran Xi."<<endm;
 if (2&useXiLanguage) gMessMgr->Info()<<"StV0FinderMaker :    BE CAREFUL : will NOT store C++ Xi, although asked."<<endm;
 if (4&useXiLanguage) gMessMgr->Info()<<"StV0FinderMaker :    BE CAREFUL : will NOT store C++ Xi, although asked."<<endm;

 return StMaker::Init();
 }

//____________________________________________________________________________

Int_t StV0FinderMaker::InitRun( int RunNumber){
	  float gufldX[3]= {0,0,0};
	  float gufldB[3];
	  gufld(gufldX,gufldB);
	  Bfield = gufldB[2]*kilogauss;
}












//_____________________________________________________________________________
Int_t StV0FinderMaker::Prepare() {

  if (prepared) return kStOk;

  unsigned short i,j,nNodes;
  StThreeVectorD p;

  // Get pars
  GetPars();
  ITTFflag=kITKalmanFitId;

  // Get event 
  event = (StEvent*) GetInputDS("StEvent");
  if (!event)
     {gMessMgr->Warning("StV0FinderMaker : no StEvent ; skipping event.");
      return kStWarn;
      }

  // Get Primary Vertex Position
  StPrimaryVertex* pvert = event->primaryVertex();
  if (!pvert)
     {gMessMgr->Warning("StV0FinderMaker : no primary vertex ; skipping event.");
      return kStWarn;
      }
  mainv = pvert->position();

  StSPtrVecTrackNode& theNodes = event->trackNodes();
  nNodes = theNodes.size();
  
  // Find which global tracks to use
  trks=0;
  ntrks=0;
  ptrks=0;
  for (i=0; i<nNodes; i++) {
    for (j=0; j<theNodes[i]->entries(global); j++) {

      StTrack* tri = theNodes[i]->track(global,j);
      ///Begin Betty
      if(useSVT){
        StTrack* svtTrack = theNodes[i]->track(estGlobal,j);
        if (svtTrack){  //if there is a track that uses an SVT point, set the track to estGlobal
          tri=svtTrack;
        }
      }
      ///End Betty
      //Cut: track type
      if ((tri->fittingMethod() != ITTFflag && (GetTrackerUsage() == kTrackerUseITTF)) ||
          (tri->fittingMethod() == ITTFflag && (GetTrackerUsage() == kTrackerUseTPT))) continue;

      //Cut: track flag
      if (tri->flag() <= 0) continue;

        // Determine detector id of track i
        const StTrackTopologyMap& map = tri->topologyMap();
        Bool_t tpcHit = map.hasHitInDetector(kTpcId);
        Bool_t silHit = map.hasHitInDetector(kSvtId) ||
                        map.hasHitInDetector(kSsdId);
        if (tpcHit) {
          if (silHit)
            detId[trks] = 3; //SVT+TPC
          else
            detId[trks] = 1; //TPC-only
        } else if (silHit)
          detId[trks] = 2; //SVT-only
        else
          //ignore this track
          continue;

        trk[trks] = tri;

        StTrackGeometry* triGeom = tri->geometry();
        heli[trks] = triGeom->helix();
       
        p = triGeom->momentum();

	pt[trks] = p.perp();
        ptot[trks] = p.mag();
	trkID[trks]=tri->key();

        // Determine number of hits (in SVT+TPC)
        hits[trks] = map.numberOfHits(kTpcId) +
                     map.numberOfHits(kSvtId) +
                     map.numberOfHits(kSsdId);

        if (triGeom->charge() > 0) ptrk[ptrks++] = trks;
        else if (triGeom->charge() < 0) ntrk[ntrks++] = trks;
        trks++;
    }
  }

  gMessMgr->Info() << "StV0FinderMaker : No. of nodes is : "
                   << nNodes << endm;
  gMessMgr->Info() << "StV0FinderMaker : No. of tracks is : "
                   << trks << endm;

  prepared = kTRUE;
  return kStOk;
}











//_____________________________________________________________________________
Int_t StV0FinderMaker::Make() {

  // Variables:
  StThreeVectorD xi,xj,pi,pj,xpp,pp,impact;      // 3D vectors of the V0
  StThreeVectorD xi1,xj1,pi1,pj1,tmp3V;          // temporary 3D vectors
  TVector2 ri,rj,xci,xcj,tmp2V;                  // 2D vectors of the tracks
  double rad_i,rad_j,separation,solution,dxc;    // helix circle params
  double dca_ij,dca_ij1,rmin,dlen,pperpsq,ppmag; // V0 params
  double alpha,ptArm_sq,pPosAlongV0,pNegAlongV0; // Armenteros params
  double cosij,sin2ij,t1,t2;                     // 3D dca calculation vars
  unsigned short i,j,ii,jj;                      // track iteration vars
  pairD paths,path2;                             // helix pathLength vars
  double temp;
  Bool_t doSecond, isPrimaryV0, usedV0;
  Int_t iRes;

  ///Julien
  pairD paths1;
  
  ///Betty
  long  keepTrack;

  if (! (2&useV0Language)) return kStOk;

  gMessMgr->Info("StV0FinderMaker::Make() : Starting...");
  
  // Prepare event and track variables
  iRes = Prepare();
  if (iRes != kStOk) return iRes;


  StSPtrVecV0Vertex& v0Vertices = event->v0Vertices();
  gMessMgr->Info()<<"StV0FinderMaker : coming in I have "<<v0Vertices.size()<<" V0s."<<endm;

  if (!(1&useV0Language)) {
    //Erase existing V0s and Xis
    // (must do Xis too as they point to the V0s!)
    gMessMgr->Info()<<"StV0FinderMaker : pre-existing V0s and Xis deleted."<<endm; 
    StSPtrVecV0Vertex v0Vertices2;
    v0Vertices = v0Vertices2;
    StSPtrVecXiVertex& xiVertices   = event->xiVertices();
    StSPtrVecXiVertex xiVertices2;
    xiVertices = xiVertices2;
  }


  // Loop over track pairs to find V0s

  //i track is positive
  for (ii=0; ii<ptrks; ii++) {
    i = ptrk[ii];

    xci.Set(heli[i].xcenter(),heli[i].ycenter());
    ri.Set(heli[i].origin().x(),heli[i].origin().y());
    ri -= xci;
    rad_i = ri.Mod();

    //j track is negative
    for (jj=0; jj<ntrks; jj++) {
      j = ntrk[jj];

      if (GetTrackerUsage() == kTrackerUseBOTH)
         {if ((trk[i]->fittingMethod() == ITTFflag) && (trk[j]->fittingMethod() != ITTFflag)) continue;
          if ((trk[i]->fittingMethod() != ITTFflag) && (trk[j]->fittingMethod() == ITTFflag)) continue;
          }

      // Determine detector id of V0 for pars
      det_id_v0 = TMath::Max(detId[i],detId[j]);

      // Primary   V0 cut parameters
      pars  = ev0par2->GetTable(det_id_v0+2);
      // Secondary V0 cut parameters
      pars2 = ev0par2->GetTable(det_id_v0-1);

      //Cut: number of hits
      if ((hits[i] < pars2->n_point) ||
          (hits[j] < pars2->n_point)) continue;

      //Cut: Initial cut on dca of tracks to primary vertex
      // (perform as early as possible)
      // V0 can't have pt larger than sum of pts of daughters
      temp = pt[i] + pt[j];
      if ((temp*temp < 0.98*ptV0sq) &&
          ((trk[i]->impactParameter() <= pars2->dcapnmin) ||
           (trk[j]->impactParameter() <= pars2->dcapnmin))) continue;

      xcj.Set(heli[j].xcenter(),heli[j].ycenter());
      rj.Set(heli[j].origin().x(),heli[j].origin().y());
      rj -= xcj;
      rad_j = rj.Mod();

      tmp2V = xci - xcj;
      dxc = tmp2V.Mod();
      separation = dxc - (rad_i + rad_j);
      doSecond = kFALSE;
      dca_ij1 = -9999;


      // ********************* START OF DETERMINATION OF V0 GEOMETRY
      if (separation < 0)
         {// Check for one helix circle completely inside the other
          if (dxc < TMath::Abs(rad_i - rad_j)) continue;
          // Helix circles are overlapping
          //FULL 3D ITERATIVE METHOD (VERY SLOW, ONE SOLUTION)
          //      paths = heli[i].pathLengths(heli[j]);
          //2D+ METHOD (GETS 3D APPROXIMATION AFTER TWO 2D SOLUTIONS)
          path2 = heli[i].pathLength(rad_j,xcj.X(),xcj.Y());
          // Two possible solutions: process ones that aren't nans
          if (!isnan(path2.first))
             {solution = path2.first;
              if ((!isnan(path2.second)) && (path2.second != path2.first))
                 {doSecond = kTRUE;
                  }
              goto ProcessSolution;
              }
             else if (isnan(path2.second))
             {// no solutions
              continue;
              }
              //else run only with the second solution
          SecondSolution:
          solution = path2.second;
          doSecond = kFALSE;
          ProcessSolution:
          // paths contains the pathlengths for this solution with
          // that for track i stored in first, and track j stored
          // in second.
          paths.first = solution;
          xi = heli[i].at(paths.first );
          paths.second = heli[j].pathLength(xi.x(),xi.y());
          xj = heli[j].at(paths.second);
          }
          else if (separation < pars2->dca)
         {// Helix circles are close, but not overlapping,
          // find dca to point halfway between circle centers
          tmp2V = (xci + xcj) * 0.5;
          paths.first  = heli[i].pathLength(tmp2V.X(),tmp2V.Y());
          paths.second = heli[j].pathLength(tmp2V.X(),tmp2V.Y());
          xi = heli[i].at(paths.first );
          xj = heli[j].at(paths.second);
          }
          else
         {// Helix circles are too far apart
          continue;
          }
      
      dca_ij = xi.z() - xj.z();
      if (doSecond) {
        // If we have two solutions, save this one and compare
        dca_ij1 = dca_ij;
        xi1=xi;
        xj1=xj;
        paths1=paths;
        goto SecondSolution;
        }
      if ((dca_ij1 != -9999) &&
          (TMath::Abs(dca_ij1) < TMath::Abs(dca_ij))) {
        // First solution was better
        dca_ij = dca_ij1;
        xi=xi1;
        xj=xj1;
        paths=paths1;
        }
      // At this point, dca_ij is *signed* for use in 3D calc
      // *********************  END  OF DETERMINATION OF V0 GEOMETRY

      pi = heli[i].momentumAt(paths.first ,Bfield);
      pj = heli[j].momentumAt(paths.second,Bfield);
      
      //Cut: check if tracks points away from prim vtx
      if ((pi.dot(xi-mainv) < 0.0) ||
          (pj.dot(xj-mainv) < 0.0)) continue;
						
						///Begin Betty 
      //Cut: check if the first point of either track is after v0vertex
      if ((pi.dot(heli[i].origin() - xi) < 0.0) ||  //if pV0 * r <0, cut
         (pj.dot(heli[j].origin() - xj) < 0.0)) continue;
      ///End Betty


      // ********************* START OF DETERMINATION OF 3D DCA
      // dca_ij will be an approximation of the 3D dca
      pi1 = pi/ptot[i];
      pj1 = pj/ptot[j];

      cosij = pi1.dot(pj1);
      sin2ij = 1.0 - cosij*cosij;
      if (sin2ij) {
        temp = dca_ij/sin2ij;
        t1 = (-pi1.z()+pj1.z()*cosij)*temp;
        t2 = ( pj1.z()-pi1.z()*cosij)*temp;
 
        temp = rad_i*(ptot[i]/pt[i]);
        temp *= sin(t1/temp);
        xi1 = xi + pi1.pseudoProduct(temp,temp,t1);

        temp = rad_j*(ptot[j]/pt[j]);
        temp *= sin(t2/temp);
        xj1 = xj + pj1.pseudoProduct(temp,temp,t2);

        dca_ij1 = (xi1 - xj1).mag2();
        dca_ij *= dca_ij; /// dca_ij no longer signed (squared)

        if (dca_ij1 < dca_ij) {
          paths.first  = heli[i].pathLength(xi1.x(),xi1.y());
          paths.second = heli[j].pathLength(xj1.x(),xj1.y());
          ///paths.first  = heli[i].pathLength(xnix,xniy);
          ///paths.second = heli[j].pathLength(xnjx,xnjy);
          xi1 = heli[i].at(paths.first);
          xj1 = heli[j].at(paths.second);
          dca_ij1 = (xi1 - xj1).mag2();
          if (dca_ij1 < dca_ij) {
            xi = xi1;
            xj = xj1;
            pi = heli[i].momentumAt(paths.first ,Bfield);
            pj = heli[j].momentumAt(paths.second,Bfield);
            dca_ij = dca_ij1;
          }
        }
        //This code is the new one.

      /*if (dca_ij1 < dca_ij) {
          paths.first  = heli[i].pathLength(xi1.x(),xi1.y());
          paths.second = heli[j].pathLength(xj1.x(),xj1.y());
          xi = xi1;
          xj = xj1;
          pi = heli[i].momentumAt(paths.first ,Bfield);
          pj = heli[j].momentumAt(paths.second,Bfield);
          dca_ij = dca_ij1;
          }*/
        //And this one is the old one (comparison with Fortran).

 
      }
      // *********************  END  OF DETERMINATION OF 3D DCA


      // Now we have the positions and momenta of our tracks
      // at the V0 vertex determined. Ready to make cuts on
      // the V0 itself. 

      //Cut: dca between tracks      
      if (dca_ij >= (pars2->dca*pars2->dca)) continue;

      pp = pi + pj;
      pperpsq = pp.perp2();
                  
      //Cut: dca of tracks to primary vertex (early as possible)
      if ((pperpsq < ptV0sq) && 
          ((trk[i]->impactParameter() <= pars2->dcapnmin) ||
           (trk[j]->impactParameter() <= pars2->dcapnmin))) continue;

      xpp = (xi + xj) * 0.5;
      impact = xpp - mainv;
      dlen = impact.mag2();

      //Cut: decay length from prim vtx
      if (dlen <= (pars2->dlen*pars2->dlen)) continue;

      //Cut: V0 momentum should be away from prim vtx 
      if (pp.dot(impact) < 0.0) continue;

      ppmag = pperpsq + pp.z()*pp.z();
      rmin = impact.cross(pp).mag2()/ppmag;

      //Cut: dca of V0 to prim vtx
      if (rmin >= (pars2->dcav0*pars2->dcav0)) continue;

      tmp3V = pp/::sqrt(ppmag);
      pPosAlongV0 = pi.dot(tmp3V);
      pNegAlongV0 = pj.dot(tmp3V);
      alpha = (pPosAlongV0-pNegAlongV0) /
              (pPosAlongV0+pNegAlongV0);

      //Cut: Armenteros alpha
      if (TMath::Abs(alpha) > pars2->alpha_max) continue;

      ptArm_sq = ptot[i]*ptot[i] - pPosAlongV0*pPosAlongV0;

      //Cut: Armenteros pt
      if (ptArm_sq > (pars2->ptarm_max*pars2->ptarm_max)) continue;

      rmin   = ::sqrt(rmin);
      dca_ij = ::sqrt(dca_ij);
      if (trk[i]->fittingMethod() == ITTFflag) dca_ij=-dca_ij;
      
      // Fill an StV0Vertex
      v0Vertex = new StV0Vertex();
      v0Vertex->setPosition(xpp);
      v0Vertex->addDaughter(trk[i]);
      v0Vertex->addDaughter(trk[j]);
      v0Vertex->setDcaDaughterToPrimaryVertex(positive,trk[i]->impactParameter());
      v0Vertex->setDcaDaughterToPrimaryVertex(negative,trk[j]->impactParameter());
      ///3VectorF vs 3VectorD???
      v0Vertex->setMomentumOfDaughter(positive,pi);
      v0Vertex->setMomentumOfDaughter(negative,pj);
      v0Vertex->setDcaDaughters(dca_ij);
      v0Vertex->setDcaParentToPrimaryVertex(rmin);
      
      ///Begin Betty
      /*Set bits 0000 0000 -> 1001 0000 if C++ used
      if SVT ran in tracking, 1001 1000;
                   SVT and +: 1001 1100
                SVT, + and -: 1001 1110; 
                 SVT, - only: 1001 1010*/
      
      ///variable keepTrack (32 bits)
      keepTrack=0;
      keepTrack |=((long)1 << 4); // c++ v0 finder ran
      if (useSVT)
         {keepTrack |=((long)1 << 3);  //sets to one the fourth digit if SVT was used
          if (detId[i]==2 || detId[i]==3){
             keepTrack |=((long)1 << 2); //sets 3rd bit to 1 if + track came from SVT  
             }
          if(detId[j]==2 || detId[j]==3){
             keepTrack |=((long)1 << 1); //sets 2nd bit to 1 if - track came from SVT
             }
          }
      keepTrack *= -1; //sets to negative the last digit 
      v0Vertex->setChiSquared((float)keepTrack);
      ///End Betty
      
      // Use primary V0 cut parameters
      isPrimaryV0 =
        (rmin < pars->dcav0) &&
        ((pperpsq >= ptV0sq) ||
         ((trk[i]->impactParameter() > pars->dcapnmin) &&
          (trk[j]->impactParameter() > pars->dcapnmin)));
      
      // Call secondary usage of V0 (such as for Xis)
      usedV0 = UseV0();
      
      // Tag used V0s to indicate if they aren't primary
      if (usedV0 && !isPrimaryV0)
        v0Vertex->setDcaParentToPrimaryVertex(-rmin);

      // If used or primary, keep it
      if (usedV0 || isPrimaryV0) {
        v0Vertices.push_back(v0Vertex);
      } else {
        delete v0Vertex;
        v0Vertex = 0;
      }

    } // j-Loop
  } // i-Loop

  gMessMgr->Info()<<"StV0FinderMaker : now I have "<<v0Vertices.size()<<" V0s."<<endm;
  gMessMgr->Info()<<"StV0FinderMaker : using magnetic field : "<<Bfield/tesla<<" T."<<endm;

  // Any cleanup involved for using KeepV0()
  
  return kStOk;
}














//_____________________________________________________________________________
void StV0FinderMaker::Trim() {
  // Loop over V0s and remove those that don't satisfy the tight V0 cuts

  gMessMgr->Info() << "StV0FinderMaker::Trim() : Starting..." << endm;

  event = (StEvent*) GetInputDS("StEvent");
  pars = ev0par2->GetTable(3);
  StSPtrVecV0Vertex& v0Vertices = event->v0Vertices();
  int iV0s = v0Vertices.size();
  for (int i=iV0s-1; i>=0; i--) {

    v0Vertex = v0Vertices[i];
    if ((v0Vertex) &&
        // Is it not a seconday V0?
        (v0Vertex->dcaParentToPrimaryVertex() >= 0) &&
        // Is it not a primary V0?
        ! ((v0Vertex->dcaParentToPrimaryVertex() < pars->dcav0) &&
           ((v0Vertex->momentum().perp2() >= ptV0sq) ||
            ((v0Vertex->dcaDaughterToPrimaryVertex(positive) > pars->dcapnmin) &&
             (v0Vertex->dcaDaughterToPrimaryVertex(negative) > pars->dcapnmin))))) {
      v0Vertex->makeZombie();
      iV0s--;
    }
  } // V0 loop

  gMessMgr->Info() << "StV0FinderMaker::Trim() : saving " << iV0s <<
                      " V0 candidates" << endm;
}
//_____________________________________________________________________________
// $Id$
// $Log$
// Revision 1.10  2003/08/22 17:47:14  caines
// Get sign AND magnitude of mag field correctly for Xi and V0 finder
//
// Revision 1.9  2003/07/17 17:01:10  faivre
// Add one causality check. Exhaustive listing of cuts. Tab->spaces.
//
// Revision 1.8  2003/07/15 17:40:36  faivre
// Fixed charge of Bfield (used for print, and XiFinder since now).
//
// Revision 1.7  2003/07/04 17:52:54  faivre
// Use SVT cuts if any dg has a SVT hit.
//
// Revision 1.6  2003/06/24 16:20:01  faivre
// Uses SVT tracks. Fixed bool calculations. Exits when bad param. Reshaping.
//
// Revision 1.5  2003/05/14 19:15:03  faivre
// Fancy choices Fortran/C++ V0's and Xi's. SVT tracks.
//
// Revision 1.4  2003/05/02 21:21:08  lbarnby
// Now identify ITTF tracks by fittingMethod() equal to  kITKalmanFitId
//
// Revision 1.3  2003/04/30 20:38:22  perev
// Warnings cleanup. Modified lines marked VP
//
// Revision 1.2  2003/04/30 19:14:27  faivre
// ITTF vs TPT V0s
//
//
