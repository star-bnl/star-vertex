/*!
 * \class StPPVertexFinder
 * \author Jan Balewski, July 2004
 *
 *  StGenericVertexFinder implementation of PPV
 * $Id: StPPVertexFinder.h,v 1.44 2017/03/04 04:50:21 smirnovd Exp $
 *
 */

#include <vector>

#include "StGenericVertexMaker/StGenericVertexFinder.h"
#include "StGenericVertexMaker/StiPPVertex/TrackData.h"
#include "StGenericVertexMaker/StiPPVertex/VertexData.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"

class TH1F;
class TH2F;
class TH1D;

class StiKalmanTrack;
class StEvent; 
class StiToolkit;
class StEEmcDb;

class StMuDst;
class EEmcGeomSimple;
class StMuTrack;
class StBTofGeometry; 

class St_db_Maker;
class BtofHitList;  
class CtbHitList;
class BemcHitList;
class EemcHitList;


class StPPVertexFinder: public StGenericVertexFinder
{
 private:

  /// Takes a list of vertex candidates/seeds and updates each vertex position
  /// by fitting tracks pointing to it
  int fitTracksToVertex(VertexData &vertex);

  /// Creates DCA states for selected tracks (mTrackData) and fills the member
  /// container mDCAs
  void createTrackDcas(const VertexData &vertex);

  /// using the ROOT's TSpectrum peak finder applied to the distribution of
  /// track DCAs along the `z` axis
  void findSeeds_TSpectrum();

  /// Searches for vertex candidates and fills private `mVertexData` container
  /// by building a likelihood distribution of track DCAs along the `z` axis
  void findSeeds_PPVLikelihood();

  enum {mxH=32};
  bool examinTrackDca(const StiKalmanTrack*, TrackData &track);
  void matchTrack2BTOF(const StiKalmanTrack*, TrackData &track, StBTofGeometry *geom);
  void matchTrack2CTB(const StiKalmanTrack*, TrackData &track);

  void matchTrack2EEMC(const StiKalmanTrack*, TrackData &track);
  void matchTrack2EEMC(const StMuTrack& muTrack, TrackData &track);
  void matchTrack2EEMC(const StPhysicalHelixD& helix, TrackData &track);

  void matchTrack2BEMC(const StiKalmanTrack*, TrackData &track);
  void matchTrack2BEMC(const StMuTrack& muTrack, TrackData &track);
  void matchTrack2BEMC(const StPhysicalHelixD& helix, TrackData &track);

  bool matchTrack2Membrane(const StiKalmanTrack*, TrackData &track);
  void matchTrack2Membrane(const StMuTrack& muTrack, TrackData &track);

  bool isPostCrossingTrack(const StiKalmanTrack* stiTrack);

  bool buildLikelihoodZ();
  bool findVertexZ(VertexData &);
  bool evalVertexZ(VertexData &);
  void exportVertices(); 
  void saveHisto(TString fname);

  /// A container with pre-selected tracks to be used in seed finding
  std::vector<TrackData>  mTrackData;
  std::vector<VertexData> mVertexData;
  int  mTotEve;
  int  eveID;
  int  nBadVertex;
  unsigned int  mAlgoSwitches; ///< binary, assign 1bit per change, use enum below
                               ///< default is 0, as for 2008 pp data production
  enum {kSwitchOneHighPT=1}; 

  TH1F *hA[mxH];
  TH2F *hACorr;
  TH1D *hL ;      // likelyhood distribution
  TH1D *hM, *hW ; // cumulative track mult & weight distribution, for better errZ calculation
  TObjArray HList;

  // params
  double mMinTrkPt;               ///< ~ pT=0.16(GeV/c) == R=2 (m )in 2001
  double mMaxTrkDcaRxy;           ///< DCA to nominal beam line for each track
  float  mMaxZradius;             ///< used in matching: tracks to zVertex
  int    mMinMatchTr;             ///< for valid vertex
  float  mMaxZrange;              ///< cut off for tracks Z_DCA
  float  mMinZBtof;               ///< BTOF local z min cut
  float  mMaxZBtof;               ///< BTOF local z max cut
  float  mMinAdcBemc;             ///< BEMC towers with MIP response
  float  mMinAdcEemc;             ///< EEMC towers with MIP response
  float  mMinFitPfrac;            ///< nFit/nPossible

  /// A flag whether to use nFit/nPossible in track weighting (ranking).
  /// Introduced in 2012 for pp510 to differentiate between global track
  /// quality, together with lowering the overall threshold from 0.7 to 0.51.
  /// Set to false prior to 2012, true thereafter
  bool   mFitPossWeighting;

  bool   mDropPostCrossingTrack;  ///< enable/disable post crossing tarck rejection
  int    mStoreUnqualifiedVertex; ///< set the max # of vertices, sorted by rank
  float  mCut_oneTrackPT;         ///< threshold for storing one track vertices.
                                  ///< In GeV, used only if coresponding algoSwitch switch is ON.

  StiToolkit     *mToolkit;
  BtofHitList    *btofList;
  CtbHitList     *ctbList;
  BemcHitList    *bemcList;
  EemcHitList    *eemcList;
  StBTofGeometry *btofGeom;
  EEmcGeomSimple *geomE;

  const StMuDst* mStMuDst;
  
  void dumpKalmanNodes(const StiKalmanTrack *stiTrack);
  void initHisto();

  virtual void  UseVertexConstraint() {}

public:

  virtual void UsePCT(bool x=true) { mDropPostCrossingTrack = !x; }
  virtual void Finish();
  virtual void Init();
  virtual void InitRun(int runumber, const St_db_Maker* db_maker);
  virtual void Clear(); 
  virtual void CalibBeamLine(); // activates saving high quality prim tracks for 3D fit of the beamLine

  StPPVertexFinder(VertexFit_t fitMode=VertexFit_t::Beamline1D);

  virtual ~StPPVertexFinder() {}
  virtual int fit(StEvent*);
  virtual int Fit(const StMuDst& muDst);
  void printInfo(std::ostream& = std::cout) const;
};
