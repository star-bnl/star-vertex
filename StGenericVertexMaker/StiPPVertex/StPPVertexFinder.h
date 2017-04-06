/*!
 * \class StPPVertexFinder
 * \author Jan Balewski, July 2004
 *
 *  StGenericVertexFinder implementation of PPV
 * $Id: StPPVertexFinder.h,v 1.53 2017/05/24 05:02:06 genevb Exp $
 *
 */

#ifndef StPPVertexFinder_h
#define StPPVertexFinder_h

#include <vector>

#include "StGenericVertexMaker/StGenericVertexFinder.h"
#include "StGenericVertexMaker/StiPPVertex/TrackData.h"
#include "StGenericVertexMaker/StiPPVertex/VertexData.h"

#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StMuDSTMaker/COMMON/StMuDst.h"

class TH1F;
class TH2F;
class TH1D;

class StiKalmanTrack;
class StEvent; 
class StiToolkit;


class St_db_Maker;
class BtofHitList;  
class CtbHitList;
class BemcHitList;
class EemcHitList;


template<class Event_t, class Track_t>
class StPPVertexFinderT: public StGenericVertexFinder
{
 protected:

  /// Takes a list of vertex candidates/seeds and updates each vertex position
  /// by fitting tracks pointing to it
  int fitTracksToVertex(VertexData &vertex);

  /// Creates DCA states for selected tracks (mTrackData) and fills the member
  /// container mDCAs
  void fillTrackDcas(const VertexData &vertex);

  /// using the ROOT's TSpectrum peak finder applied to the distribution of
  /// track DCAs along the `z` axis
  void findSeeds_TSpectrum();

  /// Searches for vertex candidates and fills private `mVertexData` container
  /// by building a likelihood distribution of track DCAs along the `z` axis
  void findSeeds_PPVLikelihood();

  enum {mxH=32};

  void matchTrack2BTOF(Track_t &track);
  void matchTrack2BTOF(const StPhysicalHelixD& helix, Track_t &track);

  void matchTrack2EEMC(Track_t &track);
  void matchTrack2EEMC(const StPhysicalHelixD& helix, Track_t &track);

  void matchTrack2BEMC(Track_t &track);
  void matchTrack2BEMC(const StPhysicalHelixD& helix, Track_t &track);

  void matchTrack2BTOF(Track_t &track);
  void matchTrack2BTOF(const StPhysicalHelixD& helix, Track_t &track);

  bool matchTrack2Membrane(Track_t &track);

  bool buildLikelihoodZ();
  bool findVertexZ(VertexData &);
  bool evalVertexZ(VertexData &);
  void exportVertices(); 
  void saveHisto(TString fname);

  /// A container with pre-selected tracks to be used in seed finding
  std::vector<Track_t>  mTrackData;
  std::vector<VertexData> mVertexData;
  int  mTotEve;
  int  eveID;
  int  nBadVertex;
  unsigned int  mAlgoSwitches; ///< binary, assign 1bit per change, use enum below
                               ///< default is 0, as for 2008 pp data production
  enum {kSwitchOneHighPT=1}; 

  TH1D *hL ;      // likelyhood distribution
  TH1D *hM, *hW ; // cumulative track mult & weight distribution, for better errZ calculation
  std::array<int, 8> ntrk;

  // params
  float  mMinZBtof;               ///< BTOF local z min cut
  float  mMaxZBtof;               ///< BTOF local z max cut
  float  mMinAdcBemc;             ///< BEMC towers with MIP response
  float  mMinAdcEemc;             ///< EEMC towers with MIP response

  /// A flag whether to use nFit/nPossible in track weighting (ranking).
  /// Introduced in 2012 for pp510 to differentiate between global track
  /// quality, together with lowering the overall threshold from 0.7 to 0.51.
  /// Set to false prior to 2012, true thereafter
  bool   mFitPossWeighting;

  bool   mDropPostCrossingTrack;  ///< enable/disable post crossing tarck rejection
  int    mStoreUnqualifiedVertex; ///< set the max # of vertices, sorted by rank
  float  mCut_oneTrackPT;         ///< threshold for storing one track vertices.
                                  ///< In GeV, used only if coresponding algoSwitch switch is ON.
  bool   mUseBTOFmatchOnly;        ///< enable/disable using only TOF-matched tracks

  StiToolkit     *mToolkit;
  BtofHitList    *btofList;
  CtbHitList     *ctbList;
  BemcHitList    *bemcList;
  EemcHitList    *eemcList;

  /// A pointer to muDST event
  const StMuDst* mStMuDst;
  
  /// A helper function to do common processing for StEvent and StMuDst cases
  void seed_fit_export();

  virtual void  UseVertexConstraint() {}

  virtual void UpdateVertexCuts(int run_number);

public:

  virtual void numVerticesToStore(int n) { mStoreUnqualifiedVertex = n; }

  virtual void UsePCT(bool x=true) { mDropPostCrossingTrack = !x; }
  virtual void Finish();
  virtual void Init();
  virtual void InitRun(int run_number, const St_db_Maker* db_maker);
  virtual void Clear(); 

  StPPVertexFinderT(VertexFit_t fitMode=VertexFit_t::BeamlineNoFit);

  virtual ~StPPVertexFinderT() {}
  virtual int fit(StEvent*) { return -1; }
  virtual int fit(const Event_t& event) { return -1; }
  virtual void SetStoreUnqualifiedVertex(int n) { mStoreUnqualifiedVertex = n; }
  virtual void UseBTOFmatchOnly(bool useBTOFmatchOnly = true) { UseBTOF(); mUseBTOFmatchOnly = useBTOFmatchOnly; }

  virtual void printInfo(std::ostream& os = std::cout) const;
};



template<class Event_t>
class StPPVertexFinder : public StPPVertexFinderT<Event_t, void> { };



template<>
class StPPVertexFinder<StEvent> : public StPPVertexFinderT<StEvent, TrackData<StiKalmanTrack> >
{
public:

   using Track_t = TrackData<StiKalmanTrack>;

   StPPVertexFinder(VertexFit_t fitMode=VertexFit_t::BeamlineNoFit) :
      StPPVertexFinderT(fitMode) { };

   virtual int fit(const StEvent& event);

   bool examinTrackDca(Track_t &track);
   void matchTrack2CTB(Track_t &track);
   bool isPostCrossingTrack(const StiKalmanTrack* stiTrack);
   void dumpKalmanNodes(const StiKalmanTrack& stiTrack);
};



template<>
class StPPVertexFinder<StMuDst> : public StPPVertexFinderT<StMuDst, TrackData<StMuTrack> >
{
public:

   using Track_t = TrackData<StMuTrack>;

   StPPVertexFinder(VertexFit_t fitMode=VertexFit_t::Beamline3D) :
      StPPVertexFinderT(fitMode) { };

   virtual int fit(const StMuDst& event);

   virtual void result(TClonesArray& stMuPrimaryVertices, TClonesArray& stMuPrimaryTracks);
};



#endif
