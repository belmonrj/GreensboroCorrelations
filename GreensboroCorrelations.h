#ifndef __GREENSBOROCORRELATIONS_H__
#define __GREENSBOROCORRELATIONS_H__

// standard includes
#include <string>
#include <vector>
#include <SubsysReco.h>
#include "TComplex.h"


class dAuBES_utils;
class TFile;
//class TTree;
class TH1D;
class TH2D;
class TProfile;


class GreensboroCorrelations: public SubsysReco
{
 public:
  GreensboroCorrelations();
  virtual ~GreensboroCorrelations();

  /// Fun4All calls...
  int  Init         (PHCompositeNode *topNode);
  int  InitRun      (PHCompositeNode *topNode);
  int  process_event(PHCompositeNode *topNode);
  //  int  ResetEvent   (PHCompositeNode *topNode);
  int  End          (PHCompositeNode *topNode);
  int  EndRun       (PHCompositeNode *topNode);
  void Verbosity    (int verbosity) {_verbosity = verbosity;}

  /// Single particle ntuple output...
  void set_output_filename(std::string filename) { _output_filename = filename; } // select output file name externally
  void SetQvectorOffsets(int runnumber);
  void SetQvectorOffsetsRBR(int runnumber);
  void set_do_double_track_cut(bool b){do_double_track_cut = b;}
  void set_zvtxcut(double z){_cut_zvtx = z;}
  void set_chi2cut(double c){_cut_chi2 = c;}
  void set_dcacut(double d){_cut_dca = d;}
  void set_nhitcut(int n){_cut_nhit = n;}

 protected:

  // --- do the analysis
  int EventRecursion();
  int EventStuff();

  // --- special event cuts
  bool PassesTracksChargeRatio(int, double);

  // --- cumulants functions
  float calc2_event(float, float, float);
  float calc4_event(float, float, float, float, float);
  float calc6_event(TComplex&, TComplex&, TComplex&, float);
  // --- acceptance correction functions
  float calccossum2_event(TComplex&, TComplex&, float);
  float calcsinsum2_event(TComplex&, TComplex&, float);
  float calccos3_event(TComplex&, TComplex&, float);
  float calcsin3_event(TComplex&, TComplex&, float);

  //static const int h1=2;
  //static const int h2=-2;
  //static const int h3=2;
  //static const int h4=-2;
  //static const int h5=2;
  //static const int h6=-2;
  //static const int h7=2;
  //static const int h8=-2; // simplifed to v2
  //static const int sum = (h1<0?-1*h1:h1)+(h2<0?-1*h2:h2)+(h3<0?-1*h3:h3)+(h4<0?-1*h4:h4)
  //  + (h5<0?-1*h5:h5)+(h6<0?-1*h6:h6)+(h7<0?-1*h7:h7)+(h8<0?-1*h8:h8);
  //static const int maxCorrelator = 8; // We will not go beyond 8-p correlations
  //static const int maxHarmonic = sum+1;
  //static const int maxPower = maxCorrelator+1;
  static const int maxCorrelator = 12; // Somewhat abusing the setup as it is...
  static const int maxHarmonic = 10; // Need to assess on case-by-case basis, but this gets you v2{8} and v3{6}
  static const int maxPower = 9;
  TComplex Qvector[maxHarmonic][maxPower]; // All needed Q-vector components
  TComplex Qvector_north[maxHarmonic][maxPower];
  TComplex Qvector_south[maxHarmonic][maxPower];
  TComplex Qoffset[maxHarmonic][maxPower];
  TComplex Qoffset_north[maxHarmonic][maxPower];
  TComplex Qoffset_south[maxHarmonic][maxPower];
  TComplex Q(int, int);
  TComplex Recursion(int, int*);
  TComplex Recursion(int, int*, int, int);

  // ---------------------------------------------------------------------
  static const int maxTracks = 650; // accept no more FVTX tracks than this

  static const double default_cut_zvtx = 10.0;
  static const double default_cut_chi2 = 5.0;
  static const double default_cut_dca = 2.0;
  static const int default_cut_nhit = 3;

  std::vector<double> fphi;
  std::vector<double> feta;
  std::vector<double> fdcax;
  std::vector<double> fdcay;
  std::vector<double> fchi2ndf;
  std::vector<int> fnhitspc;

  bool fvtx_track_passes[maxTracks];

  /// current event
  unsigned long _ievent;

  /// verbosity level
  int _verbosity;

  /// module output filename
  std::string _output_filename;
  TFile* _output_file;

  // whether to do the double track cut
  bool do_double_track_cut;

  double _cut_zvtx;
  double _cut_chi2;
  double _cut_dca;
  int _cut_nhit;



  dAuBES_utils* _utils;            ///< Utilities class
  bool use_utils;
  std::string _collsys;

  int tmp_evt;

  static const int nharm = 5;
  float d_SouthQX[nharm];
  float d_SouthQY[nharm];
  float d_SouthQW;
  float d_NorthQX[nharm];
  float d_NorthQY[nharm];
  float d_NorthQW;

  // --- Q-vector offset variables
  double qvoff_nfvtxt[maxTracks][2][maxHarmonic];
  double qvoff_nfvtxt_north[maxTracks][2][maxHarmonic];
  double qvoff_nfvtxt_south[maxTracks][2][maxHarmonic];
  double qvoff_cent[100][2][maxHarmonic];
  double qvoff_cent_north[100][2][maxHarmonic];
  double qvoff_cent_south[100][2][maxHarmonic];



  // ------------------------------------------
  // --- Histograms
  // ------------------------------------------

  TH1D* th1d_nfvtxt_combinedER;
  TH1D* th1d_nfvtxt_combined;
  TH1D* th1d_nfvtxt_north;
  TH1D* th1d_nfvtxt_south;
  TH2D* th2d_nfvtxt_northsouth;
  TH2D* th2d_nfvtxt_bbcsum;
  TH2D* th2d_nfvtxt_centrality;
  TH2D* th2d_nfvtxt_centralityA;
  TH2D* th2d_nfvtxt_bbcsumratio;
  TH1D* th1d_centrality;
  TH1D* th1d_centralityA;
  TH1D* th1d_track_deta;
  TH1D* th1d_track_dphi;
  TH2D* th2d_track_before_eta;
  TH2D* th2d_track_before_phi;
  TH2D* th2d_track_after_eta;
  TH2D* th2d_track_after_phi;
  TH2D* th2d_track_aafter_eta;
  TH2D* th2d_track_aafter_phi;
  TH2D* th2d_cent_dcax;
  TH2D* th2d_cent_dcay;
  TH2D* th2d_cent_nhitr;
  TH2D* th2d_cent_nhits;
  TH2D* th2d_cent_chisq;
  TProfile* tp1f_track_detacutpass;

  // --- centrality stuff




  // --- event plane decorrelation stuff





};

#endif /* __GREENSBOROCORRELATIONS_H__ */



