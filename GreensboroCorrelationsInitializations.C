#include <GreensboroCorrelations.h>

#include <iostream>

#include <TFile.h>
//#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>

#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <RunHeader.h>
#include <dAuBES_utils.h>




using namespace std;


// --- Init method, part of Fun4All inheriance
int GreensboroCorrelations::Init(PHCompositeNode *topNode)
{

  //  ResetEvent(topNode); // is this needed?

  if (_verbosity > 1) cout << PHWHERE << "::Init() - entered." << endl;

  _output_file = new TFile(_output_filename.c_str(), "RECREATE");

  // ---
  // --- initialize histograms
  // ---

  th1d_cnt_both_phi = new TH1D("th1d_cnt_both_phi","",640,-3.2,3.2);
  th1d_cnt_east_phi = new TH1D("th1d_cnt_east_phi","",640,-3.2,3.2);
  th1d_cnt_west_phi = new TH1D("th1d_cnt_west_phi","",640,-3.2,3.2);
  th1d_cnt_both_phi_high = new TH1D("th1d_cnt_both_phi_high","",640,-3.2,3.2);
  th1d_cnt_east_phi_high = new TH1D("th1d_cnt_east_phi_high","",640,-3.2,3.2);
  th1d_cnt_west_phi_high = new TH1D("th1d_cnt_west_phi_high","",640,-3.2,3.2);

  th1d_nfvtxt_combinedER = new TH1D("th1d_nfvtxt_combinedER","",5000, -0.5, 4999.5);
  th1d_nfvtxt_combined = new TH1D("th1d_nfvtxt_combined","",700, -0.5, 699.5);
  th1d_centrality = new TH1D("th1d_centrality","",100, -0.5, 99.5);
  th1d_centralityA = new TH1D("th1d_centralityA","",100, -0.5, 99.5);
  th2d_nfvtxt_bbcsum = new TH2D("th2d_nfvtxt_bbcsum","",700, -0.5, 699.5, 1000, 0, 4000);
  th2d_nfvtxt_centrality = new TH2D("th2d_nfvtxt_centrality","",700, -0.5, 699.5, 100, -0.5, 99.5);
  th2d_nfvtxt_centralityA = new TH2D("th2d_nfvtxt_centralityA","",700, -0.5, 699.5, 100, -0.5, 99.5);
  th2d_nfvtxt_bbcsumratio = new TH2D("th2d_nfvtxt_bbcsumratio","",700, -0.5, 699.5, 1000, 0, 5);
  th1d_nfvtxt_north = new TH1D("th1d_nfvtxt_north","",700, -0.5, 699.5);
  th1d_nfvtxt_south = new TH1D("th1d_nfvtxt_south","",700, -0.5, 699.5);
  th2d_nfvtxt_northsouth = new TH2D("th2d_nfvtxt_northsouth","",1000, -0.5, 999.5, 1000, -0.5, 999.5);
  th1d_track_deta = new TH1D("th1d_track_deta","",2000, -0.1, 0.1);
  th1d_track_dphi = new TH1D("th1d_track_dphi","",2000, -0.1, 0.1);
  th2d_track_before_eta = new TH2D("th2d_track_before_eta","",100,-0.5,99.5, 700, -3.5, 3.5);
  th2d_track_before_phi = new TH2D("th2d_track_before_phi","",100,-0.5,99.5, 640, -3.2, 3.2);
  th2d_track_after_eta  = new TH2D("th2d_track_after_eta", "",100,-0.5,99.5, 700, -3.5, 3.5);
  th2d_track_after_phi  = new TH2D("th2d_track_after_phi", "",100,-0.5,99.5, 640, -3.2, 3.2);
  th2d_track_aafter_eta = new TH2D("th2d_track_aafter_eta","",100,-0.5,99.5, 700, -3.5, 3.5);
  th2d_track_aafter_phi = new TH2D("th2d_track_aafter_phi","",100,-0.5,99.5, 640, -3.2, 3.2);
  th2d_cent_dcax = new TH2D("th2d_cent_dcax","",100,-0.5,99.5, 600,-3.0,3.0);
  th2d_cent_dcay = new TH2D("th2d_cent_dcay","",100,-0.5,99.5, 600,-3.0,3.0);
  th2d_cent_nhitr = new TH2D("th2d_cent_nhitr","",100,-0.5,99.5, 7,-0.5,6.5);
  th2d_cent_nhits = new TH2D("th2d_cent_nhits","",100,-0.5,99.5, 7,-0.5,6.5);
  th2d_cent_chisq = new TH2D("th2d_cent_chisq","",100,-0.5,99.5, 100,0.0,10.0);
  tp1f_track_detacutpass = new TProfile("tp1f_track_detacutpass","",100,-0.5,99.5,-0.1,1.1);



  th1d_fvtxs_eta = new TH1D("th1d_fvtxs_eta","",20,-3.0,-1.0);
  th1d_fvtxs_phi = new TH1D("th1d_fvtxs_phi","",20,-pi/4,3*pi/4);
  th1d_fvtxn_eta = new TH1D("th1d_fvtxn_eta","",20,1.0,3.0);
  th1d_fvtxn_phi = new TH1D("th1d_fvtxn_phi","",20,-pi/4,3*pi/4);
  th1d_fvtxs_gap_phi = new TH1D("th1d_fvtxs_gap_phi","",20,-pi/4,3*pi/4);
  th1d_fvtxn_gap_phi = new TH1D("th1d_fvtxn_gap_phi","",20,-pi/4,3*pi/4);

  th2d_fvtxs_eta = new TH2D("th2d_fvtxs_eta","",20,-3.0,-1.0,20,-3.0,-1.0);
  th2d_fvtxs_phi = new TH2D("th2d_fvtxs_phi","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);
  th2d_fvtxn_eta = new TH2D("th2d_fvtxn_eta","",20,1.0,3.0,20,1.0,3.0);
  th2d_fvtxn_phi = new TH2D("th2d_fvtxn_phi","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);
  th2d_fvtxs_gap_phi = new TH2D("th2d_fvtxs_gap_phi","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);
  th2d_fvtxn_gap_phi = new TH2D("th2d_fvtxn_gap_phi","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);

  th1d_fvtxs_eta_cent0 = new TH1D("th1d_fvtxs_eta_cent0","",20,-3.0,-1.0);
  th1d_fvtxs_phi_cent0 = new TH1D("th1d_fvtxs_phi_cent0","",20,-pi/4,3*pi/4);
  th1d_fvtxn_eta_cent0 = new TH1D("th1d_fvtxn_eta_cent0","",20,1.0,3.0);
  th1d_fvtxn_phi_cent0 = new TH1D("th1d_fvtxn_phi_cent0","",20,-pi/4,3*pi/4);
  th1d_fvtxs_gap_phi_cent0 = new TH1D("th1d_fvtxs_gap_phi_cent0","",20,-pi/4,3*pi/4);
  th1d_fvtxn_gap_phi_cent0 = new TH1D("th1d_fvtxn_gap_phi_cent0","",20,-pi/4,3*pi/4);

  th2d_fvtxs_eta_cent0 = new TH2D("th2d_fvtxs_eta_cent0","",20,-3.0,-1.0,20,-3.0,-1.0);
  th2d_fvtxs_phi_cent0 = new TH2D("th2d_fvtxs_phi_cent0","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);
  th2d_fvtxn_eta_cent0 = new TH2D("th2d_fvtxn_eta_cent0","",20,1.0,3.0,20,1.0,3.0);
  th2d_fvtxn_phi_cent0 = new TH2D("th2d_fvtxn_phi_cent0","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);
  th2d_fvtxs_gap_phi_cent0 = new TH2D("th2d_fvtxs_gap_phi_cent0","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);
  th2d_fvtxn_gap_phi_cent0 = new TH2D("th2d_fvtxn_gap_phi_cent0","",20,-pi/4,3*pi/4,20,-pi/4,3*pi/4);



  // ---------------------------------------------------------------------------------------------------------

  // ---
  // ---

  return EVENT_OK;

}


// --- InitRun, part of Fun4All inheritance
int GreensboroCorrelations::InitRun(PHCompositeNode *topNode)
{

  cout << "InitRun called " << endl;

  int runnumber = 0;

  RunHeader *rh = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if ( !rh )
  {
    cout << PHWHERE << " ERROR::RunHeader not found" << endl;
    return ABORTEVENT;
  }
  runnumber = rh->get_RunNumber();

  // --- set Q-vector offsets
  // cout << "setting Q-vector offsets..." << endl;
  // SetQvectorOffsets(runnumber);
  // SetQvectorOffsetsRBR(runnumber);



  // Setup the utility class
  // This is done in init run so that the collision system can be
  // determined from the run number
  _collsys = "Run16dAu200"; // default to 200 GeV
  use_utils = true;
  // --- Run14AuAu200
  if ( runnumber >= 405839 && runnumber <= 414988 )
    {
      _collsys = "Run14AuAu200";
      use_utils = false;
    }
  // --- Run14HeuAu200
  if ( runnumber >= 415370 && runnumber <= 416893 )
      _collsys = "Run14HeAu200";
  // --- Run15pAu200
  if ( runnumber >= 432637 && runnumber <= 436647 )
    _collsys = "Run15pAu200";
  // --- Run15pAl200
  if ( runnumber >= 436759 && runnumber <= 438422 )
    _collsys = "Run15pAl200";
  // --- Run16dAu200
  if ( runnumber >= 454774 && runnumber <= 455639 )
    _collsys = "Run16dAu200";
  // --- Run16dAu62
  if ( runnumber >= 455792 && runnumber <= 456283 )
    _collsys = "Run16dAu62";
  // --- Run16dAu20
  if ( runnumber >= 456652 && runnumber <= 457298 )
    _collsys = "Run16dAu20";
  // --- Run16dAu39
  if ( runnumber >= 457634 && runnumber <= 458167 )
    _collsys = "Run16dAu39";

  // --- delete this pointer in EndRun
  if ( use_utils )
    {
      cout << "initializing uitls..." << _utils << endl;
      _utils = new dAuBES_utils(_collsys, true);
      cout << "done initializing utils? " << _utils << endl;
      if ( _utils == NULL )
        {
          cout << "WARNING: utils class expected but not found..." << endl;
          use_utils = false;
        }
    }
  // _utils->is_sim(_is_sim);


  return EVENT_OK;
}



int GreensboroCorrelations::EndRun(PHCompositeNode *topNode)
{
  if ( _utils ) delete _utils;
  return EVENT_OK;
}



int GreensboroCorrelations::End(PHCompositeNode *topNode)
{
  if (_verbosity > 1) cout << PHWHERE << "::End() - entered." << endl;
  cout << "total events: " << _ievent << " fraction passing vtx cut: " << tmp_evt * 1.0 / _ievent << endl;
  _output_file->Write();
  _output_file->Close();
  delete _output_file;
  return EVENT_OK;
}
