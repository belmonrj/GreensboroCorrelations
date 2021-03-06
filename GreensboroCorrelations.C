#include <GreensboroCorrelations.h>

#include <iostream>
#include <fstream>



#include <TOAD.h>
#include <TMath.h>
#include <TFile.h>
//#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TComplex.h>

// Fun4All includes...
#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <PHGlobal.h>
#include <EventHeader.h>
#include <TrigLvl1.h>
#include <VtxOut.h>
#include <PHPoint.h>
#include <RunHeader.h>
#include <TFvtxCompactTrkMap.h>
#include <TFvtxCompactCoordMap.h>
#include <dAuBES_utils.h>




using namespace std;


// --- class constructor
GreensboroCorrelations::GreensboroCorrelations(): SubsysReco("GREENSBOROCORRELATIONS")
{
  //ResetEvent(NULL);

  _ievent = 0;
  _verbosity = 0;
  _output_filename = "NULL";
  _output_file = NULL;
  do_double_track_cut = false;
  _cut_zvtx = 10.0;
  _cut_chi2 = 5.0;
  _cut_nhit = 3;
  _cut_dca = 2.0;
  _utils = NULL;
  _collsys = "NULL";

  use_utils = true;

  tmp_evt = 0;


  FVTX_X = -9999.9;
  FVTX_Y = -9999.9;
  FVTX_Z = -9999.9;

  event = -9999;
  bbc_z = -9999.9;
  centrality = -9999;
  icent = 0;
  npc1 = -9999;
  trigger_scaled = -9999;
  trigger_live = -9999;
  bc_x = -9999.9;
  bc_y = -9999.9;
  vtx_z = -9999.9;
  zvtx = -9999.9;

  bbc_qn = -9999.9;
  bbc_qs = -9999.9;
  nfvtxt = 0;
  nfvtxt_south = 0;
  nfvtxt_north = 0;
  nfvtxt_raw = 0;
  for ( int i = 0; i < maxTracks; ++i ) fvtx_track_passes[i] = false;
  for ( int i = 0; i < nharm; ++i )
    {
      d_SouthQX[i] = 0;
      d_SouthQY[i] = 0;
      d_NorthQX[i] = 0;
      d_NorthQY[i] = 0;
    }
  d_SouthQW = 0;
  d_NorthQW = 0;


  th1d_cnt_both_phi = NULL;
  th1d_cnt_east_phi = NULL;
  th1d_cnt_west_phi = NULL;
  th1d_cnt_both_phi_high = NULL;
  th1d_cnt_east_phi_high = NULL;
  th1d_cnt_west_phi_high = NULL;

  th1d_nfvtxt_combinedER = NULL;
  th1d_nfvtxt_combined = NULL;
  th1d_nfvtxt_north = NULL;
  th1d_nfvtxt_south = NULL;
  th2d_nfvtxt_northsouth = NULL;
  th1d_centrality = NULL;
  th1d_centralityA = NULL;
  th2d_nfvtxt_bbcsum = NULL;
  th2d_nfvtxt_centrality = NULL;
  th2d_nfvtxt_centralityA = NULL;
  th2d_nfvtxt_bbcsumratio = NULL;
  th1d_track_deta = NULL;
  th1d_track_dphi = NULL;
  th2d_track_before_eta = NULL;
  th2d_track_before_phi = NULL;
  th2d_track_after_eta = NULL;
  th2d_track_after_phi = NULL;
  th2d_track_aafter_eta = NULL;
  th2d_track_aafter_phi = NULL;
  th2d_cent_dcax = NULL;
  th2d_cent_dcay = NULL;
  th2d_cent_nhitr = NULL;
  th2d_cent_nhits = NULL;
  th2d_cent_chisq = NULL;
  tp1f_track_detacutpass = NULL;

  th1d_fvtxs_eta = NULL;
  th1d_fvtxs_phi = NULL;
  th1d_fvtxn_eta = NULL;
  th1d_fvtxn_phi = NULL;
  th1d_fvtxs_gap_phi = NULL;
  th1d_fvtxn_gap_phi = NULL;

  th2d_fvtxs_eta = NULL;
  th2d_fvtxs_phi = NULL;
  th2d_fvtxn_eta = NULL;
  th2d_fvtxn_phi = NULL;
  th2d_fvtxs_gap_phi = NULL;
  th2d_fvtxn_gap_phi = NULL;

  th1d_fvtxs_eta_cent0 = NULL;
  th1d_fvtxs_phi_cent0 = NULL;
  th1d_fvtxn_eta_cent0 = NULL;
  th1d_fvtxn_phi_cent0 = NULL;
  th1d_fvtxs_gap_phi_cent0 = NULL;
  th1d_fvtxn_gap_phi_cent0 = NULL;

  th2d_fvtxs_eta_cent0 = NULL;
  th2d_fvtxs_phi_cent0 = NULL;
  th2d_fvtxn_eta_cent0 = NULL;
  th2d_fvtxn_phi_cent0 = NULL;
  th2d_fvtxs_gap_phi_cent0 = NULL;
  th2d_fvtxn_gap_phi_cent0 = NULL;


  for ( int i = 0; i < maxTracks; ++i )
    {
      for ( int cs = 0; cs < 2; ++cs )
        {
          for ( int c = 0; c < maxHarmonic; ++c )
            {
              qvoff_nfvtxt[i][cs][c] = 0;
              qvoff_nfvtxt_north[i][cs][c] = 0;
              qvoff_nfvtxt_south[i][cs][c] = 0;
            }
        }
    }

  for ( int i = 0; i < 100; ++i )
    {
      for ( int cs = 0; cs < 2; ++cs )
        {
          for ( int c = 0; c < maxHarmonic; ++c )
            {
              qvoff_cent[i][cs][c] = 0;
              qvoff_cent_north[i][cs][c] = 0;
              qvoff_cent_south[i][cs][c] = 0;
            }
        }
    }

} // end of class constructor


// --- class destructor
GreensboroCorrelations::~GreensboroCorrelations()
{
  //if ( _utils ) delete _utils;
}
