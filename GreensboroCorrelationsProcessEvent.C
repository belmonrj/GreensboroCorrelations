#include <GreensboroCorrelations.h>

#include <iostream>

#include <TMath.h>
//#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TComplex.h>

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



// --- process_event, inherited from Fun4All, does the main part of the analysis on an event by event basis
int GreensboroCorrelations::process_event(PHCompositeNode *topNode)
{

  if (_verbosity > 1) cout << PHWHERE << "::process_event() - entered on event #" << _ievent << endl;
  else if ((_ievent) % 10000 == 0) cout << PHWHERE << "::process_event() - event #" << _ievent << endl;

  // advance event counter
  _ievent++;

  event = _ievent;



  //--------------------------------
  // --- get info from the node tree
  //--------------------------------

  if ( _verbosity > 1 ) cout << "getting the info from the node tree" << endl;

  // --- get the run number
  if ( _verbosity > 1 ) cout << "getting the run information" << endl;
  int runnumber = 0;
  RunHeader *rh = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if (!rh)
  {
    cout << PHWHERE << " ERROR::RunHeader not found" << endl;
    return ABORTEVENT;
  }
  else
  {
    runnumber = rh->get_RunNumber();
  }
  if ( _verbosity > 1 ) cout << "run number is " << runnumber << endl;

  // --- trigger object
  TrigLvl1 *triggers = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  if (!triggers)
  {
    cout << PHWHERE << " ERROR::TrigLvl1 not found" << endl;
    return ABORTEVENT;
  }

  // --- global object (centrality, bbc charge, etc)
  PHGlobal *global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if (!global)
  {
    cout << PHWHERE << " ERROR::PHGlobal not found" << endl;
    return ABORTEVENT;
  }

  // --- event header
  EventHeader *evthead = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!evthead)
  {
    cout << PHWHERE << " ERROR::EventHeader not found" << endl;
    return ABORTEVENT;
  }

  // --- vertex object (bbc event vertex, fvtx event vertex, etc)
  VtxOut *vertexes = findNode::getClass<VtxOut>(topNode, "VtxOut");
  if (!vertexes)
  {
    cout << PHWHERE << " ERROR::VtxOut not found" << endl;
    return ABORTEVENT;
  }

  // --- fvtx track object
  TFvtxCompactTrkMap* trkfvtx_map = findNode::getClass<TFvtxCompactTrkMap>(topNode, "TFvtxCompactTrkMap");
  if (!trkfvtx_map)
  {
    cout << PHWHERE << " No TFvtxCompactTrkMap object !" << endl;
    return ABORTEVENT;
  }

  //---------------------------------------------------------//
  //
  //         Make Event Selection
  //
  //---------------------------------------------------------//

  if ( _verbosity > 1 ) cout << "applying event selection criteria" << endl;

  if ( use_utils )
    {
      if ( _verbosity > 1 ) cout << "using utils to check if event is ok " << endl;
      if (!_utils->is_event_ok(topNode)) return EVENT_OK;
      if ( _verbosity > 1 ) cout << "event passed utils check " << endl;
    }



  //---------------------------------------------------------//
  //
  //         Reading in Global Event Information into Tree
  //
  //---------------------------------------------------------//

  PHPoint precise_vertex1 = vertexes->get_Vertex("SVX_PRECISE");
  vtx_z = precise_vertex1.getZ();
  if (vtx_z != vtx_z) vtx_z = -9999; //NAN check

  PHPoint svx_fast = vertexes->get_Vertex("SVX");//seed vertex
  bc_x = svx_fast.getX();//these are actually the beam center
  bc_y = svx_fast.getY();

  // phglobal fields
  centrality  = global->getCentrality();
  icent = (int)centrality;
  if ( icent < 0 ) icent = 99; // last bin is always empty and this protects against invalid read
  bbc_qn      = global->getBbcChargeN();
  bbc_qs      = global->getBbcChargeS();
  float bbc_charge_sum = bbc_qn+bbc_qs;
  npc1        = global->getNumberPC1Hits();
  event = evthead->get_EvtSequence();
  trigger_scaled = triggers->get_lvl1_trigscaled();
  trigger_live = triggers->get_lvl1_triglive();

  if ( !use_utils && centrality < 0  ) return EVENT_OK;
  if ( !use_utils && centrality > 99 ) return EVENT_OK;



  // --- bbc_z...
  PHPoint vertex1 = vertexes->get_Vertex("BBC");
  bbc_z = vertex1.getZ();
  if ( bbc_z != bbc_z ) bbc_z = -9999; // reassign nan

  if ( !use_utils && fabs(bbc_z) > _cut_zvtx ) return EVENT_OK;

  PHPoint fvtx_vertex = vertexes->get_Vertex("FVTX");
  FVTX_X = fvtx_vertex.getX();
  FVTX_Y = fvtx_vertex.getY();
  FVTX_Z = fvtx_vertex.getZ();
  if ( FVTX_Z != FVTX_Z ) FVTX_Z = -9999; // reassign nan

  // cout << endl;
  // cout << "--- starting vertex checking ---" << endl;
  zvtx = bbc_z;
  if ( use_utils )
    {
      if ( _verbosity > 1 ) cout << "using utils to get vertex " << endl;
      zvtx = _utils->get_vrtx(topNode);
      if ( _verbosity > 1 ) cout << "got the vertex from utils" << endl;
    }

  if ( _verbosity > 1 ) cout << "FVTX vertex points: " << FVTX_X << " " << FVTX_Y << " " << FVTX_Z << endl;




  //int ntr = -1;
  //int ntr = 0;
  nfvtxt = 0;
  nfvtxt_south = 0;
  nfvtxt_north = 0;
  nfvtxt_raw = 0;

  // --- now global variables in header file
  // vector<double> fphi;
  // vector<double> feta;
  // vector<double> fdcax;
  // vector<double> fdcay;
  // vector<double> fchi2ndf;
  // vector<int> fnhitspc;
  // --- reset these for each event
  fphi.clear();
  feta.clear();
  fdcax.clear();
  fdcay.clear();
  fchi2ndf.clear();
  fnhitspc.clear();

  if ( _verbosity > 1 ) cout << "entering fvtx track loop" << endl;
  TFvtxCompactTrkMap::const_iterator trk_iter = trkfvtx_map->range();
  while ( TFvtxCompactTrkMap::const_pointer trk_ptr = trk_iter.next() )
    {
      // --- count the raw number of tracks
      ++nfvtxt_raw;

      // --- get the track object
      TFvtxCompactTrk* fvtx_trk = trk_ptr->get();
      // --- get the track patterns
      bool pattern0 = ((fvtx_trk->get_hit_pattern() & 0x3) > 0);
      bool pattern2 = ((fvtx_trk->get_hit_pattern() & (0x3 << 2)) > 0 );
      bool pattern4 = ((fvtx_trk->get_hit_pattern() & (0x3 << 4)) > 0 );
      bool pattern6 = ((fvtx_trk->get_hit_pattern() & (0x3 << 6)) > 0 );
      int nhits_special = pattern0 + pattern2 + pattern4 + pattern6;

      // --- basic track variables
      float the = fvtx_trk->get_fvtx_theta();
      float eta = fvtx_trk->get_fvtx_eta();
      float phi = fvtx_trk->get_fvtx_phi();
      // int   arm = (int)fvtx_trk->get_arm();
      float fvtx_x      = fvtx_trk->get_fvtx_vtx().getX();
      float fvtx_y      = fvtx_trk->get_fvtx_vtx().getY();
      float fvtx_z      = fvtx_trk->get_fvtx_vtx().getZ();
      float chisq       = fvtx_trk->get_chi2_ndf();
      int   nhits       = (int)fvtx_trk->get_nhits();

      // cout << "----------------------------------------------------" << endl;
      // cout << "pattern " << (int)fvtx_trk->get_hit_pattern() << endl;
      // cout << "part0 " << (fvtx_trk->get_hit_pattern() & 0x3) << endl;
      // cout << "part2 " << (fvtx_trk->get_hit_pattern() & (0x3 << 2)) << endl;
      // cout << "part4 " << (fvtx_trk->get_hit_pattern() & (0x3 << 4)) << endl;
      // cout << "part6 " << (fvtx_trk->get_hit_pattern() & (0x3 << 6)) << endl;
      // cout << "nhits_special " << nhits_special << endl;
      // cout << "nhits " << nhits << endl;
      // cout << "----------------------------------------------------" << endl;

      if ( nhits_special < default_cut_nhit ) continue; // need at least 3 hits in FVTX, excluding VTX

      // fix total momentum to 1.0 (for rotating due to beamtilt)
      double px = 1.0 * sin(the) * cos(phi);
      double py = 1.0 * sin(the) * sin(phi);
      double pz = 1.0 * cos(the);

      if ( use_utils )
	{
	  // rotate based on beamtilt, need to do both rotations with lab frame coordinates
	  double pxprime = _utils->rotate_x(px, pz);
	  double pzprime = _utils->rotate_z(px, pz);
	  // now reassign px and pz to the new rotated frame coordinates
	  px = pxprime;
	  pz = pzprime;
	  phi = TMath::ATan2(py, px);
	  the = TMath::ACos(pz / TMath::Sqrt(px * px + py * py + pz * pz));
	}

      float vertex_z = zvtx;
      if ( FVTX_Z > -999 ) vertex_z = FVTX_Z;
      float DCA_x      = fvtx_x + tan(the) * cos(phi) * (vertex_z - fvtx_z);
      float DCA_y      = fvtx_y + tan(the) * sin(phi) * (vertex_z - fvtx_z);

      th2d_track_before_eta->Fill(centrality,eta);
      th2d_track_before_phi->Fill(centrality,phi);
      th2d_cent_dcax->Fill(centrality,DCA_x);
      th2d_cent_dcay->Fill(centrality,DCA_y);
      th2d_cent_nhitr->Fill(centrality,nhits);
      th2d_cent_nhits->Fill(centrality,nhits_special);
      th2d_cent_chisq->Fill(centrality,chisq);

      // --- if it exists, use the utility class to make the track selections
      if ( use_utils )
        {
          if ( _verbosity > 2 ) cout << "using utils to check if the track is ok " << endl;
          if ( !_utils->is_fvtx_track_ok(fvtx_trk, zvtx) ) continue;
          if ( _verbosity > 2 ) cout << "track pass utils " << endl;
        }
      // --- if it doesn't, make the cuts by hand
      else
        {
          if ( fabs(DCA_x) > default_cut_dca || fabs(DCA_y) > default_cut_dca ) continue;
          if ( nhits < default_cut_nhit ) continue;
          if ( chisq > default_cut_chi2 ) continue;
        }
      // --- done with first loop, so fill after histos, push variables, and count total number of good tracks
      th2d_track_after_eta->Fill(centrality,eta);
      th2d_track_after_phi->Fill(centrality,phi);
      fphi.push_back(phi);
      feta.push_back(eta);
      fdcax.push_back(DCA_x);
      fdcay.push_back(DCA_y);
      fchi2ndf.push_back(chisq);
      fnhitspc.push_back(nhits_special);
      ++nfvtxt;
      if ( eta < 0 ) ++nfvtxt_south;
      if ( eta > 0 ) ++nfvtxt_north;
    } // end first for loop over tracks

  if ( nfvtxt > maxTracks ) return EVENT_OK;

  // ----------------------------------------
  // --- event categorization and special cut
  // ----------------------------------------

  th1d_nfvtxt_combinedER->Fill(nfvtxt);
  th1d_nfvtxt_combined->Fill(nfvtxt);
  th1d_nfvtxt_north->Fill(nfvtxt_north);
  th1d_nfvtxt_south->Fill(nfvtxt_south);
  th2d_nfvtxt_northsouth->Fill(nfvtxt_north,nfvtxt_south);

  th1d_centrality->Fill(centrality);
  th2d_nfvtxt_bbcsum->Fill(nfvtxt,bbc_charge_sum);
  th2d_nfvtxt_centrality->Fill(nfvtxt,centrality);
  th2d_nfvtxt_bbcsumratio->Fill(nfvtxt,bbc_charge_sum/(float)nfvtxt);

  bool passes = PassesTracksChargeRatio(nfvtxt,bbc_charge_sum);
  if ( _collsys == "Run14AuAu200" && !passes )
    {
      if ( _verbosity > 1 ) cout << "Making special event cut for " << _collsys << endl;
      return EVENT_OK; // now testing revised cut...
    }
  th2d_nfvtxt_centralityA->Fill(nfvtxt,centrality);
  th1d_centralityA->Fill(centrality);

  // -------------------------------------------------
  // --- now move on to do more stuff with fvtx tracks
  // -------------------------------------------------

  // --- second fvtx track loop to get the double track cut
  //bool fvtx_track_passes[nfvtxt];
  int number_of_tracks_that_pass = 0;
  if ( do_double_track_cut )
    {
      for ( int i = 0; i < nfvtxt; ++i )
        {
          // --- initialize to true
          fvtx_track_passes[i] = true;
        }
      for ( int i = 0; i < nfvtxt; ++i )
        {
          for ( int j = i+1; j < nfvtxt; ++j )
            {
              double eta1 = feta[i];
              double eta2 = feta[j];
              double phi1 = fphi[i];
              double phi2 = fphi[j];
              th1d_track_deta->Fill(eta1-eta2);
              if ( fabs(eta1-eta2) < 0.0002 )
                {
                  th1d_track_dphi->Fill(phi1-phi2);
                  fvtx_track_passes[i] = false;
                  fvtx_track_passes[j] = false;
                } // check on narrow eta
            } // inner loop
          if ( fvtx_track_passes[i] == true )
            {
              ++number_of_tracks_that_pass;
              th2d_track_aafter_eta->Fill(centrality,feta[i]);
              th2d_track_aafter_phi->Fill(centrality,fphi[i]);
            } // check on pass after nested loop
        } // outer loop
    } // check on do_double_track_cut
  double passratio = (double)number_of_tracks_that_pass/(double)nfvtxt;
  tp1f_track_detacutpass->Fill(centrality,passratio);
  if ( _verbosity > 1 )
    {
      cout << "number of tracks that pass " << number_of_tracks_that_pass
           << " total number of tracks " << nfvtxt
           << " ratio " << passratio << endl;
    }

  // --- do the analysis
  EventStuff();

  if ( _verbosity > 0 ) cout << "sucessfully processed this event, number of fvtx tracks is " << nfvtxt_raw << ", number of fvtx tracks passing cuts is " << nfvtxt << endl;

  ++tmp_evt;//to keep track of how many events pass event cuts

  return EVENT_OK;

} // end process_event





int GreensboroCorrelations::EventStuff()
{

  // ------------------------------------
  // --- initialize the Q-vectors to zero
  // ------------------------------------

  // --- fvtx tracks
  float fvtxs_tracks_qx2[3]; // both, inner, outer
  float fvtxs_tracks_qy2[3];
  float fvtxs_tracks_qx3[3];
  float fvtxs_tracks_qy3[3];
  float fvtxs_tracks_qx4[3];
  float fvtxs_tracks_qy4[3];
  float fvtxs_tracks_qx6[3];
  float fvtxs_tracks_qy6[3];
  float fvtxs_tracks_qw[3];
  float fvtxn_tracks_qx2[3]; // both, inner, outer
  float fvtxn_tracks_qy2[3];
  float fvtxn_tracks_qx3[3];
  float fvtxn_tracks_qy3[3];
  float fvtxn_tracks_qx4[3];
  float fvtxn_tracks_qy4[3];
  float fvtxn_tracks_qx6[3];
  float fvtxn_tracks_qy6[3];
  float fvtxn_tracks_qw[3];

  for ( int i = 0; i < 3; ++i )
  {
    fvtxs_tracks_qx2[i] = 0.0;
    fvtxs_tracks_qy2[i] = 0.0;
    fvtxs_tracks_qx3[i] = 0.0;
    fvtxs_tracks_qy3[i] = 0.0;
    fvtxs_tracks_qx4[i] = 0.0;
    fvtxs_tracks_qy4[i] = 0.0;
    fvtxs_tracks_qx6[i] = 0.0;
    fvtxs_tracks_qy6[i] = 0.0;
    fvtxs_tracks_qw[i] = 0.0;
    fvtxn_tracks_qx2[i] = 0.0;
    fvtxn_tracks_qy2[i] = 0.0;
    fvtxn_tracks_qx3[i] = 0.0;
    fvtxn_tracks_qy3[i] = 0.0;
    fvtxn_tracks_qx4[i] = 0.0;
    fvtxn_tracks_qy4[i] = 0.0;
    fvtxn_tracks_qx6[i] = 0.0;
    fvtxn_tracks_qy6[i] = 0.0;
    fvtxn_tracks_qw[i] = 0.0;
  } // loop over layers

  int ntrack_south_inner = 0;
  int ntrack_north_inner = 0;
  int ntrack_south_outer = 0;
  int ntrack_north_outer = 0;

  // --- initialize Q-vectors for tree
  for ( int i = 0; i < nharm; ++i )
    {
      d_NorthQX[i] = 0;
      d_NorthQY[i] = 0;
      d_SouthQX[i] = 0;
      d_SouthQY[i] = 0;
    }
  d_NorthQW = 0;
  d_SouthQW = 0;

  // --- third fvtxt track loop to calculate Q-vectors
  for ( int i = 0; i < nfvtxt; ++i )
    {
      // --- double track cut
      if ( do_double_track_cut && !fvtx_track_passes[i] ) continue;
      double eta = feta[i];
      double phi = fphi[i];
      double DCA_x = fdcax[i];
      double DCA_y = fdcay[i];
      double chisq = fchi2ndf[i];
      int nhits_special = fnhitspc[i];
      // --- need to do different cuts here
      if ( nhits_special < _cut_nhit ) continue; // need at least 3 hits in FVTX, excluding VTX
      if ( fabs(DCA_x) > _cut_dca || fabs(DCA_y) > _cut_dca ) continue;
      //if ( nhits < _cut_nhit ) continue;
      if ( chisq > _cut_chi2 ) continue;
      // --- Q-vectors for tree
      if ( eta > 0 )
        {
          for ( int i = 0; i < nharm; ++i )
            {
              d_NorthQX[i] += cos(i*phi);
              d_NorthQY[i] += sin(i*phi);
            }
          d_NorthQW += 1;
        }
      if ( eta < 0 )
        {
          for ( int i = 0; i < nharm; ++i )
            {
              d_SouthQX[i] += cos(i*phi);
              d_SouthQY[i] += sin(i*phi);
            }
          d_SouthQW += 1;
        }

      bool is_south = ( eta < 0 );
      bool is_south_inner = ( eta > -2 && eta < 0 );
      bool is_south_outer = ( eta < -2 );
      bool is_north = ( eta > 0 );
      bool is_north_inner = ( eta < 2 && eta > 0 );
      bool is_north_outer = ( eta > 2 );

      if ( is_south )
	{
	  fvtxs_tracks_qx2[0] += cos(2*phi);
	  fvtxs_tracks_qy2[0] += sin(2*phi);
	  fvtxs_tracks_qx3[0] += cos(3*phi);
	  fvtxs_tracks_qy3[0] += sin(3*phi);
	  fvtxs_tracks_qx4[0] += cos(4*phi);
	  fvtxs_tracks_qy4[0] += sin(4*phi);
	  fvtxs_tracks_qx6[0] += cos(6*phi);
	  fvtxs_tracks_qy6[0] += sin(6*phi);
	  fvtxs_tracks_qw[0] += 1;
	}
      if ( is_south_inner )
	{
	  fvtxs_tracks_qx2[1] += cos(2*phi);
	  fvtxs_tracks_qy2[1] += sin(2*phi);
	  fvtxs_tracks_qx3[1] += cos(3*phi);
	  fvtxs_tracks_qy3[1] += sin(3*phi);
	  fvtxs_tracks_qx4[1] += cos(4*phi);
	  fvtxs_tracks_qy4[1] += sin(4*phi);
	  fvtxs_tracks_qx6[1] += cos(6*phi);
	  fvtxs_tracks_qy6[1] += sin(6*phi);
	  fvtxs_tracks_qw[1] += 1;
	  ++ntrack_south_inner;
	}
      if ( is_south_outer )
	{
	  fvtxs_tracks_qx2[2] += cos(2*phi);
	  fvtxs_tracks_qy2[2] += sin(2*phi);
	  fvtxs_tracks_qx3[2] += cos(3*phi);
	  fvtxs_tracks_qy3[2] += sin(3*phi);
	  fvtxs_tracks_qx4[2] += cos(4*phi);
	  fvtxs_tracks_qy4[2] += sin(4*phi);
	  fvtxs_tracks_qx6[2] += cos(6*phi);
	  fvtxs_tracks_qy6[2] += sin(6*phi);
	  fvtxs_tracks_qw[2] += 1;
	  ++ntrack_south_outer;
	}
      if ( is_north )
	{
	  fvtxn_tracks_qx2[0] += cos(2*phi);
	  fvtxn_tracks_qy2[0] += sin(2*phi);
	  fvtxn_tracks_qx3[0] += cos(3*phi);
	  fvtxn_tracks_qy3[0] += sin(3*phi);
	  fvtxn_tracks_qx4[0] += cos(4*phi);
	  fvtxn_tracks_qy4[0] += sin(4*phi);
	  fvtxn_tracks_qx6[0] += cos(6*phi);
	  fvtxn_tracks_qy6[0] += sin(6*phi);
	  fvtxn_tracks_qw[0] += 1;
	}
      if ( is_north_inner )
	{
	  fvtxn_tracks_qx2[1] += cos(2*phi);
	  fvtxn_tracks_qy2[1] += sin(2*phi);
	  fvtxn_tracks_qx3[1] += cos(3*phi);
	  fvtxn_tracks_qy3[1] += sin(3*phi);
	  fvtxn_tracks_qx4[1] += cos(4*phi);
	  fvtxn_tracks_qy4[1] += sin(4*phi);
	  fvtxn_tracks_qx6[1] += cos(6*phi);
	  fvtxn_tracks_qy6[1] += sin(6*phi);
	  fvtxn_tracks_qw[1] += 1;
	  ++ntrack_north_inner;
	}
      if ( is_north_outer )
	{
	  fvtxn_tracks_qx2[2] += cos(2*phi);
	  fvtxn_tracks_qy2[2] += sin(2*phi);
	  fvtxn_tracks_qx3[2] += cos(3*phi);
	  fvtxn_tracks_qy3[2] += sin(3*phi);
	  fvtxn_tracks_qx4[2] += cos(4*phi);
	  fvtxn_tracks_qy4[2] += sin(4*phi);
	  fvtxn_tracks_qx6[2] += cos(6*phi);
	  fvtxn_tracks_qy6[2] += sin(6*phi);
	  fvtxn_tracks_qw[2] += 1;
	  ++ntrack_north_outer;
	}
    } // end third for loop over tracks



  //---------------------------------------------------------//
  //                 finished Get FVTX Tracks
  //---------------------------------------------------------//

  // --------------------------------------------------------------------------------------- //
  // --- calculations and histograms designed to be used with/for acceptance corrections --- //
  // --------------------------------------------------------------------------------------- //

  if ( _verbosity > 1 ) cout << "doing cumulant calculations and filling histograms" << endl;

  // --- FVTX south

  // --- FVTX north

  // --- FVTX north and south combined

  // --- scalar product, fvtxs dot fvtxn

  // --- now have a look at some 4 particle cumulants
  // --- calc4_event has the protection/requirement on the minimum number of tracks
  // --- four particle 2sub





  // --------------------------------------------------------- //
  // --- centrality
  // --------------

  // --- south only
  // --- north only
  // --- combined
  // --- scalar product
  // --- four particle
  // --- four particle 2sub
  // --- six particle

  if ( _verbosity > 2 )
    {
    }

  // ------------------------------------------------------------------------------------- //
  // --- calculations and histograms designed to be used with/for q-vector recentering --- //
  // ------------------------------------------------------------------------------------- //

  // ---
  // --- FVTX south

  // --- FVTX north

  // --- FVTX north and south combined

  // --- scalar product, fvtxs dot fvtxn

  // --- now have a look at some 4 particle cumulants
  // --- calc4_event has the protection/requirement on the minimum number of tracks
  // --- four particle 2sub


  if ( _verbosity > 2 )
    {
      cout << "offset south 2x " << " " << Qoffset_south[2][1].Re()/Qvector_south[0][1].Re() << " " << qvoff_nfvtxt_south[icent][0][2] << endl;
      cout << "offset north 2x " << " " << Qoffset_north[2][1].Re()/Qvector_north[0][1].Re() << " " << qvoff_nfvtxt_north[icent][0][2] << endl;
      cout << "offset south 2y " << " " << Qoffset_south[2][1].Im()/Qvector_south[0][1].Re() << " " << qvoff_nfvtxt_south[icent][1][2] << endl;
      cout << "offset north 2y " << " " << Qoffset_north[2][1].Im()/Qvector_north[0][1].Re() << " " << qvoff_nfvtxt_north[icent][1][2] << endl;
    }

  // --------------------------------------------------------- //
  // --- centrality
  // --------------

  // --- south only
  // --- north only
  // --- combined
  // --- scalar product
  // --- four particle
  // --- four particle 2sub
  // --- six particle

  return EVENT_OK;

} // end of EventStuff

