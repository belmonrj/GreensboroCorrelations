void Run_GreensboroCorrelationsRun14(const char *outFile = "test_train_output.root")
{

  //-----------------//
  //--- Libraries ---//
  //-----------------//

  gSystem->Load("libdAuBES_utils.so");
  gSystem->Load("libGreensboroCorrelations.so");
  gSystem->ListLibraries();

  //----------------------//
  //--- Fun4All Server ---//
  //----------------------//

  Fun4AllServer *se = Fun4AllServer::instance();

  // To get the FVTX tracks
  se->registerSubsystem(new FvtxReadbackDST());

  //--------------------//
  //--- User Module ----//
  //--------------------//

  GreensboroCorrelations *sflow = new GreensboroCorrelations();
  sflow->set_output_filename(outFile);
  sflow->set_do_double_track_cut(true);
  sflow->set_zvtxcut(10.0); // z-vertex cut in cm (default is 10)
  sflow->set_chi2cut(5.0);  // chi2/ndf cut on tracks (default is 5)
  sflow->set_nhitcut(3);    // number of hits in tracks (default is 3)
  sflow->set_dcacut(2.0);   // dca cut on tracks in cm (default is 2)
  sflow->Verbosity(0);
  se->registerSubsystem(sflow);

}

void InputData(vector<string> &indata)
{
  //indata.push_back("CNT"); // must have for Run14 (HeAu only??)
  indata.push_back("MWG"); // must have for Run14 (only)
  indata.push_back("MuonDST"); // must have for Run14 (AuAu only)
  return;
}
