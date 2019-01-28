const double pi = 3.14159265358979;

TProfile* get_profileR2(int, TH1D*, TH2D*);

void drawme(TProfile*);

void get_dphi()
{

  TFile *file = TFile::Open("test_14434.root");
  TH1D* hcent = (TH1D*)file->Get("th1d_centrality");
  int nevents = hcent->GetEntries();

  TH1D* th1d_fvtxs_phi = (TH1D*)file->Get("th1d_fvtxs_gap_phi");
  TH2D* th2d_fvtxs_phi = (TH2D*)file->Get("th2d_fvtxs_gap_phi");
  TProfile* tp1f_dphiR2 = get_profileR2(nevents,th1d_fvtxs_phi,th2d_fvtxs_phi);

  drawme(tp1f_dphiR2);

}

TProfile* get_profileR2(int nevents, TH1D* th1d_1, TH2D* th2d_12)
{

  TH1D* th1d_2 = th1d_1; // hmm...

  int nbinsx = th2d_12->GetNbinsX();
  int nbinsy = th2d_12->GetNbinsY();
  float xlo = th1d_1->GetBinLowEdge(1);
  float xhi = th1d_1->GetBinLowEdge(nbinsx+1);
  float ylo = th1d_2->GetBinLowEdge(1);
  float yhi = th1d_2->GetBinLowEdge(nbinsy+1);

  TProfile *tp1f_delta_corr = new TProfile(Form("tp1f_delta_corr"), "", nbinsx+nbinsy-1, -pi/4, 3*pi/4, -1e10, 1e10,"");
  TH1D *th1d_delta_acc = new TH1D(Form("th1d_delta_acc"), "", nbinsx+nbinsy-1, -pi/4, 3*pi/4);

  for(int i=0; i<nbinsx; i++)
    {
      for(int j=0; j<nbinsy; j++)
  	{
  	  double content = th2d_12->GetBinContent(i+1,j+1);
  	  double h1 = th1d_1->GetBinContent(i+1);
  	  double h2 = th1d_2->GetBinContent(j+1);
  	  double angle1 = th1d_1->GetBinCenter(i+1);
  	  double angle2 = th1d_2->GetBinCenter(j+1);
  	  double delta = angle2 - angle1;
          if ( delta < -pi/4 ) delta += 2*pi;
          if ( delta > 3*pi/4 ) delta -= 2*pi;
  	  content = content/nevents;
  	  h1 = h1/nevents;
  	  h2 = h2/nevents;
  	  double cumulant = (content/(h1*h2))-1; // correct definition of cumulant
  	  //double cumulant = (content/(h1*h2)); // leave off the -1 to match PHENIX method
          if ( delta == 0 ) continue; // don't understand the big spike at zero...
  	  tp1f_delta_corr->Fill(delta,cumulant);
  	  th1d_delta_acc->Fill(delta);
  	}
    }
  return tp1f_delta_corr;
}



void drawme(TProfile* histo)
{

  TCanvas* c1 = new TCanvas("c1","");

  histo->Draw();

  c1->Print(Form("figs/named_%s.png",histo->GetName()));
  c1->Print(Form("figs/named_%s.pdf",histo->GetName()));

}

