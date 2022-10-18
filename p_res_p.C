const double centerZ=23910;
const double centerY=-2384.73;

TGeoManager *geo = NULL;
bool inSTT(double X, double Y, double Z){

 TString name = geo->FindNode(X, Y, Z)->GetName();
 if (name.Contains("STT")||name.Contains("horizontalST")) return true;
 return false;
}

bool inSTT(TVector3 pos){

 TString name = geo->FindNode(pos.X(), pos.Y(), pos.Z())->GetName();
 if (name.Contains("STT")||name.Contains("horizontalST")) return true;
 return false;
}

bool inFV(double X, double Y, double Z){
 
  if (!inSTT(X, Y, Z)) return false;
  if (sqrt(pow(Y-centerY,2) + pow(Z-centerZ,2)) > 1800) return false;

  return true;

}

void p_res_p(){
	
	TCanvas *c1 = new TCanvas();
	TCanvas *c2 = new TCanvas();

	TFile *f = new TFile("reco2308.root", "READ");
	TTree *T = (TTree*)f -> Get("edep_smeared_tree");
	
	TFile *geofile = new TFile("/home/shailesh/EDEP/stt_aug/edep_files/test_v2308_0.edep.root");
	geo = (TGeoManager *)geofile->Get("EDepSimGeometry");
	
	int NPar, arbNPar = 300;
	int pdg[arbNPar], p_id[arbNPar];
	double p_true[arbNPar][4], p_reco[arbNPar][4], pvec_true, pvec_reco, res_func;
	double x, y, z;
	double p_min = 0, p_max = 10.0;
	
	T -> SetBranchAddress("NPar", &NPar);
	T -> SetBranchAddress("TrueP4", &p_true);
	T -> SetBranchAddress("RecoP4", &p_reco);
	T -> SetBranchAddress("Pdg", &pdg);
	T -> SetBranchAddress("VtxX", &x);
	T -> SetBranchAddress("VtxY", &y);
	T -> SetBranchAddress("VtxZ", &z);
	T -> SetBranchAddress("ParentId", &p_id);
	
	
	int nentry = T -> GetEntries();
	cout << nentry << endl;
	int steps = 25.0;
	TH1F *p_res[steps];
//	TF1* gausfit[steps];
	
/*	TH1F *p_reso = new TH1F("p_res", "Momentum Resolution of #mu^{+}/#mu^{-} in STT", 100, -0.2, 0.2);
	TH1F *p_sim = new TH1F("Simp_mu", "Simulated Muon 3-Momentum", 100, p_min, p_max);
	
	for(int i = 0; i < nentry; i++)
	{
		T -> GetEntry(i);
		
		if(!inFV(x, y, z)) continue;
		for(int j = 0; j < NPar; j++)
		{
			if(abs(pdg[j]) == 13 && p_id[j] == -1)
			{
				pvec_true = (sqrt(pow(p_true[j][0],2) + pow(p_true[j][1],2) + pow(p_true[j][2],2)))/1000.0;
				pvec_reco = (sqrt(pow(p_reco[j][0],2) + pow(p_reco[j][1],2) + pow(p_reco[j][2],2)))/1000.0;
				
				if(pvec_true > p_min && pvec_true < p_max)
				{
					res_func = (pvec_true - pvec_reco)/pvec_true;
					
					p_reso -> Fill(res_func);
					p_sim -> Fill(pvec_true);
				}
			}
			
		}
	
	}
	TH2F *hist_res = new TH2F("Res storage", "Binwise Resolution", steps, p_min, p_max, 100, -20, 20);
*/	
	TGraphErrors *p_res_p = new TGraphErrors();
	

	for(int num = 0; num < steps; num++)
	{
		p_res[num] = new TH1F(Form("p_res_%d", num), "Momentum resolution", 100, -0.2, 0.2);
		TF1* gausfit = new TF1("ranged gaussian", "gaus", -0.075, 0.075);
		
		for(int i = 0; i < nentry; i++)
		{
			T -> GetEntry(i);
			
			if(!inFV(x, y, z)) continue;
			for(int j = 0; j < NPar; j++)
			{
				if(abs(pdg[j]) == 13 && p_id[j] == -1)
				{
					pvec_true = (sqrt(pow(p_true[j][0],2) + pow(p_true[j][1],2) + pow(p_true[j][2],2)))/1000.0;
					pvec_reco = (sqrt(pow(p_reco[j][0],2) + pow(p_reco[j][1],2) + pow(p_reco[j][2],2)))/1000.0;
					
					res_func = (pvec_true - pvec_reco)/pvec_true;
					
					if(pvec_true >= num*(p_max - p_min)/steps && pvec_true < (num+1)*(p_max - p_min)/steps)
					{
//						hist_res -> Fill((num+0.5)*(p_max - p_min)/steps, res_func * 100);
						p_res[num] -> Fill(res_func);
					}
				}
				
			}
		
		}
		p_res[num] -> Fit("ranged gaussian", "R");
		
		p_res_p -> SetPoint(num, (num + 0.5)*(p_max - p_min)/steps, (gausfit -> GetParameter(2)) * 100);
		p_res_p -> SetPointError(num, 0.5*(p_max - p_min)/steps, (gausfit -> GetParError(2)) * 100);
	}
	
	c1 -> cd();
/*	p_sim -> Draw();
	p_sim -> GetXaxis() -> SetTitle("Momentum (in GeV)");
	p_sim -> GetYaxis() -> SetTitle("Entries");
	hist_res -> GetXaxis() -> SetTitle("True Momentum (GeV/c)");
	hist_res -> GetYaxis() -> SetTitle("(p_{true} - p_{reco})/p_{true} (in %)");
	hist_res -> GetZaxis() -> SetTitle("Entries");
	hist_res -> Draw("colz");
*/	
	c2 -> cd();
	TF1* linfit = new TF1("ranged linear", "pol0", 0.5, 10.0);
	p_res_p -> Draw("AP");
	p_res_p -> SetTitle("Momentum Resolution vs Momentum");
	p_res_p -> GetXaxis() -> SetTitle("True Momentum (GeV/c)");
	p_res_p -> GetYaxis() -> SetTitle("Sigma of (p_{true} - p_{reco})/p_{true} (in %)");
	p_res_p -> SetMinimum(0);
	p_res_p -> SetMaximum(10);
	p_res_p -> Fit("pol0");
	gStyle -> SetOptFit(0011);
//	p_res_p -> SetStats(0);
	
	TLegend *leg = new TLegend();
	leg -> AddEntry((TObject*)0, "Gaussian fit for each histogram between -0.075 to 0.075", "C");
	gStyle -> SetLegendFont(200);
	leg -> Draw();
		
/*	c2 -> cd();
	TF1* gaussfit = new TF1("ranged gaussian", "gaus", -0.075, 0.075);
	p_reso -> Draw();
	p_reso -> GetXaxis() -> SetTitle("(p_{true} - p_{reco})/p_{true}");
	p_reso -> GetYaxis() -> SetTitle("Entries");
	p_reso -> Fit("ranged gaussian", "R");
	gStyle -> SetOptFit(0011);
*/
}
