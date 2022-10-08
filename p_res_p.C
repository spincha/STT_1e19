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
	
	T -> SetBranchAddress("NPar", &NPar);
	T -> SetBranchAddress("TrueP4", &p_true);
	T -> SetBranchAddress("RecoP4", &p_reco);
	T -> SetBranchAddress("Pdg", &pdg);
	T -> SetBranchAddress("VtxX", &x);
	T -> SetBranchAddress("VtxY", &y);
	T -> SetBranchAddress("VtxZ", &z);
	T -> SetBranchAddress("ParentId", &p_id);
	
	
	int nentry = T -> GetEntries();
	
	TH1F *p_res[20]; double res[20], mom[20];
	
	TH2F *hist_res = new TH2F("Res storage", "Binwise Resolution", 100, 0, 10, 100, -20, 20);
	
	TH1F *p_sim = new TH1F("Simp_mu", "Simulated Muon 3-Momentum", 100, 0, 10);
	TH1F *p_rec = new TH1F("Rec_mu", "Reconstructed Muon 3-Momentum", 100, 0, 10);
	
//	TGraphErrors *p_res_p = new TGraphErrors();

	double k = 0;
	
	for(int num = 0; num < 100; num++)
	{	
//		p_res[num] = new TH1F(Form("p_res_%d", num), "Momentum resolution", 100, -0.5, 0.5);
		
		for(int i = 0; i < nentry; i++)
		{
			T -> GetEntry(i);
			
			if(!inFV(x, y, z)) continue;
			for(int j = 0; j < NPar; j++)
			{
				if(abs(pdg[j]) == 13 && p_id[j] == -1)
				{
					pvec_true = (sqrt(pow(p_true[j][0],2) + pow(p_true[j][1],2) + pow(p_true[j][2],2)))/1000;
					pvec_reco = (sqrt(pow(p_reco[j][0],2) + pow(p_reco[j][1],2) + pow(p_reco[j][2],2)))/1000;
					
//					cout << "Momentum: " << pvec_true << endl;
					
					res_func = (pvec_true - pvec_reco)/pvec_true;
					
					if(pvec_true >= num*0.1 && pvec_true < (num+1)*0.1)		//0.1 factor is to check in GeV units
					{
						hist_res -> Fill(num*0.1, res_func * 100);
//						p_res_p -> SetPoint(num, num*0.5, abs((pvec_true - pvec_reco)/pvec_true)*100);

					}
				}
				
			}
		
		}

//		p_res_p -> SetPoint(num, num*0.5, (p_res[num] -> GetStdDev())*100);
//		p_res_p -> SetPointError(num, 0.5, p_res[num] -> GetMean());
	}
	TGraph *std_dev_p = new TGraph(20, mom, res);
	
	c1 -> cd();
	hist_res -> GetXaxis() -> SetTitle("True Momentum (GeV/c)");
	hist_res -> GetYaxis() -> SetTitle("(p_{true} - p_{reco})/p_{true}");
	hist_res -> GetZaxis() -> SetTitle("Entries");
	hist_res -> Draw("colz");
	
//	p_res_p -> Fit("pol0");
//	p_res_p -> SetStats(0);

	c2 -> cd();
	hist_res -> Draw("surf");

}
