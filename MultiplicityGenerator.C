enum eData {
        kTOFRPC_Data                          =  0,
        kTOF_Data                                 ,
        kRPC_Data                                 ,
        kTracksMDC_Data                           ,
        kNumData              
    };

TString DataHistoName[kNumData]=
    {
        "hitsTOF+RPC_selected"                         ,
        "hitsTOF_selected"                             ,
        "hitsRPC_selected"                             ,
        "tracksMDC_selected"                       
    };

enum eNewCentralityEstimator {
        kTOFRPC                          =  0,
        kTOF                                 ,
        kRPC                                 ,
        kTOFRPCtot                           ,
        kTOFtot                              ,
        kRPCtot                              ,
        kSelectedParticleCand                ,
        kSelectedParticleCandCorr            ,
        kSelectedParticleCandNorm            ,
        kSelectedParticleCandSecNorm         ,
//        kSelectedParticleCandCorrPerWire     ,
        kPrimaryParticleCand                 ,
        kMdcWires                            ,
//        kMdcWiresOuterMod                    ,
        kFWSumChargeSpec                     ,
        kFWSumChargeZ                        ,
        kDirectivity                         ,
        kRatioEtEz                           ,
        kEt                                  ,
        kNumCentralityEstimator              
    };

TString EstimatorHistoName[kNumCentralityEstimator]=
    {
        "TOFRPC"                          ,
        "TOF"                             ,
        "RPC"                             ,
        "TOFRPCtot"                       ,
        "TOFtot"                          ,
        "RPCtot"                          ,
        "SelectedParticleCand"            ,
        "SelectedParticleCandCorr"        ,
        "SelectedParticleCandNorm"        ,
        "SelectedParticleCandSecNorm"     ,
//        "SelectedParticleCandCorrPerWire" ,
        "PrimaryParticleCand"             ,
        "MdcWires"                        ,
//        "MdcWiresOuterMod"                ,
        "FWSumChargeSpec"                 ,
        "FWSumChargeZ"                    ,
        "Directivity"                     ,
        "RatioEtEz"                       ,
        "Et"
    };

void MultiplicityGenerator()
{
	TFile *f  = new TFile("/home/vad/NIR_codes/Centrality/gmc-Au3Au3-snn23.6-md0.4-nd-1.0-rc1-smax99.0.root");
//	TFile *f  = new TFile("/home/vad/NIR_codes/Centrality/glau_ntuple.root");
	TFile *f2 = new TFile("/home/vad/NIR_codes/QAHistoBuilder/QASelectedPT2Histos.root");
	TFile *f3 = new TFile("/home/vad/NIR_codes/Centrality/centrality_epcorr_apr12_gen8_2018_07.root");
	TFile *f4 = new TFile("/home/vad/bkardan/GlauberMC/root/GlauberMC_Au3Au3_ntuple_23.6mb_r6.5541_a0.523_m0.9_withEstimator.root");
	std::unique_ptr<TTree> t {(TTree*)f->Get("nt_Au3_Au3")};
	TTree *t1=(TTree*)f->Get("nt_Au3_Au3");
	TTree *t2=(TTree*)f4->Get("nt_Au3_Au3");

	Float_t Npart, Ncoll, Nhard, Npart0, B;

	t->SetBranchAddress("Npart", &Npart);
	t->SetBranchAddress("Ncoll", &Ncoll);
	t->SetBranchAddress("Nhard", &Nhard);
	t->SetBranchAddress("Npart0", &Npart0);
	t->SetBranchAddress("B", &B);

	TH1F *hNpart       = new TH1F("Npart", ";Npart;counts", 400, 0, 400);
	TH1F *hNcoll       = new TH1F("Ncoll", ";Ncoll;counts", 750, 0, 750);
	TH1F *hNhard       = new TH1F("Nhard", ";Nhard;counts", 750, 0, 750);
	TH1F *hB           = new TH1F("B", ";B, fm;counts", 200, 0, 20);
	TH1F *hNpart0      = new TH1F("Npart0", ";Npart0;counts", 400, 0, 400);
	TH1F *hNhitstofrpc = new TH1F("Nhitstofrpc", ";Nhitstofrpc;counts", 300, 0, 300);
	TH1F hNtracks;
	TH1F hNhitstof;
	TH1F hNhitsrpc;
	
	TH2F *hNpart_vs_Ncoll  = new TH2F("Npart_vs_Ncoll", ";Npart;Ncoll", 400, 0, 400, 750, 0, 750);
	TH2F *hNpart_vs_Nhard  = new TH2F("Npart_vs_Nhard", ";Npart;Nhard", 400, 0, 400, 750, 0, 750);
	TH2F *hNpart_vs_Npart0 = new TH2F("Npart_vs_Npart0", ";Npart;Npart0", 400, 0, 400, 400, 0, 400);
	TH2F *hB_vs_Npart      = new TH2F("B_vs_Npart", ";B, fm;Npart", 200, 0, 20, 400, 0, 400);
	TH2F *hB_vs_Ncoll      = new TH2F("B_vs_Ncoll", ";B, fm;Ncoll", 200, 0, 20, 750, 0, 750);

	TH1F *CompareHistos[12];
	for (int i = 0; i < 3; i++) CompareHistos[i] = (TH1F*)f3->Get(Form("/EstimatorHist/%s", EstimatorHistoName[i].Data()));	
	CompareHistos[3] = (TH1F*)f3->Get(Form("/EstimatorHist/%s", EstimatorHistoName[6].Data()));
	for (int i = 4; i < 8; i++)  CompareHistos[i] = (TH1F*)f2->Get(Form("%s", DataHistoName[i-4].Data()));
	for (int i = 8; i < 11; i++) CompareHistos[i] = new TH1F(EstimatorHistoName[i-8],";centralityEstimator;counts",500,0.0,500.0);
	CompareHistos[11] = new TH1F(EstimatorHistoName[6],";centralityEstimator;counts",500,0.0,500.0);
	for (int i = 0; i < 4; i++)  CompareHistos[i] -> SetLineColor(1);
	for (int i = 4; i < 8; i++)  CompareHistos[i] -> SetLineColor(2);
	for (int i = 8; i < 12; i++) CompareHistos[i] -> SetLineColor(4);
	for (int i = 0; i < 12; i++) CompareHistos[i] -> SetLineWidth(4);

	Float_t fTOFRPC, fTOF, fRPC, fSelectedParticleCand, fPrimaryParticleCand;

	t2->SetBranchAddress("TOFRPC", &fTOFRPC);
	t2->SetBranchAddress("TOF", &fTOF);
	t2->SetBranchAddress("RPC", &fRPC);
	t2->SetBranchAddress("SelectedParticleCand", &fSelectedParticleCand);
	t2->SetBranchAddress("PrimaryParticleCand", &fPrimaryParticleCand);

	Int_t entries = (Int_t)t2->GetEntries();

	for (int i = 0; i < entries; i++)
	{
		t2->GetEntry(i);
		CompareHistos[8]  -> Fill(fTOFRPC);
		CompareHistos[9]  -> Fill(fTOF);
		CompareHistos[10] -> Fill(fRPC);
		CompareHistos[11] -> Fill(fSelectedParticleCand);
//		CompareHistos[11] -> Fill(fPrimaryParticleCand);
	}

	entries = (Int_t)t->GetEntries();

	for (int i = 0; i < entries; i++)
	{
		t->GetEntry(i);

		hNpart  -> Fill(Npart);
		hNcoll  -> Fill(Ncoll);
		hNhard  -> Fill(Nhard);
		hNpart0 -> Fill(Npart0);
		hB      -> Fill(B);

		hNpart_vs_Ncoll  -> Fill(Npart, Ncoll);
		hNpart_vs_Nhard  -> Fill(Npart, Nhard);
		hNpart_vs_Npart0 -> Fill(Npart, Npart0);
		hB_vs_Npart      -> Fill(B, Npart);
		hB_vs_Ncoll      -> Fill(B, Ncoll);
	}
	
	Glauber::Fitter *fitter=new Glauber::Fitter(move(t));

	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f2->Get("tracksMDC_selected")));
	fitter -> Glauber::Fitter::Init(entries);
	fitter -> Glauber::Fitter::SetGlauberFitHisto (1.0, 0.24, 20.34, entries, 1.10e-7, false);
	hNtracks  = fitter->Glauber::Fitter::GetGlauberFitHisto();

	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f2->Get("hitsTOF_selected")));
	fitter -> Glauber::Fitter::Init(entries);
	fitter -> Glauber::Fitter::SetGlauberFitHisto (1.0, 0.20, 6.36, entries, 1.64e-6, false);
	hNhitstof = fitter->Glauber::Fitter::GetGlauberFitHisto();

	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f2->Get("hitsRPC_selected")));
	fitter -> Glauber::Fitter::Init(entries);
	fitter -> Glauber::Fitter::SetGlauberFitHisto (1.0, 0.50, 29.06, entries, 1.64e-6, false);
	hNhitsrpc = fitter->Glauber::Fitter::GetGlauberFitHisto();

	float alpha=1.64e-6;
	for (int i = 0; i < entries; i++)
	{
		cout<<i<<endl;
		t1->GetEntry(i);
		const int Na = int(fitter->Glauber::Fitter::Nancestors(1.0));
		float nHits {0.};
		
		fitter->Glauber::Fitter::SetNBDhist(0.20, 6.36);
                std::unique_ptr<TH1F> htemp1 {(TH1F*)((fitter->Glauber::Fitter::GetNBDHisto()).Clone("htemp"))};
        	float nHitstof {0.};
        	for (int j=0; j<Na; j++){
//           	 	nHitstof += ((1-alpha*Npart*Npart)*int(htemp1->GetRandom()));
			if (gRandom->Rndm() < (1-alpha*Npart*Npart)) nHitstof += (int(htemp1->GetRandom()));
        		}
        	
		fitter->Glauber::Fitter::SetNBDhist(0.50, 29.06);
                std::unique_ptr<TH1F> htemp2 {(TH1F*)((fitter->Glauber::Fitter::GetNBDHisto()).Clone("htemp"))};
        	float nHitsrpc {0.};
        	for (int j=0; j<Na; j++){
//            		nHitsrpc += ((1-alpha*Npart*Npart)*int(htemp2->GetRandom()));
			if (gRandom->Rndm() < (1-alpha*Npart*Npart)) nHitsrpc += (int(htemp2->GetRandom()));
        		}
		
		nHits=nHitstof+nHitsrpc;
		hNhitstofrpc -> Fill(nHits);
	}

	hNhitstofrpc -> SetLineColor(3);
	hNtracks.SetLineColor(3);
	hNhitstof.SetLineColor(3);
	hNhitsrpc.SetLineColor(3);
	hNhitstofrpc -> SetLineWidth(4);
	hNtracks.SetLineWidth(4);
	hNhitstof.SetLineWidth(4);
	hNhitsrpc.SetLineWidth(4);

	TFile *f1 = new TFile("Generator.root", "recreate");

	Int_t bins;
	for (int i = 0; i < 4; i++) {
		bins = CompareHistos[i] -> GetNbinsX();
		CompareHistos[i] -> Scale(1/(CompareHistos[i] -> Integral(4,bins)));
		if (i == 0) {
			hNhitstofrpc -> Scale((CompareHistos[i] -> Integral(100,150))/(hNhitstofrpc -> Integral(100,150)));
			CompareHistos[i+4] -> Scale((CompareHistos[i] -> Integral(100,150))/(CompareHistos[i+4] -> Integral(100,150)));
			CompareHistos[i+8] -> Scale((CompareHistos[i] -> Integral(100,150))/(CompareHistos[i+8] -> Integral(100,150)));
			}
		if (i == 1) {
			hNhitstof.Scale((CompareHistos[i] -> Integral(20,40))/(hNhitstof.Integral(20,40)));
			CompareHistos[i+4] -> Scale((CompareHistos[i] -> Integral(20,40))/(CompareHistos[i+4] -> Integral(20,40)));
			CompareHistos[i+8] -> Scale((CompareHistos[i] -> Integral(20,40))/(CompareHistos[i+8] -> Integral(20,40)));
			}
		if (i == 2) {
			hNhitsrpc.Scale((CompareHistos[i] -> Integral(50,100))/(hNhitsrpc.Integral(50,100)));
			CompareHistos[i+4] -> Scale((CompareHistos[i] -> Integral(50,100))/(CompareHistos[i+4] -> Integral(50,100)));
			CompareHistos[i+8] -> Scale((CompareHistos[i] -> Integral(50,100))/(CompareHistos[i+8] -> Integral(50,100)));
                        }
		if (i == 3) {
			hNtracks.Scale((CompareHistos[i] -> Integral(20,50))/(hNtracks.Integral(20,50)));	
			CompareHistos[i+4] -> Scale((CompareHistos[i] -> Integral(20,50))/(CompareHistos[i+4] -> Integral(20,50)));
			CompareHistos[i+8] -> Scale((CompareHistos[i] -> Integral(20,50))/(CompareHistos[i+8] -> Integral(20,50)));
			}
		
		}
	
	TCanvas* c1 = new TCanvas("TOFRPC","TOFRPC");
	TLegend legend1(0.1, 0.2, 0.3, 0.4, "TOFRPC");
	legend1.AddEntry(CompareHistos[0], "Behruz_TOFRPC");
	legend1.AddEntry(CompareHistos[4], "Data_TOFRPC");
	legend1.AddEntry(CompareHistos[8], "Behruz'sCode_TOFRPC");
	legend1.AddEntry(hNhitstofrpc, "CentralityFramework_TOFRPC");
	THStack *TOFRPC = new THStack("TOFRPC","TOFRPC");
	TOFRPC -> Add(CompareHistos[0]);
	TOFRPC -> Add(CompareHistos[4]);
	TOFRPC -> Add(CompareHistos[8]);
	TOFRPC -> Add(hNhitstofrpc);
	TOFRPC -> Draw("nostack");
	legend1.Draw("same");
	c1->Write();
	
	TCanvas* c2 = new TCanvas("TOF","TOF");
	TLegend legend2(0.1, 0.2, 0.3, 0.4, "TOF");
	legend2.AddEntry(CompareHistos[1], "Behruz_TOF");
	legend2.AddEntry(CompareHistos[5], "Data_TOF");
	legend2.AddEntry(CompareHistos[9], "Behruz'sCode_TOF");
	legend2.AddEntry(&hNhitstof, "CentralityFramework_TOF");
	THStack *TOF = new THStack("TOF","TOF");
	TOF -> Add(CompareHistos[1]);
	TOF -> Add(CompareHistos[5]);
	TOF -> Add(CompareHistos[9]);
	TOF -> Add(&hNhitstof);
	TOF -> Draw("nostack");
	legend2.Draw("same");
	c2->Write(); 

	TCanvas* c3 = new TCanvas("RPC","RPC");
	TLegend legend3(0.1, 0.2, 0.3, 0.4, "RPC");
	legend3.AddEntry(CompareHistos[2], "Behruz_RPC");
	legend3.AddEntry(CompareHistos[6], "Data_RPC");
	legend3.AddEntry(CompareHistos[10], "Behruz'sCode_RPC");
	legend3.AddEntry(&hNhitsrpc, "CentralityFramework_RPC");
	THStack *RPC = new THStack("RPC","RPC");
	RPC -> Add(CompareHistos[2]);
	RPC -> Add(CompareHistos[6]);
	RPC -> Add(CompareHistos[10]);
	RPC -> Add(&hNhitsrpc);
	RPC -> Draw("nostack");
	legend3.Draw("same");
	c3->Write();
	
	TCanvas* c4 = new TCanvas("tracksMDC","tracksMDC");
	TLegend legend4(0.1, 0.2, 0.3, 0.4, "tracksMDC");
	legend4.AddEntry(CompareHistos[3], "Behruz_tracksMDC");
	legend4.AddEntry(CompareHistos[7], "Data_tracksMDC");
	legend4.AddEntry(CompareHistos[11], "Behruz'sCode_tracksMDC");
	legend4.AddEntry(&hNtracks, "CentralityFramework_tracksMDC");
	THStack *tracksMDC = new THStack("tracksMDC","tracksMDC");
	tracksMDC -> Add(CompareHistos[3]);
	tracksMDC -> Add(CompareHistos[7]);
	tracksMDC -> Add(CompareHistos[11]);
	tracksMDC -> Add(&hNtracks);
	tracksMDC -> Draw("nostack");
	legend4.Draw("same");
	c4->Write();

	for (int j = 0; j < 12; j++) CompareHistos[j] -> Write();

	hNpart    -> Write();
	hNcoll    -> Write();
	hNhard    -> Write();
	hNpart0   -> Write();
	hB        -> Write();

	hNtracks.Write();
	hNhitstof.Write();
	hNhitsrpc.Write();
	hNhitstofrpc -> Write();

	hNpart_vs_Ncoll  -> Write();
	hNpart_vs_Nhard  -> Write();
	hNpart_vs_Npart0 -> Write();
	hB_vs_Npart      -> Write();
	hB_vs_Ncoll      -> Write();

	f1->Close();
	f2->Close();
	f3->Close();
	f ->Close();
	

}
