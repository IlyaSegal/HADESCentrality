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

enum eCentralityClass { 
        k0010TOF                          =  0,
        k1020TOF                              ,
        k2030TOF                              ,
        k3040TOF                              ,
        k4050TOF                              ,
        k5060TOF                              ,
        k6070TOF                              ,
        k7080TOF                              ,
        k8090TOF                              ,
        k90100TOF                             ,
        k0010FW                               ,
        k1020FW                               ,
        k2030FW                               ,
        k3040FW                               ,
        k4050FW                               ,
        k5060FW                               ,
        k6070FW                               ,
        k7080FW                               ,
        k8090FW                               ,
        k90100FW                                    
    };

enum eTriggers {
		kGoodVertexClust = 0,
		kGoodVertexCand,
		kGoodSTART,
		kNoPileUpSTART, 
		kNoPileUpMETA,
		kNoPileUpMDC,
		kNoFlashMDC,
		kGoodMDCMult,
		kGoodMDCMIPSMult,
		kGoodLepMult, 
		kGoodTRIGGER,
		kGoodSTART2,
		kNoVETO,
		kGoodSTARTVETO,
		kGoodSTARTMETA, 
		kPT1,
		kPT2,
		kPT3,
		kPT4,
		kNtriggers
};

enum event_cuts{
        cVeretexPositionZ = 0,
	cVertexQuality,
        cVeretexPositionXY,
        cTriggerVertexClust,
        cTriggerVertexCand,
        cTriggerGoodStart,
        cTriggerNoPileUp,
        cTriggerGoodStartVeto,
        cTriggerGoodStartMeta,
        cTriggerNoVeto,
	cNumOfEventCuts,
        cPT2,
        cPT3
};

enum track_cuts{
        cDCA = 0,
        cTrackHitMatchX,
        cTrackHitMatchY,
        cChi2,
        cBeta,
        cNumOfTrackCuts
};

enum eDetectors {
		kMDC_in = 0,
		kMDC_out,
		kMDC_all,
		kMETA,
		kNdetectors
};

enum eCentralityEstimators {
		kNhitsTOF_RPC = 0,
		kNhitsTOF_RPC_cut,
		kNtracks, 
		kNselectedTracks, 
		kNCentralityEstimators
};

Float_t GetCentralityPercentile(Int_t CentralityEstimator, DataTreeEvent* fEvent, TH1F* EstimatorHisto);
Float_t GetNewCentralityEstimator(Int_t CentralityEstimator, DataTreeEvent* fEvent);

void NewCentrality()
{
	TChain* fChain=new TChain("DataTree");
	fChain->Add("/home/vad/NIR codes/Centrality/AuAu_1_23AGev_gen9_108.root");

	DataTreeEvent* fEvent=new DataTreeEvent;
	fChain->SetBranchAddress("DTEvent", &fEvent);	

	TFile f("/home/vad/NIR codes/Centrality/centrality_epcorr_apr12_gen8_2018_07.root","read");
	TH1F* EstimatorHisto[kNumCentralityEstimator];
	for (int i = 0; i < kNumCentralityEstimator; i++) EstimatorHisto[i] = (TH1F*)f.Get(Form("/Centrality/%s_percentile", EstimatorHistoName[i].Data()));	
	
	TFile *f1 = new TFile("NewCentralityPT3Histos.root", "recreate");
	TH1F* Histo1D_newCentrality[20];
	for (int i = 0; i < 10; i++)   Histo1D_newCentrality[i] = new TH1F(Form("%s %.1f%%-%.1f%%", EstimatorHistoName[0].Data(), i*10.0, (i+1)*10.0),";centralityEstimator percentile;counts",100,0.0,100.0);
	for (int i = 10; i < 20; i++)  Histo1D_newCentrality[i] = new TH1F(Form("%s %.1f%%-%.1f%%", EstimatorHistoName[12].Data(), (i-10)*10.0, (i-9)*10.0),";centralityEstimator percentile;counts",100,0.0,100.0); 
	for (int i = 0; i < 9; i++)   Histo1D_newCentrality[i] -> SetLineColor(i+1);
	Histo1D_newCentrality[9]  -> SetLineColor(28);
	for (int i = 10; i < 19; i++) Histo1D_newCentrality[i] -> SetLineColor(i-9);
	Histo1D_newCentrality[19] -> SetLineColor(28);
	for (int i = 0; i < 20; i++) Histo1D_newCentrality[i]  -> SetLineWidth(4);
	
	Long64_t lNEvents = fChain->GetEntries();
	Float_t fNHitsTOF;
        Float_t fNTracksMDC;
        Float_t fChargeFW;
        Float_t fVertexPosition[3];
        DataTreeTrack* fTrack;
	DataTreeTOFHit* fHit;
	
	for (long i = 0; i < lNEvents; i++)
	{
		fChain->GetEntry(i);
		fNTracksMDC =    fEvent->GetCentralityEstimator(kNselectedTracks);
        	fNHitsTOF =      fEvent->GetCentralityEstimator(kNhitsTOF_RPC_cut);
        	fChargeFW =      fEvent->GetPSDEnergy();
		for(int j=0;j<3;j++) fVertexPosition[j] = fEvent->GetVertexPositionComponent(j);

		//start of event selection
		Bool_t EVENT_FLAG=false;		

		if (!fEvent->GetTrigger(16)->GetIsFired())                                         EVENT_FLAG=true;   
	    	if (fVertexPosition[2]>0 || fVertexPosition[2]<-60)                                EVENT_FLAG=true;
	    	if (fVertexPosition[0]*fVertexPosition[0]+fVertexPosition[1]*fVertexPosition[1]>9) EVENT_FLAG=true;
	    	if (fEvent->GetVertexQuality()<0.5 || fEvent->GetVertexQuality()>40)               EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodVertexClust)->GetIsFired())                           EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodVertexCand)->GetIsFired())                            EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodSTART)->GetIsFired())                                 EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kNoPileUpSTART)->GetIsFired())                             EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodSTARTVETO)->GetIsFired())                             EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodSTARTMETA)->GetIsFired())                             EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kNoVETO)->GetIsFired())                                    EVENT_FLAG=true;

		if(EVENT_FLAG) continue;
		//end of event selection
	
		Float_t CentralityPercentileTOF = GetCentralityPercentile(kTOFRPC, fEvent, EstimatorHisto[kTOFRPC]);
		Float_t CentralityPercentileFW  = GetCentralityPercentile(kFWSumChargeSpec, fEvent, EstimatorHisto[kFWSumChargeSpec]);

		for (int j = 0; j < 10; j++)  {
			if(CentralityPercentileTOF>=(j*10.0) && CentralityPercentileTOF<((j+1)*10.0)) Histo1D_newCentrality[j] -> Fill(CentralityPercentileTOF);
			else continue;
			}
		for (int j = 10; j < 20; j++) {
			if(CentralityPercentileFW>=((j-10)*10.0) && CentralityPercentileFW<((j-9)*10.0)) Histo1D_newCentrality[j] -> Fill(CentralityPercentileFW);
			else continue;
			}
	}	
	f.Close();
		
	TCanvas* c1 = new TCanvas("newCentralityTOF_PT3","newCentralityTOF_PT3");
	TLegend legend1(0.1, 0.2, 0.3, 0.4, "newCentralityTOF_PT3");
	for (int i = 0; i < 10; i++) legend1.AddEntry(Histo1D_newCentrality[i], Form("%s %.1f%%-%.1f%%", EstimatorHistoName[0].Data(), i*10.0, (i+1)*10.0));
	THStack *newCentralityTOF = new THStack("newCentralityTOF_PT3","newCentralityTOF_PT3");
	for (int i = 0; i < 10; i++) newCentralityTOF -> Add(Histo1D_newCentrality[i]);
	newCentralityTOF -> Draw("nostack");
	legend1.Draw("same");
	c1->Write();
	for (int i = 0; i < 10; i++) Histo1D_newCentrality[i] -> Write();
	
	TCanvas* c2 = new TCanvas("newCentralityFW_PT3","newCentralityFW_PT3");
	TLegend legend2(0.1, 0.2, 0.3, 0.4, "newCentralityFW_PT3");
	for (int i = 10; i < 20; i++) legend2.AddEntry(Histo1D_newCentrality[i], Form("%s %.1f%%-%.1f%%", EstimatorHistoName[12].Data(), (i-10)*10.0, (i-9)*10.0));
	THStack *newCentralityFW = new THStack("newCentralityFW_PT3","newCentralityFW_PT3");
	for (int i = 10; i < 20; i++) newCentralityFW -> Add(Histo1D_newCentrality[i]);
	newCentralityFW -> Draw("nostack");
	legend2.Draw("same");
	c2->Write(); 
	for (int i = 10; i < 20; i++) Histo1D_newCentrality[i] -> Write();
	
	f1->Close();
	
}

Float_t GetCentralityPercentile(Int_t CentralityEstimator, DataTreeEvent* fEvent, TH1F* EstimatorHisto)
{

    if(CentralityEstimator>=kNumCentralityEstimator){
        return 101.;
    }
    else if(!EstimatorHisto){
        return 101.;
    }
    else {
        Int_t meaning = GetNewCentralityEstimator(CentralityEstimator, fEvent);
        Int_t bin = EstimatorHisto->FindBin(meaning);
        return EstimatorHisto->GetBinContent(bin);
    }

    return 101.;
}

Float_t GetNewCentralityEstimator(Int_t CentralityEstimator, DataTreeEvent* fEvent)
{
    if(CentralityEstimator==kTOFRPC) return fEvent->GetCentralityEstimator(kNhitsTOF_RPC_cut);
    else if(CentralityEstimator==kFWSumChargeSpec) return fEvent->GetPSDEnergy();
    else return 0;
}
