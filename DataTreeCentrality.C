enum eCebtralityClass { 
        k0005TOF                          =  0,
        k0510TOF                              ,
        k1015TOF                              ,
        k1520TOF                              ,
        k2025TOF                              ,
        k2530TOF                              ,
        k3035TOF                              ,
        k3540TOF                              ,
        k4045TOF                              ,
        k4550TOF                              ,
        k0005MDC                              ,
        k0510MDC                              ,
        k1015MDC                              ,
        k1520MDC                              ,
        k2025MDC                              ,
        k2530MDC                              ,
        k3035MDC                              ,
        k3540MDC                              ,
        k4045MDC                              ,
        k4550MDC                                   
    };

TString EstimatorName [2] = {"TOF", "MDC"};

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

void DataTreeCentrality()
{
	TChain* fChain=new TChain("DataTree");
	fChain->Add("/home/vad/NIR codes/Centrality/AuAu_1_23AGev_gen9_108.root");
	DataTreeEvent* fEvent=new DataTreeEvent;
	fChain->SetBranchAddress("DTEvent", &fEvent);
	
	TH1F* Histo1D_datatreeCentrality[20];
	for (int i = 0; i < 10; i++)   Histo1D_datatreeCentrality[i] = new TH1F(Form("%s %.1f%%-%.1f%%", EstimatorName[0].Data(), i*5.0, i*5.0+5.0), ";centralityEstimator;counts",50,0.0,100.0);
	for (int i = 10; i < 20; i++)  Histo1D_datatreeCentrality[i] = new TH1F(Form("%s %.1f%%-%.1f%%", EstimatorName[1].Data(), (i-10)*5.0, (i-9)*5.0), ";centralityEstimator;counts",50,0.0,100.0);
	for (int i = 0; i < 9; i++)   Histo1D_datatreeCentrality[i] -> SetLineColor(i+1);
	Histo1D_datatreeCentrality[9]  -> SetLineColor(28);
	for (int i = 10; i < 19; i++) Histo1D_datatreeCentrality[i] -> SetLineColor(i-9);
	Histo1D_datatreeCentrality[19] -> SetLineColor(28);
	for (int i = 0; i < 20; i++) Histo1D_datatreeCentrality[i] -> SetLineWidth(4);

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

		if (!fEvent->GetTrigger(17) -> GetIsFired())                                         EVENT_FLAG=true;   
	    	if (fVertexPosition[2]>0 || fVertexPosition[2]<-60)                                  EVENT_FLAG=true;
	    	if (fVertexPosition[0]*fVertexPosition[0]+fVertexPosition[1]*fVertexPosition[1]>9)   EVENT_FLAG=true;
	    	if (fEvent->GetVertexQuality()<0.5 || fEvent->GetVertexQuality()>40)                 EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodVertexClust) -> GetIsFired())                           EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodVertexCand) -> GetIsFired())                            EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodSTART) -> GetIsFired())                                 EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kNoPileUpSTART) -> GetIsFired())                             EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodSTARTVETO) -> GetIsFired())                             EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kGoodSTARTMETA) -> GetIsFired())                             EVENT_FLAG=true;
	    	if (!fEvent->GetTrigger(kNoVETO) -> GetIsFired())                                    EVENT_FLAG=true;

		if(EVENT_FLAG) continue;
		//end of event selection
				
		Float_t CentralityPercentileTOF = fEvent->GetCentralityEstimator(kNhitsTOF_RPC_cut);
		Float_t CentralityPercentileMDC = fEvent->GetCentralityEstimator(kNselectedTracks);
		Float_t Centrality              = fEvent->GetCentrality();

		for (int i = 0; i < 10; i++)  {
			if(Centrality>=(i*5.0) && Centrality<((i+1)*5.0)) Histo1D_datatreeCentrality[i] -> Fill(CentralityPercentileTOF);
			else continue;
			}
		for (int i = 10; i < 20; i++) {
			if(Centrality>=((i-10)*5.0) && Centrality<((i-9)*5.0)) Histo1D_datatreeCentrality[i] -> Fill(CentralityPercentileMDC);
			else continue;
			}

	}

	TFile *f1 = new TFile("DataTreeCentralityPT3Histos.root", "recreate");

	TCanvas* c1 = new TCanvas("datatreeCentralityTOF_PT3","datatreeCentralityTOF_PT3");
	TLegend legend1(0.1, 0.2, 0.3, 0.4, "datatreeCentralityTOF_PT3");
	for (int i = 0; i < 10; i++) legend1.AddEntry(Histo1D_datatreeCentrality[i], Form ("%s %.1f%%-%.1f%%", EstimatorName[0].Data(), i*5.0, i*5.0+5.0));
	THStack *datatreeCentralityTOF = new THStack("datatreeCentralityTOF_PT3","datatreeCentralityTOF_PT3");
	for (int i = 0; i < 10; i++) datatreeCentralityTOF -> Add(Histo1D_datatreeCentrality[i]);
	datatreeCentralityTOF -> Draw("nostack");
	legend1.Draw("same");
	c1->Write();
	for (int i = 0; i < 10; i++) Histo1D_datatreeCentrality[i] -> Write();

	TCanvas* c2 = new TCanvas("datatreeCentralityMDC_PT3","datatreeCentralityMDC_PT3");
	TLegend legend2(0.1, 0.2, 0.3, 0.4, "datatreeCentralityMDC_PT3");
	for (int i = 10; i < 20; i++) legend2.AddEntry(Histo1D_datatreeCentrality[i], Form ("%s %.1f%%-%.1f%%", EstimatorName[1].Data(), (i-10)*5.0, (i-9)*5.0));
	THStack *datatreeCentralityMDC = new THStack("datatreeCentralityMDC_PT3","datatreeCentralityMDC_PT3");
	for (int i = 10; i < 20; i++) datatreeCentralityMDC -> Add(Histo1D_datatreeCentrality[i]);
	datatreeCentralityMDC -> Draw("nostack");
	legend2.Draw("same");
	c2->Write();
	for (int i = 10; i < 20; i++) Histo1D_datatreeCentrality[i] -> Write();
	
	f1->Close();
	
}
