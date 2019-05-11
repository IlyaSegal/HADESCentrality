enum eCentralityEstimator {
        kTOFRPC                          =  0,
        kTOF                                 ,
        kRPC                                 ,
        kMDC                                 ,
	kFW                                  ,
        kCentralityEstimator              
    };

TString CentralityEstimatorName[kCentralityEstimator]=
    {
        "hitsTOF+RPC"                         ,
        "hitsTOF"                             ,
        "hitsRPC"                             ,
        "tracksMDC"                           ,
	"FWSumChargeZ"
    };

void CentralityClasses(Int_t CentralityClasses){
	TFile *f1 = new TFile("/home/vad/centrality/build/Result/HistoCutResult_TOF.root");
	TFile *f2 = new TFile("/home/vad/centrality/build/Result/HistoCutResult_RPC.root");
	TFile *f3 = new TFile("/home/vad/centrality/build/Result/HistoCutResult_MDC.root");
	TFile *f4 = new TFile("/home/vad/centrality/build/Result/HistoCutResult_FW.root");
	TFile *f5 = new TFile("/home/vad/NIR_codes/Centrality/gmc-Au3Au3-snn23.6-md0.4-nd-1.0-rc1-smax99.0_1000000.root");
	TFile *f6 = new TFile("/home/vad/NIR_codes/QAHistoBuilder/QASelectedPT2Histos.root");
	TFile *f7 = new TFile("/home/vad/NIR_codes/Centrality/2.root");
	TFile *f8 = new TFile("/home/vad/NIR_codes/Centrality/1.root");
	TFile *f  = new TFile("FINAL.root", "recreate");
	TTree *Borders_TOF=(TTree*)f1->Get(Form("Borders_%s", CentralityEstimatorName[1].Data()));
	TTree *Borders_RPC=(TTree*)f2->Get(Form("Borders_%s", CentralityEstimatorName[2].Data()));
	TTree *Borders_MDC=(TTree*)f3->Get(Form("Borders_%s", CentralityEstimatorName[3].Data()));
	TTree *Borders_FW=(TTree*)f4->Get(Form("Borders_%s", CentralityEstimatorName[4].Data()));
	TTree *tGlauber=(TTree*)f8->Get("nt_Au3_Au3");
	TTree *tGlauber_FW=(TTree*)f8->Get("nt_Au3_Au3");
	TTree *HADES_FW=(TTree*)f7->Get("HADES_FW");
	TTree *Result=new TTree("Result", "Result");

	Int_t Ncc, MinBorder_TOF, MaxBorder_TOF,  MinBorder_RPC, MaxBorder_RPC, MinBorder_MDC, MaxBorder_MDC, MinBorder_TOF_R, MaxBorder_TOF_R,  MinBorder_RPC_R, MaxBorder_RPC_R, MinBorder_MDC_R, MaxBorder_MDC_R, MinBorder_FW_R, MaxBorder_FW_R; 
	Int_t MinBorder_FW, MaxBorder_FW;
	Float_t MinPercent, MaxPercent, B, Npart, Ncoll, NprotonsA; 
	double B_av, Npart_av, Ncoll_av, B_TOF_av, Npart_TOF_av, Ncoll_TOF_av, B_RPC_av, Npart_RPC_av, Ncoll_RPC_av, B_MDC_av, Npart_MDC_av, Ncoll_MDC_av, B_FW_av, Npart_FW_av, Ncoll_FW_av;
	double RMS_B_av, RMS_Npart_av, RMS_Ncoll_av, RMS_B_TOF_av, RMS_Npart_TOF_av, RMS_Ncoll_TOF_av, RMS_B_RPC_av, RMS_Npart_RPC_av, RMS_Ncoll_RPC_av, RMS_B_MDC_av, RMS_Npart_MDC_av, RMS_Ncoll_MDC_av, RMS_B_FW_av, RMS_Npart_FW_av, RMS_Ncoll_FW_av;
	
	Borders_TOF -> SetBranchAddress("MinBorder", &MinBorder_TOF);
	Borders_TOF -> SetBranchAddress("MaxBorder", &MaxBorder_TOF);
	Int_t entries_TOF=(Int_t)Borders_TOF->GetEntries();
	Borders_RPC -> SetBranchAddress("MinBorder", &MinBorder_RPC);
	Borders_RPC -> SetBranchAddress("MaxBorder", &MaxBorder_RPC);
	Int_t entries_RPC=(Int_t)Borders_RPC->GetEntries();
	Borders_MDC -> SetBranchAddress("MinBorder", &MinBorder_MDC);
	Borders_MDC -> SetBranchAddress("MaxBorder", &MaxBorder_MDC);
	Int_t entries_MDC=(Int_t)Borders_MDC->GetEntries();
	Borders_FW  -> SetBranchAddress("MinBorder", &MinBorder_FW);
	Borders_FW  -> SetBranchAddress("MaxBorder", &MaxBorder_FW);
	Int_t entries_FW=(Int_t)Borders_FW->GetEntries();
	tGlauber    -> SetBranchAddress("Npart", &Npart);
	tGlauber    -> SetBranchAddress("Ncoll", &Ncoll);
	tGlauber    -> SetBranchAddress("B", &B);
	Int_t entries_glauber=(Int_t)tGlauber->GetEntries();
	tGlauber_FW -> SetBranchAddress("Npart", &Npart);
	tGlauber_FW -> SetBranchAddress("Ncoll", &Ncoll);
	tGlauber_FW -> SetBranchAddress("B", &B);
	Int_t entries_glauber_FW=(Int_t)tGlauber_FW->GetEntries();
	HADES_FW    -> SetBranchAddress("NprotonsA", &NprotonsA);
	Int_t entries_fw=(Int_t)HADES_FW->GetEntries();

	Result -> Branch("Ncc", &Ncc);
	Result -> Branch("MinPercent", &MinPercent);
	Result -> Branch("MaxPercent", &MaxPercent);
	Result -> Branch("MinBorder_TOF", &MinBorder_TOF_R);
	Result -> Branch("MaxBorder_TOF", &MaxBorder_TOF_R);
	Result -> Branch("MinBorder_RPC", &MinBorder_RPC_R);
	Result -> Branch("MaxBorder_RPC", &MaxBorder_RPC_R); 
	Result -> Branch("MinBorder_MDC", &MinBorder_MDC_R);
	Result -> Branch("MaxBorder_MDC", &MaxBorder_MDC_R); 
	Result -> Branch("MinBorder_FW", &MinBorder_FW_R);
	Result -> Branch("MaxBorder_FW", &MaxBorder_FW_R);
	Result -> Branch("B_average", &B_av); 
	Result -> Branch("Npart_average", &Npart_av);
	Result -> Branch("Ncoll_average", &Ncoll_av);
	Result -> Branch("B_TOF_average", &B_TOF_av);
	Result -> Branch("Npart_TOF_average", &Npart_TOF_av);
	Result -> Branch("Ncoll_TOF_average", &Ncoll_TOF_av);
	Result -> Branch("B_RPC_average", &B_RPC_av);
	Result -> Branch("Npart_RPC_average", &Npart_RPC_av);
	Result -> Branch("Ncoll_RPC_average", &Ncoll_RPC_av);
	Result -> Branch("B_MDC_average", &B_MDC_av);
	Result -> Branch("Npart_MDC_average", &Npart_MDC_av);
	Result -> Branch("Ncoll_MDC_average", &Ncoll_MDC_av);
	Result -> Branch("B_FW_average", &B_FW_av); 
	Result -> Branch("Npart_FW_average", &Npart_FW_av);
	Result -> Branch("Ncoll_FW_average", &Ncoll_FW_av);
	Result -> Branch("RMS_B_average", &RMS_B_av); 
	Result -> Branch("RMS_Npart_average", &RMS_Npart_av);
	Result -> Branch("RMS_Ncoll_average", &RMS_Ncoll_av);
	Result -> Branch("RMS_B_TOF_average", &RMS_B_TOF_av);
	Result -> Branch("RMS_Npart_TOF_average", &RMS_Npart_TOF_av);
	Result -> Branch("RMS_Ncoll_TOF_average", &RMS_Ncoll_TOF_av);
	Result -> Branch("RMS_B_RPC_average", &RMS_B_RPC_av);
	Result -> Branch("RMS_Npart_RPC_average", &RMS_Npart_RPC_av);
	Result -> Branch("RMS_Ncoll_RPC_average", &RMS_Ncoll_RPC_av);
	Result -> Branch("RMS_B_MDC_average", &RMS_B_MDC_av);
	Result -> Branch("RMS_Npart_MDC_average", &RMS_Npart_MDC_av);
	Result -> Branch("RMS_Ncoll_MDC_average", &RMS_Ncoll_MDC_av);
	Result -> Branch("RMS_B_FW_average", &RMS_B_FW_av); 
	Result -> Branch("RMS_Npart_FW_average", &RMS_Npart_FW_av);
	Result -> Branch("RMS_Ncoll_FW_average", &RMS_Ncoll_FW_av);

	TH1F* B_VS_CentralityHisto_TOF[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   B_VS_CentralityHisto_TOF[i] = new TH1F(Form("B_VS_CentralityClass_TOF %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;B, fm;counts", 200, 0, 20);
	int j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		B_VS_CentralityHisto_TOF[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) B_VS_CentralityHisto_TOF[i] -> SetLineWidth(4);

	TH1F* Npart_VS_CentralityHisto_TOF[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Npart_VS_CentralityHisto_TOF[i] = new TH1F(Form("Npart_VS_CentralityClass_TOF %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Npart;counts", 400, 0, 400);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Npart_VS_CentralityHisto_TOF[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Npart_VS_CentralityHisto_TOF[i] -> SetLineWidth(4);

	TH1F* Ncoll_VS_CentralityHisto_TOF[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Ncoll_VS_CentralityHisto_TOF[i] = new TH1F(Form("Ncoll_VS_CentralityClass_TOF %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Ncoll;counts", 750, 0, 750);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Ncoll_VS_CentralityHisto_TOF[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Ncoll_VS_CentralityHisto_TOF[i] -> SetLineWidth(4);

	TH1F* B_VS_CentralityHisto_RPC[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   B_VS_CentralityHisto_RPC[i] = new TH1F(Form("B_VS_CentralityClass_RPC %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;B, fm;counts", 200, 0, 20);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		B_VS_CentralityHisto_RPC[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) B_VS_CentralityHisto_RPC[i] -> SetLineWidth(4);

	TH1F* Npart_VS_CentralityHisto_RPC[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Npart_VS_CentralityHisto_RPC[i] = new TH1F(Form("Npart_VS_CentralityClass_RPC %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Npart;counts", 400, 0, 400);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Npart_VS_CentralityHisto_RPC[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Npart_VS_CentralityHisto_RPC[i] -> SetLineWidth(4);

	TH1F* Ncoll_VS_CentralityHisto_RPC[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Ncoll_VS_CentralityHisto_RPC[i] = new TH1F(Form("Ncoll_VS_CentralityClass_RPC %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Ncoll;counts", 750, 0, 750);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Ncoll_VS_CentralityHisto_RPC[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Ncoll_VS_CentralityHisto_RPC[i] -> SetLineWidth(4);

	TH1F* B_VS_CentralityHisto_MDC[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   B_VS_CentralityHisto_MDC[i] = new TH1F(Form("B_VS_CentralityClass_MDC %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;B, fm;counts", 200, 0, 20);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		B_VS_CentralityHisto_MDC[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) B_VS_CentralityHisto_MDC[i] -> SetLineWidth(4);

	TH1F* Npart_VS_CentralityHisto_MDC[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Npart_VS_CentralityHisto_MDC[i] = new TH1F(Form("Npart_VS_CentralityClass_MDC %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Npart;counts", 400, 0, 400);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Npart_VS_CentralityHisto_MDC[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Npart_VS_CentralityHisto_MDC[i] -> SetLineWidth(4);

	TH1F* Ncoll_VS_CentralityHisto_MDC[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Ncoll_VS_CentralityHisto_MDC[i] = new TH1F(Form("Ncoll_VS_CentralityClass_MDC %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Ncoll;counts", 750, 0, 750);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Ncoll_VS_CentralityHisto_MDC[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Ncoll_VS_CentralityHisto_MDC[i] -> SetLineWidth(4);

	TH1F* B_VS_CentralityHisto_FW[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   B_VS_CentralityHisto_FW[i] = new TH1F(Form("B_VS_CentralityClass_FW %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;B, fm;counts", 200, 0, 20);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		B_VS_CentralityHisto_FW[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) B_VS_CentralityHisto_FW[i] -> SetLineWidth(4);

	TH1F* Npart_VS_CentralityHisto_FW[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Npart_VS_CentralityHisto_FW[i] = new TH1F(Form("Npart_VS_CentralityClass_FW %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Npart;counts", 400, 0, 400);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Npart_VS_CentralityHisto_FW[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Npart_VS_CentralityHisto_FW[i] -> SetLineWidth(4);

	TH1F* Ncoll_VS_CentralityHisto_FW[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   Ncoll_VS_CentralityHisto_FW[i] = new TH1F(Form("Ncoll_VS_CentralityClass_FW %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;Ncoll;counts", 750, 0, 750);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		Ncoll_VS_CentralityHisto_FW[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) Ncoll_VS_CentralityHisto_FW[i] -> SetLineWidth(4);

	TH1F* hitsTOF_VS_CentralityHisto[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   hitsTOF_VS_CentralityHisto[i] = new TH1F(Form("hitsTOF_VS_CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;hitsTOF;counts", 100, 0, 100);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		hitsTOF_VS_CentralityHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) hitsTOF_VS_CentralityHisto[i] -> SetLineWidth(4);

	TH1F* hitsRPC_VS_CentralityHisto[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   hitsRPC_VS_CentralityHisto[i] = new TH1F(Form("hitsRPC_VS_CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;hitsRPC;counts", 200, 0, 200);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		hitsRPC_VS_CentralityHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) hitsRPC_VS_CentralityHisto[i] -> SetLineWidth(4);

	TH1F* tracksMDC_VS_CentralityHisto[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   tracksMDC_VS_CentralityHisto[i] = new TH1F(Form("tracksMDC_VS_CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;tracksMDC;counts", 150, 0, 150);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		tracksMDC_VS_CentralityHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) tracksMDC_VS_CentralityHisto[i] -> SetLineWidth(4);

	TH1F* FWSumChargeZ_VS_CentralityHisto[CentralityClasses+1];
	for (int i = 0; i < CentralityClasses+1; i++)   FWSumChargeZ_VS_CentralityHisto[i] = new TH1F(Form("FWSumChargeZ_VS_CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";;FWSumChargeZ;counts", 100, 0, 100);
	j=1; 
	for (int i = 0; i < CentralityClasses+1; i++) {
		FWSumChargeZ_VS_CentralityHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses+1; i++) FWSumChargeZ_VS_CentralityHisto[i] -> SetLineWidth(4);

	TH1F* B_average_VS_Centrality         = new TH1F("B_average_VS_Centrality", ";Centrality Percent;B, fm", CentralityClasses, 0, 100);
	TH1F* Npart_average_VS_Centrality     = new TH1F("Npart_average_VS_Centrality", ";Centrality Percent;Npart", CentralityClasses, 0, 100);
	TH1F* Ncoll_average_VS_Centrality     = new TH1F("Ncoll_average_VS_Centrality", ";Centrality Percent;Ncoll", CentralityClasses, 0, 100);	
	TH1F* B_average_VS_Centrality_TOF     = new TH1F("B_average_VS_Centrality_TOF", "Centrality Percent;B, fm;", CentralityClasses, 0, 100);
	TH1F* Npart_average_VS_Centrality_TOF = new TH1F("Npart_average_VS_Centrality_TOF", ";Centrality Percent;Npart", CentralityClasses, 0, 100);
	TH1F* Ncoll_average_VS_Centrality_TOF = new TH1F("Ncoll_average_VS_Centrality_TOF", ";Centrality Percent;Ncoll", CentralityClasses, 0, 100);
	TH1F* B_average_VS_Centrality_RPC     = new TH1F("B_average_VS_Centrality_RPC", ";Centrality Percent;B, fm", CentralityClasses, 0, 100);
	TH1F* Npart_average_VS_Centrality_RPC = new TH1F("Npart_average_VS_Centrality_RPC", ";Centrality Percent;Npart", CentralityClasses, 0, 100);
	TH1F* Ncoll_average_VS_Centrality_RPC = new TH1F("Ncoll_average_VS_Centrality_RPC", "l;Centrality Percent;Ncol", CentralityClasses, 0, 100);
	TH1F* B_average_VS_Centrality_MDC     = new TH1F("B_average_VS_Centrality_MDC", ";Centrality Percent;B, fm", CentralityClasses, 0, 100);
	TH1F* Npart_average_VS_Centrality_MDC = new TH1F("Npart_average_VS_Centrality_MDC", ";Centrality Percent;Npart", CentralityClasses, 0, 100);
	TH1F* Ncoll_average_VS_Centrality_MDC = new TH1F("Ncoll_average_VS_Centrality_MDC", ";Centrality Percent;Ncoll", CentralityClasses, 0, 100);
	TH1F* B_average_VS_Centrality_FW      = new TH1F("B_average_VS_Centrality_FW ", ";Centrality Percent;B, fm", CentralityClasses, 0, 100);
	TH1F* Npart_average_VS_Centrality_FW  = new TH1F("Npart_average_VS_Centrality_FW ", ";Centrality Percent;Npart", CentralityClasses, 0, 100);
	TH1F* Ncoll_average_VS_Centrality_FW  = new TH1F("Ncoll_average_VS_Centrality_FW ", ";Centrality Percent;Ncoll", CentralityClasses, 0, 100);

	B_average_VS_Centrality         -> SetLineColor(1);
	Npart_average_VS_Centrality     -> SetLineColor(1);
	Ncoll_average_VS_Centrality     -> SetLineColor(1);	
	B_average_VS_Centrality_TOF     -> SetLineColor(2);
	Npart_average_VS_Centrality_TOF -> SetLineColor(2);
	Ncoll_average_VS_Centrality_TOF -> SetLineColor(2);
	B_average_VS_Centrality_RPC     -> SetLineColor(3);
	Npart_average_VS_Centrality_RPC -> SetLineColor(3);
	Ncoll_average_VS_Centrality_RPC -> SetLineColor(3);
	B_average_VS_Centrality_MDC     -> SetLineColor(4);
	Npart_average_VS_Centrality_MDC -> SetLineColor(4);
	Ncoll_average_VS_Centrality_MDC -> SetLineColor(4);
	B_average_VS_Centrality_FW      -> SetLineColor(5);
	Npart_average_VS_Centrality_FW  -> SetLineColor(5);
	Ncoll_average_VS_Centrality_FW  -> SetLineColor(5);

	B_average_VS_Centrality         -> SetLineWidth(4);
	Npart_average_VS_Centrality     -> SetLineWidth(4);
	Ncoll_average_VS_Centrality     -> SetLineWidth(4);	
	B_average_VS_Centrality_TOF     -> SetLineWidth(4);
	Npart_average_VS_Centrality_TOF -> SetLineWidth(4);
	Ncoll_average_VS_Centrality_TOF -> SetLineWidth(4);
	B_average_VS_Centrality_RPC     -> SetLineWidth(4);
	Npart_average_VS_Centrality_RPC -> SetLineWidth(4);
	Ncoll_average_VS_Centrality_RPC -> SetLineWidth(4);
	B_average_VS_Centrality_MDC     -> SetLineWidth(4);
	Npart_average_VS_Centrality_MDC -> SetLineWidth(4);
	Ncoll_average_VS_Centrality_MDC -> SetLineWidth(4);
	B_average_VS_Centrality_FW      -> SetLineWidth(4);
	Npart_average_VS_Centrality_FW  -> SetLineWidth(4);
	Ncoll_average_VS_Centrality_FW  -> SetLineWidth(4);

	B_average_VS_Centrality         -> SetMarkerStyle(20);
	Npart_average_VS_Centrality     -> SetMarkerStyle(20);
	Ncoll_average_VS_Centrality     -> SetMarkerStyle(20);	
	B_average_VS_Centrality_TOF     -> SetMarkerStyle(20);
	Npart_average_VS_Centrality_TOF -> SetMarkerStyle(20);
	Ncoll_average_VS_Centrality_TOF -> SetMarkerStyle(20);
	B_average_VS_Centrality_RPC     -> SetMarkerStyle(20);
	Npart_average_VS_Centrality_RPC -> SetMarkerStyle(20);
	Ncoll_average_VS_Centrality_RPC -> SetMarkerStyle(20);
	B_average_VS_Centrality_MDC     -> SetMarkerStyle(20);
	Npart_average_VS_Centrality_MDC -> SetMarkerStyle(20);
	Ncoll_average_VS_Centrality_MDC -> SetMarkerStyle(20);
	B_average_VS_Centrality_FW      -> SetMarkerStyle(20);
	Npart_average_VS_Centrality_FW  -> SetMarkerStyle(20);
	Ncoll_average_VS_Centrality_FW  -> SetMarkerStyle(20);

	B_average_VS_Centrality         -> SetMarkerSize(2);
	Npart_average_VS_Centrality     -> SetMarkerSize(2);
	Ncoll_average_VS_Centrality     -> SetMarkerSize(2);	
	B_average_VS_Centrality_TOF     -> SetMarkerSize(2);
	Npart_average_VS_Centrality_TOF -> SetMarkerSize(2);
	Ncoll_average_VS_Centrality_TOF -> SetMarkerSize(2);
	B_average_VS_Centrality_RPC     -> SetMarkerSize(2);
	Npart_average_VS_Centrality_RPC -> SetMarkerSize(2);
	Ncoll_average_VS_Centrality_RPC -> SetMarkerSize(2);
	B_average_VS_Centrality_MDC     -> SetMarkerSize(2);
	Npart_average_VS_Centrality_MDC -> SetMarkerSize(2);
	Ncoll_average_VS_Centrality_MDC -> SetMarkerSize(2);
	B_average_VS_Centrality_FW      -> SetMarkerSize(2);
	Npart_average_VS_Centrality_FW  -> SetMarkerSize(2);
	Ncoll_average_VS_Centrality_FW  -> SetMarkerSize(2);

	B_average_VS_Centrality         -> SetMarkerColor(1);
	Npart_average_VS_Centrality     -> SetMarkerColor(1);
	Ncoll_average_VS_Centrality     -> SetMarkerColor(1);	
	B_average_VS_Centrality_TOF     -> SetMarkerColor(2);
	Npart_average_VS_Centrality_TOF -> SetMarkerColor(2);
	Ncoll_average_VS_Centrality_TOF -> SetMarkerColor(2);
	B_average_VS_Centrality_RPC     -> SetMarkerColor(3);
	Npart_average_VS_Centrality_RPC -> SetMarkerColor(3);
	Ncoll_average_VS_Centrality_RPC -> SetMarkerColor(3);
	B_average_VS_Centrality_MDC     -> SetMarkerColor(4);
	Npart_average_VS_Centrality_MDC -> SetMarkerColor(4);
	Ncoll_average_VS_Centrality_MDC -> SetMarkerColor(4);
	B_average_VS_Centrality_FW      -> SetMarkerColor(5);
	Npart_average_VS_Centrality_FW  -> SetMarkerColor(5);
	Ncoll_average_VS_Centrality_FW  -> SetMarkerColor(5);

	std::unique_ptr<TTree> t {(TTree*)f8->Get("nt_Au3_Au3")};
	Glauber::Fitter *fitter=new Glauber::Fitter(move(t), "Default");
	tGlauber    -> SetBranchAddress("Npart", &Npart);
	tGlauber    -> SetBranchAddress("Ncoll", &Ncoll);
	tGlauber    -> SetBranchAddress("B", &B);

	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f6->Get("hitsTOF_selected")));
	fitter -> Glauber::Fitter::Init(entries_glauber, "Default");
	double alpha=1.64e-6;
	for (int i = 0; i < entries_glauber; i++) {
		cout<<i<<endl;		
		tGlauber -> GetEntry(i);
		
		double v=1.0;
		int Na = int(v*Npart+(1-v)*Ncoll);
		fitter->Glauber::Fitter::SetNBDhist(0.21, 8.0);
                std::unique_ptr<TH1F> htemp {(TH1F*)((fitter->Glauber::Fitter::GetNBDHisto()).Clone("htemp"))};
        	float nHits {0.};
        	for (int j=0; j<Na; j++) nHits += (1-alpha*Npart*Npart)*(htemp->GetRandom());	
		cout<<"alpha="<<(1-alpha*Npart*Npart)<<"      nhits_TOF="<<nHits<<endl;
        	for (int j=0; j<entries_TOF; j++) {
			Borders_TOF -> GetEntry(j);
			if (nHits>=MinBorder_TOF && nHits<=MaxBorder_TOF) {
				cout<<"true"<<endl;
				B_VS_CentralityHisto_TOF[j]     -> Fill(B);
				Npart_VS_CentralityHisto_TOF[j] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_TOF[j] -> Fill(Ncoll);
				hitsTOF_VS_CentralityHisto[j] -> Fill(nHits);
				B_VS_CentralityHisto_TOF[CentralityClasses]     -> Fill(B);
				Npart_VS_CentralityHisto_TOF[CentralityClasses] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_TOF[CentralityClasses] -> Fill(Ncoll);
				hitsTOF_VS_CentralityHisto[CentralityClasses] -> Fill(nHits);
				break;
				}
			}
		}

	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f6->Get("hitsRPC_selected")));
	fitter -> Glauber::Fitter::Init(entries_glauber, "Default");
	for (int i = 0; i < entries_glauber; i++) {
		cout<<i<<endl;		
		tGlauber -> GetEntry(i);
		
		double v=1.0;
		int Na = int(v*Npart+(1-v)*Ncoll);
		fitter->Glauber::Fitter::SetNBDhist(0.49, 31.0);
                std::unique_ptr<TH1F> htemp {(TH1F*)((fitter->Glauber::Fitter::GetNBDHisto()).Clone("htemp"))};
        	float nHits {0.};
        	for (int j=0; j<Na; j++) nHits += (1-alpha*Npart*Npart)*(htemp->GetRandom());	
		cout<<"alpha="<<(1-alpha*Npart*Npart)<<"      nhits_RPC="<<nHits<<endl;
        	for (int j=0; j<entries_RPC; j++) {
			Borders_RPC -> GetEntry(j);
			if (nHits>=MinBorder_RPC && nHits<=MaxBorder_RPC) {
				cout<<"true"<<endl;
				B_VS_CentralityHisto_RPC[j]     -> Fill(B);
				Npart_VS_CentralityHisto_RPC[j] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_RPC[j] -> Fill(Ncoll);
				hitsRPC_VS_CentralityHisto[j] -> Fill(nHits);
				B_VS_CentralityHisto_RPC[CentralityClasses]     -> Fill(B);
				Npart_VS_CentralityHisto_RPC[CentralityClasses] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_RPC[CentralityClasses] -> Fill(Ncoll);
				hitsRPC_VS_CentralityHisto[CentralityClasses] -> Fill(nHits);
				break;
				}
			}
		}

	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f6->Get("tracksMDC_selected")));
	fitter -> Glauber::Fitter::Init(entries_glauber, "Default");
	alpha=1.10e-7;
	for (int i = 0; i < entries_glauber; i++) {
		cout<<i<<endl;		
		tGlauber -> GetEntry(i);
		
		double v=1.0;
		int Na = int(v*Npart+(1-v)*Ncoll);
		fitter->Glauber::Fitter::SetNBDhist(0.22, 22.0);
                std::unique_ptr<TH1F> htemp {(TH1F*)((fitter->Glauber::Fitter::GetNBDHisto()).Clone("htemp"))};
        	float nHits {0.};
        	for (int j=0; j<Na; j++) nHits += (1-alpha*Npart*Npart)*(htemp->GetRandom());
		cout<<"alpha="<<(1-alpha*Npart*Npart)<<"      nhits_MDC="<<nHits<<endl;
        	for (int j=0; j<entries_MDC; j++) {
			Borders_MDC -> GetEntry(j);
			if (nHits>=MinBorder_MDC && nHits<=MaxBorder_MDC) {
				cout<<"true"<<endl;
				B_VS_CentralityHisto_MDC[j]     -> Fill(B);
				Npart_VS_CentralityHisto_MDC[j] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_MDC[j] -> Fill(Ncoll);
				tracksMDC_VS_CentralityHisto[j] -> Fill(nHits);
				B_VS_CentralityHisto_MDC[CentralityClasses]     -> Fill(B);
				Npart_VS_CentralityHisto_MDC[CentralityClasses] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_MDC[CentralityClasses] -> Fill(Ncoll);
				tracksMDC_VS_CentralityHisto[CentralityClasses] -> Fill(nHits);
				break;
				}
			}
		}
	
	tGlauber_FW -> SetBranchAddress("Npart", &Npart);
	tGlauber_FW -> SetBranchAddress("Ncoll", &Ncoll);
	tGlauber_FW -> SetBranchAddress("B", &B);
	fitter -> Glauber::Fitter::SetInputHisto (*((TH1F*)f6->Get("FWSumChargeZ_selected")));
	fitter -> Glauber::Fitter::Init(entries_glauber, "PSD");
	for (int i = 0; i < entries_fw; i++) {
		cout<<i<<endl;		
		tGlauber_FW -> GetEntry(i);
		HADES_FW    -> GetEntry(i);
		
		const int Na = int(NprotonsA);
		fitter->Glauber::Fitter::SetNBDhist(1.4, 20.0);
                std::unique_ptr<TH1F> htemp {(TH1F*)((fitter->Glauber::Fitter::GetNBDHisto()).Clone("htemp"))};
        	float nHits {0.};
        	for (int j=0; j<Na; j++) nHits += (htemp->GetRandom());
		
		cout<<"Na="<<Na<<"   nhits_FW="<<nHits<<endl;
        	for (int j=0; j<entries_FW; j++) {
			Borders_FW -> GetEntry(j);
			if (nHits>=MinBorder_FW && nHits<=MaxBorder_FW) {
				cout<<"true"<<endl;
				B_VS_CentralityHisto_FW[j]     -> Fill(B);
				Npart_VS_CentralityHisto_FW[j] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_FW[j] -> Fill(Ncoll);
				FWSumChargeZ_VS_CentralityHisto[j] -> Fill(nHits);
				B_VS_CentralityHisto_FW[CentralityClasses]     -> Fill(B);
				Npart_VS_CentralityHisto_FW[CentralityClasses] -> Fill(Npart);
				Ncoll_VS_CentralityHisto_FW[CentralityClasses] -> Fill(Ncoll);
				FWSumChargeZ_VS_CentralityHisto[CentralityClasses] -> Fill(nHits);
				break;
				}
			}
		}
	
	double sum;
	Int_t bins;
	for (int i = 0; i < CentralityClasses; i++) {
		cout<<i<<endl;
		Ncc=i+1;
		MinPercent=(i)*100/CentralityClasses;
		MaxPercent=(i+1)*100/CentralityClasses;
		
		Borders_TOF -> GetEntry(i);
		Borders_RPC -> GetEntry(i);
		Borders_MDC -> GetEntry(i);
		Borders_FW  -> GetEntry(i);

		MinBorder_TOF_R=MinBorder_TOF;
		MaxBorder_TOF_R=MaxBorder_TOF;
		MinBorder_RPC_R=MinBorder_RPC;
		MaxBorder_RPC_R=MaxBorder_RPC;
		MinBorder_MDC_R=MinBorder_MDC;
		MaxBorder_MDC_R=MaxBorder_MDC;
		MinBorder_FW_R=MinBorder_FW;
		MaxBorder_FW_R=MaxBorder_FW;

		B_TOF_av=B_VS_CentralityHisto_TOF[i]->GetMean();
		Npart_TOF_av=Npart_VS_CentralityHisto_TOF[i]->GetMean();
		Ncoll_TOF_av=Ncoll_VS_CentralityHisto_TOF[i]->GetMean();
		B_RPC_av=B_VS_CentralityHisto_RPC[i]->GetMean();
		Npart_RPC_av=Npart_VS_CentralityHisto_RPC[i]->GetMean();
		Ncoll_RPC_av=Ncoll_VS_CentralityHisto_RPC[i]->GetMean();
		B_MDC_av=B_VS_CentralityHisto_MDC[i]->GetMean();
		Npart_MDC_av=Npart_VS_CentralityHisto_MDC[i]->GetMean();
		Ncoll_MDC_av=Ncoll_VS_CentralityHisto_MDC[i]->GetMean();
		B_FW_av=B_VS_CentralityHisto_FW[i]->GetMean();
		Npart_FW_av=sum/Npart_VS_CentralityHisto_FW[i]->GetMean();
		Ncoll_FW_av=sum/Ncoll_VS_CentralityHisto_FW[i]->GetMean();
		B_av=(B_TOF_av+B_RPC_av+B_MDC_av+B_FW_av)/4.0;
		Npart_av=(Npart_TOF_av+Npart_RPC_av+Npart_MDC_av+Npart_FW_av)/4.0;
		Ncoll_av=(Ncoll_TOF_av+Ncoll_RPC_av+Ncoll_MDC_av+Ncoll_FW_av)/4.0;

		RMS_B_TOF_av=B_VS_CentralityHisto_TOF[i]->GetRMS();
		RMS_Npart_TOF_av=Npart_VS_CentralityHisto_TOF[i]->GetRMS();
		RMS_Ncoll_TOF_av=Ncoll_VS_CentralityHisto_TOF[i]->GetRMS();
		RMS_B_RPC_av=B_VS_CentralityHisto_RPC[i]->GetRMS();
		RMS_Npart_RPC_av=Npart_VS_CentralityHisto_RPC[i]->GetRMS();
		RMS_Ncoll_RPC_av=Ncoll_VS_CentralityHisto_RPC[i]->GetRMS();
		RMS_B_MDC_av=B_VS_CentralityHisto_MDC[i]->GetRMS();
		RMS_Npart_MDC_av=Npart_VS_CentralityHisto_MDC[i]->GetRMS();
		RMS_Ncoll_MDC_av=Ncoll_VS_CentralityHisto_MDC[i]->GetRMS();
		RMS_B_FW_av=B_VS_CentralityHisto_FW[i]->GetRMS();
		RMS_Npart_FW_av=Npart_VS_CentralityHisto_FW[i]->GetRMS();
		RMS_Ncoll_FW_av=Ncoll_VS_CentralityHisto_FW[i]->GetRMS();
		RMS_B_av=(TMath::Power((TMath::Power(RMS_B_TOF_av,2)+TMath::Power(RMS_B_RPC_av,2)+TMath::Power(RMS_B_MDC_av,2)+TMath::Power(RMS_B_FW_av,2)),0.5))/4.0;
		RMS_B_av=(TMath::Power((TMath::Power(RMS_B_TOF_av,2)+TMath::Power(RMS_B_RPC_av,2)+TMath::Power(RMS_B_MDC_av,2)),0.5))/3.0;
		RMS_Npart_av=(TMath::Power((TMath::Power(RMS_Npart_TOF_av,2)+TMath::Power(RMS_Npart_RPC_av,2)+TMath::Power(RMS_Npart_MDC_av,2)+TMath::Power(RMS_Npart_FW_av,2)),0.5))/4.0;
		RMS_Npart_av=(TMath::Power((TMath::Power(RMS_Npart_TOF_av,2)+TMath::Power(RMS_Npart_RPC_av,2)+TMath::Power(RMS_Npart_MDC_av,2)),0.5))/3.0;
		RMS_Ncoll_av=(TMath::Power((TMath::Power(RMS_Ncoll_TOF_av,2)+TMath::Power(RMS_Ncoll_RPC_av,2)+TMath::Power(RMS_Ncoll_MDC_av,2)+TMath::Power(RMS_Ncoll_FW_av,2)),0.5))/4.0;
		RMS_Ncoll_av=(TMath::Power((TMath::Power(RMS_Ncoll_TOF_av,2)+TMath::Power(RMS_Ncoll_RPC_av,2)+TMath::Power(RMS_Ncoll_MDC_av,2)),0.5))/3.0;

		B_average_VS_Centrality         -> SetBinContent(i+1,B_av);
		Npart_average_VS_Centrality     -> SetBinContent(i+1,Npart_av);
		Ncoll_average_VS_Centrality     -> SetBinContent(i+1,Ncoll_av);	
		B_average_VS_Centrality_TOF     -> SetBinContent(i+1,B_TOF_av);
		Npart_average_VS_Centrality_TOF -> SetBinContent(i+1,Npart_TOF_av);
		Ncoll_average_VS_Centrality_TOF -> SetBinContent(i+1,Ncoll_TOF_av);
		B_average_VS_Centrality_RPC     -> SetBinContent(i+1,B_RPC_av);
		Npart_average_VS_Centrality_RPC -> SetBinContent(i+1,Npart_RPC_av);
		Ncoll_average_VS_Centrality_RPC -> SetBinContent(i+1,Ncoll_RPC_av);
		B_average_VS_Centrality_MDC     -> SetBinContent(i+1,B_MDC_av);
		Npart_average_VS_Centrality_MDC -> SetBinContent(i+1,Npart_MDC_av);
		Ncoll_average_VS_Centrality_MDC -> SetBinContent(i+1,Ncoll_MDC_av);
		B_average_VS_Centrality_FW      -> SetBinContent(i+1,B_FW_av);
		Npart_average_VS_Centrality_FW  -> SetBinContent(i+1,Npart_FW_av);
		Ncoll_average_VS_Centrality_FW  -> SetBinContent(i+1,Ncoll_FW_av);
	
		B_average_VS_Centrality         -> SetBinError(i+1,RMS_B_av);
		Npart_average_VS_Centrality     -> SetBinError(i+1,RMS_Npart_av);
		Ncoll_average_VS_Centrality     -> SetBinError(i+1,RMS_Ncoll_av);	
		B_average_VS_Centrality_TOF     -> SetBinError(i+1,RMS_B_TOF_av);
		Npart_average_VS_Centrality_TOF -> SetBinError(i+1,RMS_Npart_TOF_av);
		Ncoll_average_VS_Centrality_TOF -> SetBinError(i+1,RMS_Ncoll_TOF_av);
		B_average_VS_Centrality_RPC     -> SetBinError(i+1,RMS_B_RPC_av);
		Npart_average_VS_Centrality_RPC -> SetBinError(i+1,RMS_Npart_RPC_av);
		Ncoll_average_VS_Centrality_RPC -> SetBinError(i+1,RMS_Ncoll_RPC_av);
		B_average_VS_Centrality_MDC     -> SetBinError(i+1,RMS_B_MDC_av);
		Npart_average_VS_Centrality_MDC -> SetBinError(i+1,RMS_Npart_MDC_av);
		Ncoll_average_VS_Centrality_MDC -> SetBinError(i+1,RMS_Ncoll_MDC_av);
		B_average_VS_Centrality_FW      -> SetBinError(i+1,RMS_B_FW_av);
		Npart_average_VS_Centrality_FW  -> SetBinError(i+1,RMS_Npart_FW_av);
		Ncoll_average_VS_Centrality_FW  -> SetBinError(i+1,RMS_Ncoll_FW_av);

		Result->Fill();
		}

	TCanvas* c1 = new TCanvas("B_VS_Centrality_TOF","B_VS_Centrality_TOF");
	TLegend legend1(0.1, 0.2, 0.3, 0.4, "B_VS_Centrality_TOF");
	for (int i = 0; i < CentralityClasses; i++) legend1.AddEntry(B_VS_CentralityHisto_TOF[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend1.AddEntry(B_VS_CentralityHisto_TOF[CentralityClasses], "CentralityClass 0%-100%");
	THStack *B_VS_Centrality_TOF = new THStack("B_VS_Centrality_TOF","B_VS_Centrality_TOF");
	for (int i = 0; i <= CentralityClasses; i++) B_VS_Centrality_TOF -> Add(B_VS_CentralityHisto_TOF[i]);
	B_VS_Centrality_TOF -> Draw("nostack");
	legend1.Draw("same");
	c1->Write();
	for (int i = 0; i <= CentralityClasses; i++) B_VS_CentralityHisto_TOF[i] -> Write();

	TCanvas* c2 = new TCanvas("Npart_VS_Centrality_TOF","Npart_VS_Centrality_TOF");
	TLegend legend2(0.1, 0.2, 0.3, 0.4, "Npart_VS_Centrality_TOF");
	for (int i = 0; i < CentralityClasses; i++) legend2.AddEntry(Npart_VS_CentralityHisto_TOF[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend2.AddEntry(Npart_VS_CentralityHisto_TOF[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Npart_VS_Centrality_TOF = new THStack("Npart_VS_Centrality_TOF","Npart_VS_Centrality_TOF");
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_Centrality_TOF -> Add(Npart_VS_CentralityHisto_TOF[i]);
	Npart_VS_Centrality_TOF -> Draw("nostack");
	legend2.Draw("same");
	c2->Write();
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_CentralityHisto_TOF[i] -> Write();

	TCanvas* c3 = new TCanvas("Ncoll_VS_Centrality_TOF","Ncoll_VS_Centrality_TOF");
	TLegend legend3(0.1, 0.2, 0.3, 0.4, "Ncoll_VS_Centrality_TOF");
	for (int i = 0; i < CentralityClasses; i++) legend3.AddEntry(Ncoll_VS_CentralityHisto_TOF[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend3.AddEntry(Ncoll_VS_CentralityHisto_TOF[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Ncoll_VS_Centrality_TOF = new THStack("Ncoll_VS_Centrality_TOF","Ncoll_VS_Centrality_TOF");
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_Centrality_TOF -> Add(Ncoll_VS_CentralityHisto_TOF[i]);
	Ncoll_VS_Centrality_TOF -> Draw("nostack");
	legend3.Draw("same");
	c3->Write();
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_CentralityHisto_TOF[i] -> Write();

	TCanvas* c4 = new TCanvas("B_VS_Centrality_RPC","B_VS_Centrality_RPC");
	TLegend legend4(0.1, 0.2, 0.3, 0.4, "B_VS_Centrality_RPC");
	for (int i = 0; i < CentralityClasses; i++) legend4.AddEntry(B_VS_CentralityHisto_RPC[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend4.AddEntry(B_VS_CentralityHisto_RPC[CentralityClasses], "CentralityClass 0%-100%");
	THStack *B_VS_Centrality_RPC = new THStack("B_VS_Centrality_RPC","B_VS_Centrality_RPC");
	for (int i = 0; i <= CentralityClasses; i++) B_VS_Centrality_RPC -> Add(B_VS_CentralityHisto_RPC[i]);
	B_VS_Centrality_RPC -> Draw("nostack");
	legend4.Draw("same");
	c4->Write();
	for (int i = 0; i <= CentralityClasses; i++) B_VS_CentralityHisto_RPC[i] -> Write();

	TCanvas* c5 = new TCanvas("Npart_VS_Centrality_RPC","Npart_VS_Centrality_RPC");
	TLegend legend5(0.1, 0.2, 0.3, 0.4, "Npart_VS_Centrality_RPC");
	for (int i = 0; i < CentralityClasses; i++) legend5.AddEntry(Npart_VS_CentralityHisto_RPC[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend5.AddEntry(Npart_VS_CentralityHisto_RPC[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Npart_VS_Centrality_RPC = new THStack("Npart_VS_Centrality_RPC","Npart_VS_Centrality_RPC");
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_Centrality_RPC -> Add(Npart_VS_CentralityHisto_RPC[i]);
	Npart_VS_Centrality_RPC -> Draw("nostack");
	legend5.Draw("same");
	c5->Write();
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_CentralityHisto_RPC[i] -> Write();

	TCanvas* c6 = new TCanvas("Ncoll_VS_Centrality_RPC","Ncoll_VS_Centrality_RPC");
	TLegend legend6(0.1, 0.2, 0.3, 0.4, "Ncoll_VS_Centrality_RPC");
	for (int i = 0; i < CentralityClasses; i++) legend6.AddEntry(Ncoll_VS_CentralityHisto_RPC[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend6.AddEntry(Ncoll_VS_CentralityHisto_RPC[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Ncoll_VS_Centrality_RPC = new THStack("Ncoll_VS_Centrality_RPC","Ncoll_VS_Centrality_RPC");
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_Centrality_RPC -> Add(Ncoll_VS_CentralityHisto_RPC[i]);
	Ncoll_VS_Centrality_RPC -> Draw("nostack");
	legend6.Draw("same");
	c6->Write();
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_CentralityHisto_RPC[i] -> Write();

	TCanvas* c7 = new TCanvas("B_VS_Centrality_MDC","B_VS_Centrality_MDC");
	TLegend legend7(0.1, 0.2, 0.3, 0.4, "B_VS_Centrality_MDC");
	for (int i = 0; i < CentralityClasses; i++) legend7.AddEntry(B_VS_CentralityHisto_MDC[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend7.AddEntry(B_VS_CentralityHisto_MDC[CentralityClasses], "CentralityClass 0%-100%");
	THStack *B_VS_Centrality_MDC = new THStack("B_VS_Centrality_MDC","B_VS_Centrality_MDC");
	for (int i = 0; i <= CentralityClasses; i++) B_VS_Centrality_MDC -> Add(B_VS_CentralityHisto_MDC[i]);
	B_VS_Centrality_MDC -> Draw("nostack");
	legend7.Draw("same");
	c7->Write();
	for (int i = 0; i <= CentralityClasses; i++) B_VS_CentralityHisto_MDC[i] -> Write();

	TCanvas* c8 = new TCanvas("Npart_VS_Centrality_MDC","Npart_VS_Centrality_MDC");
	TLegend legend8(0.1, 0.2, 0.3, 0.4, "Npart_VS_Centrality_MDC");
	for (int i = 0; i < CentralityClasses; i++) legend8.AddEntry(Npart_VS_CentralityHisto_MDC[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend8.AddEntry(Npart_VS_CentralityHisto_MDC[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Npart_VS_Centrality_MDC = new THStack("Npart_VS_Centrality_MDC","Npart_VS_Centrality_MDC");
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_Centrality_MDC -> Add(Npart_VS_CentralityHisto_MDC[i]);
	Npart_VS_Centrality_MDC -> Draw("nostack");
	legend8.Draw("same");
	c8->Write();
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_CentralityHisto_MDC[i] -> Write();

	TCanvas* c9 = new TCanvas("Ncoll_VS_Centrality_MDC","Ncoll_VS_Centrality_MDC");
	TLegend legend9(0.1, 0.2, 0.3, 0.4, "Ncoll_VS_Centrality_MDC");
	for (int i = 0; i < CentralityClasses; i++) legend9.AddEntry(Ncoll_VS_CentralityHisto_MDC[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend9.AddEntry(Ncoll_VS_CentralityHisto_MDC[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Ncoll_VS_Centrality_MDC = new THStack("Ncoll_VS_Centrality_MDC","Ncoll_VS_Centrality_MDC");
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_Centrality_MDC -> Add(Ncoll_VS_CentralityHisto_MDC[i]);
	Ncoll_VS_Centrality_MDC -> Draw("nostack");
	legend9.Draw("same");
	c9->Write();
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_CentralityHisto_MDC[i] -> Write();

	TCanvas* c10 = new TCanvas("B_VS_Centrality_FW","B_VS_Centrality_FW");
	TLegend legend10(0.1, 0.2, 0.3, 0.4, "B_VS_Centrality_FW");
	for (int i = 0; i < CentralityClasses; i++) legend10.AddEntry(B_VS_CentralityHisto_FW[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend10.AddEntry(B_VS_CentralityHisto_FW[CentralityClasses], "CentralityClass 0%-100%");
	THStack *B_VS_Centrality_FW = new THStack("B_VS_Centrality_FW","B_VS_Centrality_FW");
	for (int i = 0; i <= CentralityClasses; i++) B_VS_Centrality_FW -> Add(B_VS_CentralityHisto_FW[i]);
	B_VS_Centrality_FW -> Draw("nostack");
	legend10.Draw("same");
	c10->Write();
	for (int i = 0; i <= CentralityClasses; i++) B_VS_CentralityHisto_FW[i] -> Write();

	TCanvas* c11 = new TCanvas("Npart_VS_Centrality_FW","Npart_VS_Centrality_FW");
	TLegend legend11(0.1, 0.2, 0.3, 0.4, "Npart_VS_Centrality_FW");
	for (int i = 0; i < CentralityClasses; i++) legend11.AddEntry(Npart_VS_CentralityHisto_FW[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend11.AddEntry(Npart_VS_CentralityHisto_FW[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Npart_VS_Centrality_FW = new THStack("Npart_VS_Centrality_FW","Npart_VS_Centrality_FW");
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_Centrality_FW -> Add(Npart_VS_CentralityHisto_FW[i]);
	Npart_VS_Centrality_FW -> Draw("nostack");
	legend11.Draw("same");
	c11->Write();
	for (int i = 0; i <= CentralityClasses; i++) Npart_VS_CentralityHisto_FW[i] -> Write();

	TCanvas* c12 = new TCanvas("Ncoll_VS_Centrality_FW","Ncoll_VS_Centrality_FW");
	TLegend legend12(0.1, 0.2, 0.3, 0.4, "Ncoll_VS_Centrality_FW");
	for (int i = 0; i < CentralityClasses; i++) legend12.AddEntry(Ncoll_VS_CentralityHisto_FW[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	legend12.AddEntry(Ncoll_VS_CentralityHisto_FW[CentralityClasses], "CentralityClass 0%-100%");
	THStack *Ncoll_VS_Centrality_FW = new THStack("Ncoll_VS_Centrality_FW","Ncoll_VS_Centrality_FW");
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_Centrality_FW -> Add(Ncoll_VS_CentralityHisto_FW[i]);
	Ncoll_VS_Centrality_FW -> Draw("nostack");
	legend12.Draw("same");
	c12->Write();
	for (int i = 0; i <= CentralityClasses; i++) Ncoll_VS_CentralityHisto_FW[i] -> Write();

	TCanvas* c13 = new TCanvas("B_average_VS_Centrality","B_average_VS_Centrality");
	TLegend legend13(0.1, 0.2, 0.3, 0.4, "B_average_VS_Centrality");
	legend13.AddEntry(B_average_VS_Centrality_TOF, "TOF");
	legend13.AddEntry(B_average_VS_Centrality_RPC, "RPC");
	legend13.AddEntry(B_average_VS_Centrality_MDC, "MDC");
	legend13.AddEntry(B_average_VS_Centrality_FW, "FW");
	B_average_VS_Centrality_TOF -> Draw("E1, P");
	B_average_VS_Centrality_RPC -> Draw("same, E1, P");
	B_average_VS_Centrality_MDC -> Draw("same, E1, P");
	B_average_VS_Centrality_FW  -> Draw("same, E1, P");
	legend13.Draw("same");
	c13->Write();

	TCanvas* c14 = new TCanvas("Npart_average_VS_Centrality","Npart_average_VS_Centrality");
	TLegend legend14(0.1, 0.2, 0.3, 0.4, "Npart_average_VS_Centrality");
	legend14.AddEntry(Npart_average_VS_Centrality_TOF, "TOF");
	legend14.AddEntry(Npart_average_VS_Centrality_RPC, "RPC");
	legend14.AddEntry(Npart_average_VS_Centrality_MDC, "MDC");
	legend14.AddEntry(Npart_average_VS_Centrality_FW, "FW");
	Npart_average_VS_Centrality_TOF -> Draw("E1, P");
	Npart_average_VS_Centrality_RPC -> Draw("same, E1, P");
	Npart_average_VS_Centrality_MDC -> Draw("same, E1, P");
	Npart_average_VS_Centrality_FW  -> Draw("same, E1, P");
	legend14.Draw("same");
	c14->Write();

	TCanvas* c15 = new TCanvas("Ncoll_average_VS_Centrality","Ncoll_average_VS_Centrality");
	TLegend legend15(0.1, 0.2, 0.3, 0.4, "Ncoll_average_VS_Centrality");
	legend15.AddEntry(Ncoll_average_VS_Centrality_TOF, "TOF");
	legend15.AddEntry(Ncoll_average_VS_Centrality_RPC, "RPC");
	legend15.AddEntry(Ncoll_average_VS_Centrality_MDC, "MDC");
	legend15.AddEntry(Ncoll_average_VS_Centrality_FW, "FW");
	Ncoll_average_VS_Centrality_TOF -> Draw("E1, P");
	Ncoll_average_VS_Centrality_RPC -> Draw("same, E1, P");
	Ncoll_average_VS_Centrality_MDC -> Draw("same, E1, P");
	Ncoll_average_VS_Centrality_FW  -> Draw("same, E1, P");
	legend15.Draw("same");
	c15->Write();

	TCanvas* c16 = new TCanvas("hitsTOF_VS_Centrality","hitsTOF_VS_Centrality");
	TLegend legend16(0.1, 0.2, 0.3, 0.4, "hitsTOF_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) legend16.AddEntry(hitsTOF_VS_CentralityHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
//	legend16.AddEntry(hitsTOF_VS_CentralityHisto[CentralityClasses], "CentralityClass 0%-100%");
	THStack *hitsTOF_VS_Centrality = new THStack("hitsTOF_VS_Centrality","hitsTOF_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) hitsTOF_VS_Centrality -> Add(hitsTOF_VS_CentralityHisto[i]);
	hitsTOF_VS_Centrality -> Draw("nostack");
	legend16.Draw("same");
	c16->Write();
	for (int i = 0; i <= CentralityClasses; i++) hitsTOF_VS_CentralityHisto[i] -> Write();

	TCanvas* c17 = new TCanvas("hitsRPC_VS_Centrality","hitsRPC_VS_Centrality");
	TLegend legend17(0.1, 0.2, 0.3, 0.4, "hitsRPC_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) legend17.AddEntry(hitsRPC_VS_CentralityHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
//	legend17.AddEntry(hitsRPC_VS_CentralityHisto[CentralityClasses], "CentralityClass 0%-100%");
	THStack *hitsRPC_VS_Centrality = new THStack("hitsRPC_VS_Centrality","hitsRPC_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) hitsRPC_VS_Centrality -> Add(hitsRPC_VS_CentralityHisto[i]);
	hitsRPC_VS_Centrality -> Draw("nostack");
	legend17.Draw("same");
	c17->Write();
	for (int i = 0; i <= CentralityClasses; i++) hitsRPC_VS_CentralityHisto[i] -> Write();

	TCanvas* c18 = new TCanvas("tracksMDC_VS_Centrality","tracksMDC_VS_Centrality");
	TLegend legend18(0.1, 0.2, 0.3, 0.4, "tracksMDC_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) legend18.AddEntry(tracksMDC_VS_CentralityHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
//	legend18.AddEntry(tracksMDC_VS_CentralityHisto[CentralityClasses], "CentralityClass 0%-100%");
	THStack *tracksMDC_VS_Centrality = new THStack("tracksMDC_VS_Centrality","tracksMDC_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) tracksMDC_VS_Centrality -> Add(tracksMDC_VS_CentralityHisto[i]);
	tracksMDC_VS_Centrality -> Draw("nostack");
	legend18.Draw("same");
	c18->Write();
	for (int i = 0; i <= CentralityClasses; i++) tracksMDC_VS_CentralityHisto[i] -> Write();

	TCanvas* c19 = new TCanvas("FWSumChargeZ_VS_Centrality","FWSumChargeZ_VS_Centrality");
	TLegend legend19(0.1, 0.2, 0.3, 0.4, "FWSumChargeZ_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) legend19.AddEntry(FWSumChargeZ_VS_CentralityHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
//	legend19.AddEntry(FWSumChargeZ_VS_CentralityHisto[CentralityClasses], "CentralityClass 0%-100%");
	THStack *FWSumChargeZ_VS_Centrality = new THStack("FWSumChargeZ_VS_Centrality","FWSumChargeZ_VS_Centrality");
	for (int i = 0; i < CentralityClasses; i++) FWSumChargeZ_VS_Centrality -> Add(FWSumChargeZ_VS_CentralityHisto[i]);
	FWSumChargeZ_VS_Centrality -> Draw("nostack");
	legend19.Draw("same");
	c19->Write();
	for (int i = 0; i <= CentralityClasses; i++) FWSumChargeZ_VS_CentralityHisto[i] -> Write();

	B_average_VS_Centrality         -> Write();
	Npart_average_VS_Centrality     -> Write();
	Ncoll_average_VS_Centrality     -> Write();	
	B_average_VS_Centrality_TOF     -> Write();
	Npart_average_VS_Centrality_TOF -> Write();
	Ncoll_average_VS_Centrality_TOF -> Write();
	B_average_VS_Centrality_RPC     -> Write();
	Npart_average_VS_Centrality_RPC -> Write();
	Ncoll_average_VS_Centrality_RPC -> Write();
	B_average_VS_Centrality_MDC     -> Write();
	Npart_average_VS_Centrality_MDC -> Write();
	Ncoll_average_VS_Centrality_MDC -> Write();
	B_average_VS_Centrality_FW      -> Write();
	Npart_average_VS_Centrality_FW  -> Write();
	Ncoll_average_VS_Centrality_FW  -> Write();
	
	Result->Write();

	f->Close();
}
