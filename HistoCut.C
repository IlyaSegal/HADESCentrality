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

void HistoCut(Int_t CentralityClasses, Int_t CE){
	TFile *f1 = new TFile("/home/vad/centrality/build/Result/glauber_qa_RPC.root");
	TFile *f2 = new TFile("/home/vad/NIR_codes/QAHistoBuilder/QASelectedPT2Histos.root");
	TH1F *FitHisto=(TH1F*)f1->Get("glaub_fit_histo");
	TH1F *DataHisto=(TH1F*)f2->Get(Form("%s_selected", CentralityEstimatorName[CE].Data()));

	Int_t bins = FitHisto->GetNbinsX();
	cout<<"bins="<<bins<<endl;
	Double_t max = FitHisto->GetXaxis()->GetBinCenter(bins);
	cout<<"max="<<max<<endl;
	Double_t min = FitHisto->GetYaxis()->GetXmin();
	cout<<"min="<<min<<endl;
	Double_t integral = FitHisto->Integral(2,bins);
	cout<<"integral="<<integral<<endl;

	TFile *f = new TFile("HistoCutResult.root", "recreate");
	TH1F* ResultHisto[CentralityClasses];
	for (int i = 0; i < CentralityClasses; i++)   ResultHisto[i] = new TH1F(Form("CentralityClass_Fit %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), Form(";%s;counts", CentralityEstimatorName[CE].Data()),bins, min, max);
	Int_t j=1;
	for (int i = 0; i < CentralityClasses; i++) {
		ResultHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses; i++) ResultHisto[i] -> SetLineWidth(4);

	TH1F* CentralityHisto[CentralityClasses];
	for (int i = 0; i < CentralityClasses; i++)   CentralityHisto[i] = new TH1F(Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), Form(";%s;counts", CentralityEstimatorName[CE].Data()), bins, min, max);
	j=1; 
	for (int i = 0; i < CentralityClasses; i++) {
		CentralityHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses; i++) CentralityHisto[i] -> SetLineWidth(4);

	TH1F* Centrality_vs_Multiplisity = new TH1F("Centrality_vs_Multiplisity", Form(";%s;CentralityPercent", CentralityEstimatorName[CE].Data()),bins, min, max);
	Centrality_vs_Multiplisity -> SetLineColor(1);
	Centrality_vs_Multiplisity -> SetLineWidth(4);
	
	j=CentralityClasses-1;
	Double_t sum=0.0;
	for (Int_t i=2;i<=bins;i++) {
        	sum = sum+FitHisto->GetBinContent(i);
		cout<<sum<<endl;
		cout<<j<<endl;
      		if (sum < (integral/CentralityClasses) && j<CentralityClasses) ResultHisto[j] -> SetBinContent(i, FitHisto->GetBinContent(i));
		else {
			if (j>0) {j--; sum=0.0;}
			if (j>0)ResultHisto[j] -> SetBinContent(i, FitHisto->GetBinContent(i));
			}
   		}

	j=CentralityClasses-1;
	for (Int_t i=2;i<=bins;i++) {
        	if ((ResultHisto[j]->GetBinContent(i))>0) {
			CentralityHisto[j] -> SetBinContent(i, DataHisto->GetBinContent(i));
			Centrality_vs_Multiplisity -> SetBinContent(i, (j+1)*100/CentralityClasses);
			}
		else {
			if (j>0) j--;
			if (j>0) {
				CentralityHisto[j]         -> SetBinContent(i, DataHisto->GetBinContent(i));
				Centrality_vs_Multiplisity -> SetBinContent(i, (j+1)*100/CentralityClasses);
				}
			}
   		}

	

	TCanvas* c1 = new TCanvas("HistoCutResult","HistoCutResult");
	TLegend legend1(0.1, 0.2, 0.3, 0.4, "HistoCutResult");
	for (int i = 0; i < CentralityClasses; i++) legend1.AddEntry(ResultHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	THStack *ResultHistos = new THStack("HistoCutResult","HistoCutResult");
	for (int i = 0; i < CentralityClasses; i++) ResultHistos -> Add(ResultHisto[i]);
	ResultHistos -> Draw("nostack");
	legend1.Draw("same");
	c1->Write();
	for (int i = 0; i < CentralityClasses; i++) ResultHisto[i] -> Write();

	TCanvas* c2 = new TCanvas("CentraliryClasses","CentralityClasses");
	TLegend legend2(0.1, 0.2, 0.3, 0.4, "CentraliryClasses");
	for (int i = 0; i < CentralityClasses; i++) legend2.AddEntry(CentralityHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	THStack *CentralityHistos = new THStack("HistoCutResult","HistoCutResult");
	for (int i = 0; i < CentralityClasses; i++) CentralityHistos -> Add(CentralityHisto[i]);
	CentralityHistos -> Draw("nostack");
	legend2.Draw("same");
	c2->Write();
	for (int i = 0; i < CentralityClasses; i++) CentralityHisto[i] -> Write();
	Centrality_vs_Multiplisity -> Write();
	

	f->Close();
}
