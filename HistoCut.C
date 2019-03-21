void HistoCut(Int_t CentralityClasses, TH1F* Histo){
	
	Int_t bins = Histo->GetNbinsX();
	cout<<"bins="<<bins<<endl;
	Double_t max = Histo->GetXaxis()->GetBinCenter(bins);
	cout<<"max="<<max<<endl;
	Double_t min = Histo->GetYaxis()->GetXmin();
	cout<<"min="<<min<<endl;
	Double_t integral = Histo->Integral(2,bins);
	cout<<"integral="<<integral<<endl;

	TFile *f = new TFile("HistoCutResult.root", "recreate");
	TH1F* ResultHisto[CentralityClasses];
	for (int i = 0; i < CentralityClasses; i++)   ResultHisto[i] = new TH1F(Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses), ";centralityEstimator;counts",bins,min,max);
	Int_t j=1;
	for (int i = 0; i < CentralityClasses; i++) {
		ResultHisto[i] -> SetLineColor(j);
		j++;
		if (j==9) j=1;
		}
	for (int i = 0; i < CentralityClasses; i++) ResultHisto[i] -> SetLineWidth(4);
	
	j=0;
	Double_t sum=0.0;

	for (Int_t i=2;i<=bins;i++) {
        	sum = sum+Histo->GetBinContent(i);
		cout<<sum<<endl;
		cout<<j<<endl;
      		if (sum < (integral/CentralityClasses) && j<CentralityClasses) ResultHisto[j] -> SetBinContent(i, Histo->GetBinContent(i));
		else {
			if (j<CentralityClasses-1) {j++; sum=0.0;}
			if (j<CentralityClasses)ResultHisto[j] -> SetBinContent(i, Histo->GetBinContent(i));
			}
   		}

	

	TCanvas* c1 = new TCanvas("HistoCutResult","HistoCutResult");
	TLegend legend(0.1, 0.2, 0.3, 0.4, "HistoCutResult");
	for (int i = 0; i < CentralityClasses; i++) legend.AddEntry(ResultHisto[i], Form("CentralityClass %.1f%%-%.1f%%", i*100.0/CentralityClasses, (i+1)*100.0/CentralityClasses));
	THStack *ResultHistos = new THStack("HistoCutResult","HistoCutResult");
	for (int i = 0; i < CentralityClasses; i++) ResultHistos -> Add(ResultHisto[i]);
	ResultHistos -> Draw("nostack");
	legend.Draw("same");
	c1->Write();
	for (int i = 0; i < CentralityClasses; i++) ResultHisto[i] -> Write();
	

	f->Close();
}
