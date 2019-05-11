#include <iostream>

#include "Fitter.h"
#include "FitterHelper.h"

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLegend.h"
#include "TH2.h"

using namespace Glauber;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cout << "Wrong number of parameters! Executable usage:" << std::endl;
        std::cout << "   ./glauber f0 k0" << std::endl;
        return -1;
    }
    const Float_t f0 = atof( argv[1]);
    std::cout << "f=" << f0 << std::endl;
    const Float_t f1 = atof( argv[2]);
    std::cout << "f=" << f1 << std::endl;
    const Int_t k0 = atoi( argv[3] );
    std::cout << "k0=" << k0 << std::endl;
    const Int_t k1 = atoi( argv[4] );
    std::cout << "k1=" << k1 << std::endl;
    
    // *****************************************
    // Modify this part according to your needs
    // *****************************************

    ///  |   mode    |   function for Na      |
    ///  |   Default | f*Npart + (1-f)*Ncoll  |
    ///  |   PSD     | Nprotons in spectators |
    ///  |   Npart   |     Npart^f            |
    ///  |   Ncoll   |     Ncoll^f            |
    const TString mode = "Default";
//    const TString mode = "PSD";
    
    const TString glauber_filename = "/home/vad/NIR_codes/Centrality/gmc-Au3Au3-snn23.6-md0.4-nd-1.0-rc1-smax99.0_1000000.root";   // input files
//    const TString glauber_filename = "/home/vad/NIR_codes/Centrality/2.root";
    TString glauber_treename;
    glauber_treename = "nt_Au3_Au3";
    if (mode == "PSD") glauber_treename = "HADES_FW";
    const TString in_filename = "/home/vad/NIR_codes/QAHistoBuilder/QASelectedPT2Histos.root";

//    const TString histoname = "tracksMDC_selected";
//    const TString histoname = "hitsTOF_selected";
    const TString histoname = "hitsRPC_selected";
//    const TString histoname = "hitsTOF+RPC_selected";
//    const TString histoname = "FWSumChargeZ_selected";


    const Int_t min_bin = 100;      // not fitting low multiplicity region due to trigger bias, etc
    const Int_t max_bin = 188;   // very large number to fit the whole histo
    double alpha = 1.64e-6;
    

    const TString outdir = ".";
    // *****************************************
    // *****************************************

    std::unique_ptr<TFile> glauber_file{ TFile::Open(glauber_filename, "read") };
    std::unique_ptr<TTree> glauber_tree{ (TTree*) glauber_file->Get(glauber_treename) };
    
    std::unique_ptr<TFile> f{TFile::Open(in_filename)};    
    TH1F *hdata = (TH1F*)f->Get(histoname);
    const Int_t nevents = 10*(int(hdata->Integral(min_bin,max_bin)));
//    const Int_t nevents = 99999;

    Fitter fitter ( std::move(glauber_tree), mode );

    fitter.SetMode(mode);
    fitter.SetMassNumber(f0/2);
    fitter.SetInputHisto(*hdata);
    fitter.SetBinSize(1);
    fitter.Init(nevents, mode);
    
    fitter.SetFitMinBin(min_bin);
    fitter.SetFitMaxBin(max_bin);
    fitter.SetOutDirName(outdir);

    float par[4];
//    const float chi2 = fitter.FitGlauber(par, f0, f1, k0, k1, nevents);
	
    float chi2=1e10;
    chi2 = fitter.FitGlauber(par, f0, f1, k0, k1, nevents, alpha);

    std::cout << "f = " << par[0] << "    mu = " << par[1] << "    k = " << par[2] << "    chi2 = " << chi2 << "    chi2_error = " << par[3] << std::endl; 
    
//    Glauber::DrawHistos(fitter, true, true, true, true);

    DrawHistos(fitter, true, true, true, true);

    const float range[2] = {300, 350.};
    std::unique_ptr<TH1F> hB(fitter.GetModelHisto (range, "B", par, 100000));
    hB->SaveAs( "b_test.root" );
    
    std::cout << "END!" << std::endl;

    return 0;
}
