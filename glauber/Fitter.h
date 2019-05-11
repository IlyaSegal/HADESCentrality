/** @file   Fitter.h
    @class  Glauber::Fitter
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to fit histo with Glauber based function
*/

#ifndef GlauberFitter_H
#define GlauberFitter_H 1

#include <vector>
#include "TString.h"
#include "TNamed.h"
#include "TH1F.h"
#include "TTree.h"
// #include "TMinuit.h"


namespace Glauber
{
    class Fitter
    {
        
    public:
        
        /**   Default constructor   **/
        Fitter() {};
        Fitter(std::unique_ptr<TTree> tree, TString fmode) ;
        /**   Destructor   **/
        virtual ~Fitter(){};
        
        void Init(int nEntries, TString fmode);
        void SetGlauberFitHisto (Float_t f, Float_t mu, Float_t k, Int_t n = 10000, double alpha = 0, Bool_t Norm2Data = true);
//	void SetGlauberFitHistoFW (Float_t f, Float_t mu, Float_t k, Int_t n = 10000, Float_t alpha = 0, Bool_t Norm2Data = true);
        void NormalizeGlauberFit ();
        void DrawHistos (Bool_t isSim = true, Bool_t isData = true, Bool_t isGlauber = false, Bool_t isNBD = false);
        
        float FitGlauber (float *par, Float_t f0, Float_t f1, Int_t k0, Int_t k1, Int_t nEvents, double alpha=0);
        void FindMuGoldenSection (Float_t *mu, Float_t *chi2, float*chi2_error, Float_t mu_min, Float_t mu_max, Float_t f, Float_t k, Int_t nEvents = 10000, Int_t nIter = 5, double alpha=0, int n=0);
	
//        float FitFW (float *par, Float_t f0, Float_t f1, Int_t k0, Int_t k1, Int_t nEvents);
//        void FindMuFWGoldenSection (Float_t *mu, Float_t *chi2, Float_t mu_min, Float_t mu_max, Float_t f, Float_t k, Int_t nEvents = 10000, Int_t nIter = 5, int n=0);
        
        Float_t GetChi2 (void) const;
	Float_t GetChi2Error (void) const;
        
        Float_t NBD(Float_t n, Float_t mu, Float_t k) const;
        void SetNBDhist(Float_t mu, Float_t k);
//	Double_t FW(Float_t n, Float_t mu, Float_t k) const;
//      void SetFWhist(Float_t mu, Float_t k, Float_t f);

        float Nancestors(float f) const;
        float NancestorsMax(float f) const;
        
        std::unique_ptr<TH1F> GetModelHisto (const Float_t range[2], TString name, const Float_t par[4], Int_t nEvents);
        
//         
//         Setters
//         
        void SetInputHisto (const TH1F &h)   { fDataHisto = h; }
        void SetFitMinBin  (Int_t min)      { fFitMinBin = min; }
        void SetFitMaxBin  (Int_t min)      { fFitMaxBin = min; }
        void SetNormMinBin  (Int_t min)     { fNormMinBin = min; }
        void SetBinSize  (Int_t size)        { fBinSize = size; }
        void SetOutDirName (TString name)    { fOutDirName = name; }
        void SetMode (const TString mode) { fMode = mode; }
	void SetMassNumber (Float_t A) { fA = A; }
 
//         
//         Getters
//         
        TH1F GetGlauberFitHisto () const { return fGlauberFitHisto; }
        TH1F GetDataHisto ()       const { return fDataHisto;  }
        TH1F GetNBDHisto ()        const { return fNbdHisto;   }
//	TH1F GetFWHisto ()         const { return fFWHisto;   }
        TH1F GetNpartHisto ()      const { return fNpartHisto; }
	TH1F GetNpartAHisto ()     const { return fNpartAHisto; }
        TH1F GetNcollHisto ()      const { return fNcollHisto;  }
	TH1F GetNprotonsAHisto ()  const { return fNprotonsAHisto;  }
        TH1F GetNspecAHisto ()     const { return fNspecAHisto;  }
        TH1F GetBestFiHisto ()     const { return fBestFitHisto;  }

        
    private:
        
        /**   Data members  **/
        TH1F fNpartHisto;
        TH1F fNprotonsAHisto;
        TH1F fNspecAHisto;
	TH1F fNpartAHisto; 
        TH1F fNcollHisto; 
        TH1F fDataHisto; 
        TH1F fNbdHisto;
//	TH1F fFWHisto; 
        TH1F fGlauberFitHisto; 
        TH1F fBestFitHisto;
        
        /* MC data */
        std::unique_ptr<TTree> fSimTree{nullptr};
        
	Float_t fA{-1.}; //mass number
        Float_t fNpart{-1.};
	Float_t fNpartA{-1.};
        Float_t fNcoll{-1.};
        Float_t fNprotonsA{-1.};
	Float_t fNspecA{-1.};

        Float_t fMaxValue{-1.};
        
        Int_t fNbins{-1};
        Int_t fBinSize{1};
        
        Int_t fFitMinBin{-1};
        Int_t fFitMaxBin{-1};

        Int_t fNormMinBin{-1};
        
        TString fMode{"Default"};
        
        TString fOutDirName{""};
        ClassDef(Fitter, 2);
        
    };
}

#endif
