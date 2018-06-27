#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom3.h"
#include "TVirtualFFT.h"
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TFitResult.h"
#include <TMatrixDSym.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include<TDirectory.h>
using namespace std;

Double_t myfunction(Double_t *x, Double_t *par){

  Double_t time   = x[0];
  Double_t A     = par[0];
  Double_t omega = par[1];
  Double_t phase = par[2];
  Double_t delta = par[3];

  Double_t Ratio =  par[0] * cos((par[1] * x[0]) + par[2]) + par[3];

  return Ratio;
  
}

int main(){

  double w = 1500;
  double h = 1000;
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);
  
  // Reading Root file and extracting histograms (choosen by user)
  
  int N_Pseudo_Exp;
  
  cout << "Please enter the number of histograms you want to use for the three parameter fit: " << endl;
  cin >>  N_Pseudo_Exp;
  
  TFile *file = TFile::Open("RatioMethodData.root", "READ");
  
  if (file == 0) {
    
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open PseudoExp.root");
    return 1;
    
  }
  
  TFile* file3 = new TFile("ThreeParameterFit.root", "RECREATE");  // Create a ROOT file containg three parameter fits of the Ratio Method plots

  /////// Histos outside the pseudo directories
  
  TH1D* h_amp = new TH1D("amp","", 1000, 0.049, 0.051);
  TH1D* h_omega = new TH1D("omega","", 1000, 0.0014955, 0.0014965);
  TH1D* h_phase = new TH1D("phase","", 1000, 0, 6.2);
  
  TH1D* h_prob = new TH1D("p-value","", 110, 0, 1.1);
  TH1D* h_chiSqFit = new TH1D("chiSqFit","", 1000, 0, 1500.);
  TH1D* h_chiSqPerNDFFit = new TH1D("chiSqPerNDFFit","", 100, 0, 5.);
  TH1D* h_NDFFit = new TH1D("ndfFit","", 1000, 0 -0.5, 1000. - 0.5);
  
  TH1D* h_chiSq = new TH1D("chiSq","", 1000, 0, 1500.);
  TH1D* h_chiSqPerNDF = new TH1D("chiSqPerNDF","", 100, 0, 5.);
  TH1D* h_NDF = new TH1D("ndf","", 1000, 0 - 0.5, 1000. - 0.5);

     
  string PullFittedADistName = "Pull_Distribution_Fitted_A_Pseudo";    
  TH1D* h5 = new TH1D(PullFittedADistName.c_str(),"", 100, -10, 10);
  
  
  string PullFittedOmegaDistName = "Pull_Distribution_Fitted_Omega_Pseudo";    
  TH1D* h6 = new TH1D(PullFittedOmegaDistName.c_str(),"", 100, -10, 10);
  
  
  string PullFittedPhiDistName = "Pull_Distribution_Fitted_Phi_Pseudo";    
  TH1D* h7 = new TH1D(PullFittedPhiDistName.c_str(),"", 100, -10, 10);
  
  string DiffA = "Diff_A";
  TH1D* h8 = new TH1D(DiffA.c_str(),"", 100, -0.03E-3, 0.03E-3);

  string DiffOmegaA = "Diff_Omega_A";
  TH1D* h9 = new TH1D(DiffOmegaA.c_str(),"", 100, -40E-9, 40E-9);

  string DiffPhi = "Diff_Phi";
  TH1D* h10 = new TH1D(DiffPhi.c_str(),"", 100, -0.001, 0.001);

  
  ///////////////////////////////////////////////////////////////////////////////////
  
  TF1 *fit = new TF1 ("fit", myfunction, 2100, 124000, 4);   // Three parameter fit

  fit -> SetParameter(0, 1.01*0.4);
  fit -> SetParameter(1, 0.99*(2*M_PI/4200.0));
  fit -> SetParameter(2, 1.02*3*M_PI/2.0);
  //  fit -> SetParameter(3, 2.87E-4);
    
  fit -> SetParLimits(0, 0.1, 1.0);
  fit -> SetParLimits(1, 0.0005, 0.0030);
  fit -> SetParLimits(2, 0.0, 2*M_PI);
  //  fit -> SetParLimits(3, 1E-4, 3E-4);

  fit-> FixParameter(3, 2.87E-4);
  
  fit -> SetNpx(10000);
  fit -> SetNDF(4);
    
  fit -> SetParName(0, "A");
  fit -> SetParName(1, "Omega");
  fit -> SetParName(2, "Phi");
  fit -> SetParName(3, "Delta");
  
  fit -> SetLineColor(kRed);
  fit -> SetLineStyle(1);
  fit -> SetLineWidth(1);

  for (int i(0); i < (N_Pseudo_Exp); i++){                         // Loop over number of pseudo exp

    string Pseudo_Exp_;
    ostringstream conv1;
    conv1 << i;
    string pseudo_num = "Pseudo_Exp_" + conv1.str();
          
    //Create a main directory
    TDirectory *MainDirec = file3->mkdir(pseudo_num.c_str());
    MainDirec->cd();

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string RatDist_Pseudo_num;
    ostringstream conv2;
    conv2 << i;
    string name = "RatDist_Pseudo_" + conv2.str();    
    
    TH1D* h1 = (TH1D*)file -> Get(name.c_str());                // Gets the ratio method plot
      
    double BinWidth = h1 -> GetBinWidth(1);                     // Gets bin width of ratio method plot

    double BinNum = h1 -> GetNbinsX();                          // Gets bin number of ratio method plot

    double Total_Time = h1->GetBinLowEdge(BinNum) + BinWidth;  //2100 + BinWidth * BinNum;                      // Find total time
    
    //h1-> GetXaxis() -> SetRangeUser(0,4400);
    
    h1 -> Fit(fit, "EMRI+");
    
    gStyle->SetOptFit(1111);

    h1 -> Write();

    //ResDirec ->cd();

    /////// Histos inside pseudo directories
    ostringstream count;
    count << i;
    string ResName = "Residual_Pseudo_" + count.str();    
    TH1D* h2 = new TH1D(ResName.c_str(),"", BinNum , h1->GetBinLowEdge(1),  h1->GetBinLowEdge(1) + Total_Time);

    ostringstream count2;
    count2 << i;
    string PullName = "Pull_Pseudo_" + count2.str();    
    TH1D* h3 = new TH1D(PullName.c_str(),"", BinNum , h1->GetBinLowEdge(1),  h1->GetBinLowEdge(1) + Total_Time);

    ostringstream count3;
    count3 << i;
    string PullDistName = "Pull_Distribution_Pseudo_" + count3.str();    
    TH1D* h4 = new TH1D(PullDistName.c_str(),"", 64, -15, 15);

    double chiSq = 0.0;
    int binsInChiSq = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //get our histogram and loop over bins
    for (int ibin(1); ibin < BinNum; ibin++){
      
      double time = h1->GetBinCenter(ibin);
      double measured = h1->GetBinContent(ibin);
      double error = h1->GetBinError(ibin);
      double fitted = fit->Eval(time);
      
      //cout << "time: " << time << ", meas: " << measured << ", err: " << error << ", fitted: " << fitted << "\n";

    
	 
      if (error > 0){
	
	h2 -> SetBinContent(ibin,(measured-fitted));        //plot residuals

	h3 -> SetBinContent(ibin,(measured-fitted)/error);  //Plot pull

	h4 -> Fill((measured-fitted)/error);                //Plot pull distribution

	//
	

		
	//
	chiSq += pow((measured-fitted),2) / pow(error,2);

	binsInChiSq++;
      }
      
    }


    ////FFT////
    TH1 *hm =0;
    TVirtualFFT::SetTransform(0);
    hm = h2->FFT(hm, "MAG");
    //hm->SetTitle("Fourier Transform of Residuals");
    //hm->Draw();
    //hm->Write();

    int n = hm -> GetNbinsX();
    double binwidth = h2 -> GetXaxis() -> GetBinWidth(0);
    double samFreq = 1/binwidth;
    double xmax = hm -> GetXaxis() -> GetXmax();    //binwidth;
    double xmin = hm -> GetXaxis() -> GetXmin();    //binwidth;

    TH1F *hm2 = new TH1F("hm", "title", (n/2) + 1, xmin, xmax);

    for (int i = 1 ; i <= ((n/2) + 1) ; i++) {
      double y0 = (hm -> GetBinContent(i) - hm -> GetBinContent((n+1) - i));
      double y = sqrt(y0*y0);
      double ynew = y/(sqrt(n));
      double x = hm -> GetXaxis() -> GetBinCenter(i);
      hm2 -> Fill(x,ynew);

    }


    hm2 -> GetXaxis() -> SetLimits(0, samFreq);
    hm2 -> SetNameTitle("fft_residuals", "Fourier Transform of Residuals");
    hm2 -> GetXaxis() -> SetTitle("Frequency (GHz)");
    hm2 -> GetYaxis() -> SetTitle("Magnitude");
    hm2 -> SetFillColor(kBlue);

    //hm2 -> GetXaxis() -> SetRangeUser(0.0001, 0.001);
    hm2 -> GetXaxis() -> SetRangeUser(0, 0.001);
    hm2 -> Draw();

    TF1 *FFTGaussFit1 = new TF1("FFTGaussFit1", "gaus", 0.00021, 0.00026);
    int FFTbinmax1 = hm2 -> GetMaximumBin();
    double FFTxbinmax1 = hm2 -> GetXaxis() -> GetBinCenter(FFTbinmax1);
    Double_t FFTsigmabinmax1 = hm2 -> GetRMS(1);
    FFTGaussFit1 -> SetParameter(0, FFTbinmax1);
    FFTGaussFit1 -> SetParameter(1, FFTxbinmax1);
    FFTGaussFit1 -> SetParameter(2, FFTsigmabinmax1);
    FFTGaussFit1 -> SetLineColor(2);
    FFTGaussFit1 -> SetLineWidth(3);
    FFTGaussFit1 -> SetLineStyle(1);
    FFTGaussFit1 -> SetNpx(10000);
    FFTGaussFit1 -> SetNDF(3);

    hm2->Fit("FFTGaussFit1", "MR");
    gStyle->SetOptFit(1);
    
    
    double fitted_A      = fit -> GetParameter(0);
    double fitted_omega  = fit -> GetParameter(1);
    double fitted_phi    = fit -> GetParameter(2);
    
    double fitted_AError      = fit -> GetParError(0);
    double fitted_omegaError  = fit -> GetParError(1);
    double fitted_phiError    = fit -> GetParError(2);
    
    //Inputted (ideal) values
    
    double ideal_A = 0.4;
    double ideal_omega = 2 * M_PI / 4200;
    double ideal_phi = 3*M_PI / 2;
    
    cout << "fitted A " << fitted_A << " ideal A " << ideal_A << endl;
    cout << "fitted omega " << fitted_omega << " ideal omega " << ideal_omega << endl;
    cout << "fitted phi " << fitted_phi << " ideal phi " << ideal_phi << endl; 
    
    h5->Fill((fitted_A-ideal_A)/fitted_AError);
    
    h6->Fill((fitted_omega-ideal_omega)/fitted_omegaError);
    
    h7->Fill((fitted_phi-ideal_phi)/fitted_phiError);
    
    h8->Fill((ideal_A - fitted_A));

    h9->Fill((ideal_omega - fitted_omega));

    h10->Fill((ideal_phi - fitted_phi));
    
    
  
    
      
    h2 -> Write();
    hm2->Write();
    h3 -> Write();
    h4 -> Write();

   
    h2->Delete();
    hm2->Delete();
    h3->Delete();
    h4->Delete();

    //h5->Delete();
    //h6->Delete();
    //h7->Delete();
    
    //measured
    h_chiSq->Fill(chiSq);
    h_chiSqPerNDF->Fill(chiSq / double(binsInChiSq));
    h_NDF -> Fill(binsInChiSq);

    //from fit
    h_prob->Fill(fit->GetProb());
    h_chiSqFit->Fill(fit->GetChisquare());
    h_chiSqPerNDFFit->Fill( fit->GetChisquare() / double(fit->GetNDF()));
    h_NDFFit->Fill(fit->GetNDF());

    h_amp->Fill(fit->GetParameter(0));
    h_omega->Fill(fit->GetParameter(1));
    h_phase->Fill(fit->GetParameter(2));

    
    
    //cout << "chi sq test, us: " << chiSq << " perndf: " << chiSq/ double(binsInChiSq) << " from fit: " << fit->GetChisquare() << ", per ndf " << fit->GetChisquare() / fit->GetNDF() << " with prob " << fit->GetProb() << ", total time: " << Total_Time << "\n"; 
    
    MainDirec->Close();
    

  }

  file3->cd();

  string fitOpt = "MR";
  //string fitOpt = "MRL";


  TF1 *GaussFit1 = new TF1("GaussFit1", "gaus", -4, 4);
  int binmax1 = h5 -> GetMaximumBin();
  double xbinmax1 = h5 -> GetXaxis() -> GetBinCenter(binmax1);
  Double_t sigmabinmax1 = h5 -> GetRMS(1);
  GaussFit1 -> SetParameter(0, binmax1);
  GaussFit1 -> SetParameter(1, xbinmax1);
  GaussFit1 -> SetParameter(2, sigmabinmax1);
  GaussFit1 -> SetLineColor(2);
  GaussFit1 -> SetLineWidth(3);
  GaussFit1 -> SetLineStyle(1);

    
  h5->Fit("GaussFit1", fitOpt.c_str());
  gStyle->SetOptFit(111);

  TF1 *GaussFit2 = new TF1("GaussFit2", "gaus", -4, 4);
  int binmax2 = h6 -> GetMaximumBin();
  double xbinmax2 = h6 -> GetXaxis() -> GetBinCenter(binmax2);
  Double_t sigmabinmax2 = h6 -> GetRMS(1);
  GaussFit2 -> SetParameter(0, binmax2);
  GaussFit2 -> SetParameter(1, xbinmax2);
  GaussFit2 -> SetParameter(2, sigmabinmax2);
  GaussFit2 -> SetLineColor(2);
  GaussFit2 -> SetLineWidth(3);
  GaussFit2 -> SetLineStyle(1);
    
  h6->Fit("GaussFit2", fitOpt.c_str());
  gStyle->SetOptFit(111);
  
  TF1 *GaussFit3 = new TF1("GaussFit3", "gaus", -4, 4);
  int binmax3 = h7 -> GetMaximumBin();
  double xbinmax3 = h7 -> GetXaxis() -> GetBinCenter(binmax3);
  Double_t sigmabinmax3 = h7 -> GetRMS(1);
  GaussFit3 -> SetParameter(0, binmax3);
  GaussFit3 -> SetParameter(1, xbinmax3);
  GaussFit3 -> SetParameter(2, sigmabinmax3);
  GaussFit3 -> SetLineColor(2);
  GaussFit3 -> SetLineWidth(3);
  GaussFit3 -> SetLineStyle(1);
  
  h7->Fit("GaussFit3", fitOpt.c_str());
  gStyle->SetOptFit(1);

 
  h_chiSq->Write();
  h_chiSqPerNDF->Write();
  h_NDF -> Write();
  h_prob->Write();
  h_chiSqFit->Write();
  h_chiSqPerNDFFit->Write();
  h_NDFFit->Write();
  h_amp->Write();
  h_omega->Write();
  h_phase->Write();
  h5 -> Write();
  h6 -> Write();
  h7 -> Write();
  h8 -> Write();
  h9 -> Write();
  h10 -> Write();
  
  h_chiSq->Delete();
  h_chiSqPerNDF->Delete();
  h_NDF->Delete();
  h_prob->Delete();
  h_chiSqFit->Delete();
  h_chiSqPerNDFFit->Delete();
  h_NDFFit->Delete();
  h_amp->Delete();
  h_omega->Delete();
  h_phase->Delete();
  h5->Delete();
  h6->Delete();
  h7->Delete();
  h8->Delete();
  h9->Delete();
  h10->Delete();
  
  
}


