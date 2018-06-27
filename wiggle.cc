#include <iostream>
#include <fstream>
#include <iomanip>
#include "TH1F.h"
#include "TH1D.h"
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
#include"TDirectory.h"
using namespace std;

// Function for producing wiggle plot:
// All units are in ns

Double_t Total_Time = 30.0 * 4200.0;

Double_t BinWidth = 150;

//int rebinFactor = 1;

Double_t myfunction(Double_t *x, Double_t *par){

  Double_t time     = x[0];
  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.4;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double phase     = TMath::Pi()/2;
  double tau_gamma = tau * gamma;       //Making life easier
  double N = par[0];
  
  //Double_t Npositrons =  N * exp(- time / tau_gamma ) * (1 + A * cos ((omega * time) + phase));
  //          N(t)    =  N(0) exp (-t/(tau * gamma )) * (1 + A * cos(omega * t + phase))
  //return Npositrons;

  //double Ntot      = 2E11;               // Total number of events 
  //double N         = Ntot *  ( (1 - exp( - (1 / tau_gamma) * BinWidth )) / (1 - exp( - (1 / tau_gamma) * Total_Time)) );
  
  Double_t Npositrons =  N * exp(- time / tau_gamma ) * (1 + A * cos ((omega * time) + phase));
  //          N(t)    =  N(0) exp (-t/(tau * gamma )) * (1 + A * cos(omega * t + phase))
  
  //cout << "Npositrons = " << Npositrons << ", time = " << time <<  endl;

  return Npositrons;
  
}

Double_t integral(Double_t minTime, Double_t maxTime){

  double tau       = 2200;              // Lifetime of the muon at rest (ns)
  double gamma     = 29.3;              // Magic gamma
  double A         = 0.4;              // Amplitude
  double omega     = 2 * TMath::Pi() / 4200;   // Time of a single wiggle ~ 4200 ns (omega=2pi/t)
  double phase     = TMath::Pi()/2;
  double tau_gamma = tau * gamma;       //Making life easier

  double exp_min = exp(-minTime/tau_gamma);
  double exp_max = exp(-maxTime/tau_gamma);
  double arg_min = omega * minTime + phase;
  double arg_max = omega * maxTime + phase;
  double firstPart = tau_gamma * (exp_min - exp_max);
  double secondPart_const = A / ( (1/(tau_gamma*tau_gamma)) + omega*omega );
  double secondPart_min = exp_min * (omega * sin(arg_min) - (1/tau_gamma)*cos(arg_min));
  double secondPart_max = exp_max * (omega * sin(arg_max) - (1/tau_gamma)*cos(arg_max));
  double secondPart = secondPart_const * (secondPart_max - secondPart_min);

  
  double integral = firstPart + secondPart;

  return integral;
  
}

// Function for producing a slow cosine oscillation

Double_t oscillate2(Double_t x){

  Float_t time = x;
  double A     = 0.05;               // Amplitude
  double omega = 2 * M_PI / (30*4200);  // long wavelength

  Double_t oscillate2 = 0.8 * ( 1 + A * cos(omega * time));

  return oscillate2;

}

int main(){

  TFile* file = new TFile("PseudoExp.root", "RECREATE"); // Create a ROOT file containing all pseudo experiments
  
  Double_t w = 1500;
  Double_t h = 1000;

  double Ntot = 2E11;
  
  TF1* f1 = new TF1("f1", myfunction, 0, Total_Time, 1); // Book histogram for ideal wiggle plot
  f1->SetParameter(0,1);
  double functionIntegral = f1->Integral(0, Total_Time);
  double analyticIntegral = integral(0, Total_Time);
  double N0 = Ntot/functionIntegral;
  f1->SetParameter(0,N0);
  
  //TF1 *f1 = new TF1("f1", myfunction, 0, Total_Time, 1); // Book histogram for ideal wiggle plot
  //f1->SetParameter(0, N0);
  
  // Loop over full time of experiment
  
  Double_t par[] = {1};
  Double_t Nfunc[int(Total_Time / BinWidth)];
  double integralSum = 0;
  double analyticIntegralSum = 0;

  Double_t Oscfunc2[int(Total_Time / BinWidth)];

  
  for (int ibin = 0; ibin < (Total_Time / BinWidth); ibin++){

    double time = BinWidth * ibin + BinWidth/2;
    //Nfunc[ibin] = N0 * BinWidth * myfunction(&time, par); // Time distribution for decay positrons

    //cout << "time = " << time << endl;
    
    double binIntegral = f1->Integral(time-(BinWidth/2), time+(BinWidth/2));
    double analyticBinIntegral = integral(time-(BinWidth/2), time+(BinWidth/2));
    integralSum += binIntegral;
    analyticIntegralSum += analyticBinIntegral;
    Nfunc[ibin] = analyticBinIntegral * N0;

    double centreValue = myfunction(&time, par);

    //cout << "bin integral = " << binIntegral << ", centre value = " << centreValue * BinWidth << ", difference = " << binIntegral - centreValue*BinWidth << endl;
    
    //cout << " bin integral = " << binIntegral << ", analytic bin integral = " << analyticBinIntegral << endl;

    Oscfunc2[ibin] = oscillate2(time);        // Fast Cosine oscillation (short wavelength)
      
  }

  cout << " *************** Function integral = " << functionIntegral << ", sum = " << integralSum << endl;
  cout << " *************** Analytic function integral = " << analyticIntegral << ", sum = " << analyticIntegralSum << endl;
    
  TRandom3 *r = new TRandom3(54321); // Set seed to get same random numbers (0 means no seed)

  // Loop over number of psuedo experiments (choosen by the user)

  int N_Pseudo_Exp;
  
  cout << "Please enter the number of pseudo experiements: " << endl;
  cin >>  N_Pseudo_Exp;   

  for (int i=0; i <= (N_Pseudo_Exp); i++){

    string Pseudo_Exp_;
    ostringstream conv1;
    conv1 << i;
    string pseudo_num = "Pseudo_Exp_" + conv1.str();

    //Create a main directory
    TDirectory *MainDirec = file->mkdir(pseudo_num.c_str());
    MainDirec->cd();



    
    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    
    string Pseudo_num;
    ostringstream conv;
    conv << i;
    string t1name  = "PseudoA_" + conv.str();
    string t1title = "PseudoA Wiggle Plot " + conv.str();
    string t2name  = "PseudoB_" + conv.str();
    string t2title = "PseudoB Wiggle Plot " + conv.str();
    string t3name  = "PseudoC_" + conv.str();
    string t3title = "PseudoC Wiggle Plot " + conv.str();
    string t4name  = "PseudoD_" + conv.str();
    string t4title = "PseudoD Wiggle Plot " + conv.str();
    string tTotname  = "PseudoTot_" + conv.str();
    string tTottitle = "PseudoTot Wiggle Plot " + conv.str();


    string tSlowOscname  = "osc_long_" + conv.str();
    string tSlowOsctitle = "Slow Cosine Oscillation Plot " + conv.str();

    string tWiggleASlowOscname  = "PseudoA_Slow_Osc_" + conv.str();
    string tWiggleASlowOsctitle = "PseudoA Wiggle with Slow Oscillation " + conv.str();

    string tWiggleBSlowOscname  = "PseudoB_Slow_Osc_" + conv.str();
    string tWiggleBSlowOsctitle = "PseudoB Wiggle with Slow Oscillation " + conv.str();

    string tWiggleCSlowOscname  = "PseudoC_Slow_Osc_" + conv.str();
    string tWiggleCSlowOsctitle = "PseudoC Wiggle with Slow Oscillation " + conv.str();

    string tWiggleDSlowOscname  = "PseudoD_Slow_Osc_" + conv.str();
    string tWiggleDSlowOsctitle = "PseudoD Wiggle with Slow Oscillation " + conv.str();

    string tWiggleTotSlowOscname  = "PseudoTot_Slow_Osc_" + conv.str();
    string tWiggleTotSlowOsctitle = "PseudoTot Wiggle with Slow Oscillation " + conv.str();
    

    
    string tDiffname  = "Pseudo_diff_" + conv.str();
    string tDifftitle = "Wiggle Difference Plot " + conv.str();

    // Book a histogram for pseudo wiggle plot

    //TH1D* t1 = new TH1D("t1", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    //TH1D* t2 = new TH1D("t2", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    //TH1D* t3 = new TH1D("t3", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    //TH1D* t4 = new TH1D("t4", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    //TH1D* tTot = new TH1D("tT0", "", int(Total_Time / BinWidth), -(BinWidth/2), Total_Time - (BinWidth/2));
    
    TH1D* t1 = new TH1D("t1", "", int(Total_Time / BinWidth), 0, Total_Time);
    TH1D* t2 = new TH1D("t2", "", int(Total_Time / BinWidth), 0, Total_Time);
    TH1D* t3 = new TH1D("t3", "", int(Total_Time / BinWidth), 0, Total_Time);
    TH1D* t4 = new TH1D("t4", "", int(Total_Time / BinWidth), 0, Total_Time);
    TH1D* tTot = new TH1D("tT0", "", int(Total_Time / BinWidth), 0, Total_Time);


    
        
    t1 -> SetNameTitle(t1name.c_str(), t1title.c_str());
    t1 -> GetXaxis() -> SetTitle("Time (ns)");
    t1 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t2 -> SetNameTitle(t2name.c_str(), t2title.c_str());
    t2 -> GetXaxis() -> SetTitle("Time (ns)");
    t2 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t3 -> SetNameTitle(t3name.c_str(), t3title.c_str());
    t3 -> GetXaxis() -> SetTitle("Time (ns)");
    t3 -> GetYaxis() -> SetTitle("Number of Positrons N");

    t4 -> SetNameTitle(t4name.c_str(), t4title.c_str());
    t4 -> GetXaxis() -> SetTitle("Time (ns)");
    t4 -> GetYaxis() -> SetTitle("Number of Positrons N");

    tTot -> SetNameTitle(tTotname.c_str(), tTottitle.c_str());
    tTot -> GetXaxis() -> SetTitle("Time (ns)");
    tTot -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Books a histogram for slow cosine oscillation

    TH1F* tSlowOsc = new TH1F("", "", int(Total_Time / BinWidth), 0, Total_Time);

    tSlowOsc -> SetNameTitle(tSlowOscname.c_str(), tSlowOsctitle.c_str());
    tSlowOsc -> GetXaxis() -> SetTitle("Time (ns)");
    tSlowOsc -> GetYaxis() -> SetTitle("Acceptance");

    // Books a histogram for pseudo wiggleA with slow cosine oscillation

    TH1F* tWiggleASlowOsc = new TH1F("", "", int(Total_Time / BinWidth), 0, Total_Time);

    tWiggleASlowOsc -> SetNameTitle(tWiggleASlowOscname.c_str(), tWiggleASlowOsctitle.c_str());
    tWiggleASlowOsc -> GetXaxis() -> SetTitle("Time (ns)");
    tWiggleASlowOsc -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Books a histogram for pseudo wiggleB with slow cosine oscillation

    TH1F* tWiggleBSlowOsc = new TH1F("", "", int(Total_Time / BinWidth), 0, Total_Time);

    tWiggleBSlowOsc -> SetNameTitle(tWiggleBSlowOscname.c_str(), tWiggleBSlowOsctitle.c_str());
    tWiggleBSlowOsc -> GetXaxis() -> SetTitle("Time (ns)");
    tWiggleBSlowOsc -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Books a histogram for pseudo wiggleC with slow cosine oscillation

    TH1F* tWiggleCSlowOsc = new TH1F("", "", int(Total_Time / BinWidth), 0, Total_Time);

    tWiggleCSlowOsc -> SetNameTitle(tWiggleCSlowOscname.c_str(), tWiggleCSlowOsctitle.c_str());
    tWiggleCSlowOsc -> GetXaxis() -> SetTitle("Time (ns)");
    tWiggleCSlowOsc -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Books a histogram for pseudo wiggleD with slow cosine oscillation

    TH1F* tWiggleDSlowOsc = new TH1F("", "", int(Total_Time / BinWidth), 0, Total_Time);

    tWiggleDSlowOsc -> SetNameTitle(tWiggleDSlowOscname.c_str(), tWiggleDSlowOsctitle.c_str());
    tWiggleDSlowOsc -> GetXaxis() -> SetTitle("Time (ns)");
    tWiggleDSlowOsc -> GetYaxis() -> SetTitle("Number of Positrons N");

    // Books a histogram for pseudo wiggleTot with slow cosine oscillation

    TH1F* tWiggleTotSlowOsc = new TH1F("", "", int(Total_Time / BinWidth), 0, Total_Time);

    tWiggleTotSlowOsc -> SetNameTitle(tWiggleTotSlowOscname.c_str(), tWiggleTotSlowOsctitle.c_str());
    tWiggleTotSlowOsc -> GetXaxis() -> SetTitle("Time (ns)");
    tWiggleTotSlowOsc -> GetYaxis() -> SetTitle("Number of Positrons N");
    
    
    // Book a histogram for difference between ideal wiggle and pseudo wiggle plots
    
    TH1D* tDiff = new TH1D("", "", 100, -10, 10);
    
    tDiff -> SetNameTitle(tDiffname.c_str(), tDifftitle.c_str());
    tDiff -> GetXaxis() -> SetTitle("N/#sqrt{N}");
    tDiff -> GetYaxis() -> SetTitle("Number of Entries");    
    
    for (int ibin = 0; ibin < (Total_Time / BinWidth); ibin++){
      
      double nHits = Nfunc[ibin]/4.0;
      double error = sqrt( Nfunc[ibin]/4.0 );
      
      Double_t osc2 = Oscfunc2[ibin];
      
      Double_t PseudoWiggle1 = r->Gaus(nHits, error);  // Gaussian smearing
      Double_t PseudoWiggle2 = r->Gaus(nHits, error);  // Gaussian smearing
      Double_t PseudoWiggle3 = r->Gaus(nHits, error);  // Gaussian smearing
      Double_t PseudoWiggle4 = r->Gaus(nHits, error);  // Gaussian smearing
      
      Double_t PseudoWiggle = PseudoWiggle1 + PseudoWiggle2 + PseudoWiggle3 + PseudoWiggle4;
      t1 -> SetBinContent(ibin+1, PseudoWiggle1);
      t2 -> SetBinContent(ibin+1, PseudoWiggle2);
      t3 -> SetBinContent(ibin+1, PseudoWiggle3);
      t4 -> SetBinContent(ibin+1, PseudoWiggle4);
      
      tTot -> SetBinContent(ibin+1, PseudoWiggle);

      tSlowOsc -> SetBinContent(ibin+1, osc2);

      tWiggleASlowOsc -> SetBinContent(ibin+1, osc2 * PseudoWiggle1);
      tWiggleBSlowOsc -> SetBinContent(ibin+1, osc2 * PseudoWiggle2);
      tWiggleCSlowOsc -> SetBinContent(ibin+1, osc2 * PseudoWiggle3);
      tWiggleDSlowOsc -> SetBinContent(ibin+1, osc2 * PseudoWiggle4);
      tWiggleTotSlowOsc -> SetBinContent(ibin+1, osc2 * PseudoWiggle);
      //t1 -> SetBinContent(ibin+1, nHits);
      //t2 -> SetBinContent(ibin+1, nHits);
      //t3 -> SetBinContent(ibin+1, nHits);
      //t4 -> SetBinContent(ibin+1, nHits);
      
      //tTot -> SetBinContent(ibin+1, Nfunc[ibin]);

      //if (ibin == 100 ) {
      //	cout << "p1: " << PseudoWiggle1 << " p2: " << PseudoWiggle2 << "p3: " << PseudoWiggle3 << " p4: " << PseudoWiggle4 << "\n";
      //	cout << "and p1: " << t1->GetBinContent(ibin+1) << " p2: " << t2->GetBinContent(ibin+1) << "p3: " << t3->GetBinContent(ibin+1) << " p4: " << t4->GetBinContent(ibin+1) << "\n";	
      //}
      
      Double_t diff =  (Nfunc[ibin] - PseudoWiggle) / (sqrt(Nfunc[ibin])); // Difference between ideal and pseudo wiggle
      tDiff -> Fill(diff);


      
    }

    
    // Write TH1F to ROOT file
    t1 -> Write();
    t2 -> Write();
    t3 -> Write();
    t4 -> Write();
    tTot -> Write();
    tDiff -> Write();

    tSlowOsc -> Write();
    tWiggleASlowOsc -> Write();
    tWiggleBSlowOsc -> Write();
    tWiggleCSlowOsc -> Write();
    tWiggleDSlowOsc -> Write();
    tWiggleTotSlowOsc -> Write();


      
    
    delete t1;
    delete t2;
    delete t3;
    delete t4;
    delete tTot;
    delete tDiff;

    delete tSlowOsc;
    delete tWiggleASlowOsc;
    delete tWiggleBSlowOsc;
    delete tWiggleCSlowOsc;
    delete tWiggleDSlowOsc;
    delete tWiggleTotSlowOsc;


    MainDirec->Close();
  }



    

  
  // Draw ideal wiggle plot
  
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);
  f1 -> SetTitle("Ideal Wiggle Plot");
  f1 -> GetXaxis() -> SetTitle("Time (ns)");
  f1 -> GetYaxis() -> SetTitle("Number of Positrons N");
  //f1 -> GetYaxis() -> SetRangeUser(0, 2E11* 1.1);
  f1 -> Draw();
  c2 -> SaveAs("IdealFit.eps");

  delete c2;
  delete f1;
  file ->Close();
  
}
