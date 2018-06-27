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
#include "TDirectory.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

// Function which does Ratio Method and produces R(t), D(t), S(t):

//TH1F* GetRatioMethodHistFromPseudoExp(TH1F* PseudoExp, double Total_Time, double BinWidth, int BinNum, double binshift TH1F* RatDist_hist, TH1F* Ddist_hist, TH1F* Sdist_hist){

vector<TH1D*> GetRatioMethodHistFromPseudoExp(TH1D* PseudoExpA, TH1D* PseudoExpB, TH1D* PseudoExpC, TH1D* PseudoExpD, double Total_Time, double BinWidth, int BinNum, int choice = 0){
  
  Double_t Toff    = 2100.0;   // Toff = pi/omega -> omega = 2pi/4200 -> Toff would be 2100 ns
  
  int binshift = Toff / BinWidth;

  // Book histograms for R(t), D(t), S(t)
  
  //TH1D* RatDist_hist = new TH1D("" , "" , BinNum - (2 * binshift) , Toff - (BinWidth/2.0) , Total_Time - Toff - (BinWidth/2.0));
  //TH1D* Ddist_hist   = new TH1D("" , "" , BinNum - (2 * binshift) , Toff - (BinWidth/2.0) , Total_Time - Toff - (BinWidth/2.0));
  //TH1D* Sdist_hist   = new TH1D("" , "" , BinNum - (2 * binshift) , Toff - (BinWidth/2.0) , Total_Time - Toff - (BinWidth/2.0));
  TH1D* RatDist_hist = new TH1D("" , "" , BinNum , 0 , Total_Time );
  TH1D* Ddist_hist   = new TH1D("" , "" , BinNum , 0 , Total_Time );
  TH1D* Sdist_hist   = new TH1D("" , "" , BinNum , 0 , Total_Time );

  // Loop over bins of pseudo exp
  
  //for (int ibin(BinWidth/2.0 + binshift); ibin <= (BinNum  - binshift); ibin++ ){
  for (int ibin = binshift+1; ibin < (BinNum  - binshift); ibin++ ){
    
    int ibin_up = ibin + binshift;
    int ibin_do = ibin - binshift;
    
    // Divide bin content by 4, apply shifts and add to correct array
    double d1 = PseudoExpA -> GetBinContent(ibin) ;
    double d2 = PseudoExpB -> GetBinContent(ibin) ;   
    double d3 = PseudoExpC -> GetBinContent(ibin_up) ;
    double d4 = PseudoExpD -> GetBinContent(ibin_do) ;

    double Ddist = d3 + d4 - d1 - d2;
    double Sdist = d1 + d2 + d3 + d4;
    double RatioValue = Ddist/Sdist;

    RatDist_hist -> SetBinContent(ibin, RatioValue);
    Ddist_hist   -> SetBinContent(ibin, Ddist);
    Sdist_hist   -> SetBinContent(ibin, Sdist);

    // Find error on R(t), D(t), S(t)
    double RatioError = sqrt( (1. - pow(RatioValue,2)) / Sdist );
    double DError = sqrt(Sdist);
    
    // Set error on R(t), D(t), S(t)
    
    RatDist_hist -> SetBinError(ibin, RatioError);
    gStyle       -> SetEndErrorSize(2);
    RatDist_hist -> SetMarkerStyle(20);
    RatDist_hist -> SetMarkerSize(0.5);
    RatDist_hist -> SetMarkerColor(kBlack);
    RatDist_hist -> Draw("E1");

    Ddist_hist -> SetBinError(ibin, DError);
    gStyle     -> SetEndErrorSize(2);
    Ddist_hist -> SetMarkerStyle(20);
    Ddist_hist -> SetMarkerSize(0.5);
    Ddist_hist -> SetMarkerColor(kBlack);
    Ddist_hist -> Draw("E1");

    Sdist_hist -> SetBinError(ibin, DError);
    gStyle     -> SetEndErrorSize(2);
    Sdist_hist -> SetMarkerStyle(1);
    Sdist_hist -> SetMarkerSize(0.5);
    Sdist_hist -> SetMarkerColor(kBlack);
    Sdist_hist -> Draw("E1");
  }
  
  vector<TH1D*> tmp = {RatDist_hist, Ddist_hist, Sdist_hist};

  return tmp;

}

int main(){
  
  Double_t w = 1500;
  Double_t h = 1000;
  TCanvas* c2 = new TCanvas("c2", "c2", w, h);
  
  // Reading ROOT file and extracting histograms (choosen by user)
  
  int N_Pseudo_Exp;
  
  cout << "Please enter the number of histograms you want to use for the Ratio Method: " << endl;
  cin >>  N_Pseudo_Exp;
  
  TFile *file = TFile::Open("PseudoExp.root", "READ");
  
  if (file == 0) {
    
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open PseudoExp.root");
    return 1;
    
  }

  TFile* file1 = new TFile("RatioMethodData.root", "RECREATE"); // Create a ROOT file containg Ratio Method data

  for (int i(0); i < (N_Pseudo_Exp); i++){                      // Loop over number of pseudo exp

    // Using stringstream to accommodate the number of pseudo experiements when naming and giving titles to plots
    string Pseudo_Exp_;
    ostringstream conv1;
    conv1 << i;
    string pseudo_num = "Pseudo_Exp_" + conv1.str();
    
    //file.cd(pseudo_num.c_str());
    //file->ls();
    
    //TDirectory* MainDirec = (TDirectory*)file ->Get(pseudo_num.c_str()) ; 
    
    string Pseudo_num;
    ostringstream conv;
    conv << i;
    string name1 = "PseudoA_Slow_Osc_" + conv.str();
    string name2 = "PseudoB_Slow_Osc_" + conv.str();
    string name3 = "PseudoC_Slow_Osc_" + conv.str();
    string name4 = "PseudoD_Slow_Osc_" + conv.str();

    string outname1 = "RatDist_Pseudo_" + conv.str();
    string outname2 = "DDist_Pseudo_" + conv.str();
    string outname3 = "SDist_Pseudo_" + conv.str();
    string outtitle1 = "R(t) for Pseudo Experiement " + conv.str();
    string outtitle2 = "D(t) for Pseudo Experiement " + conv.str();
    string outtitle3 = "S(t) for Pseudo Experiement " + conv.str();
    
    TH1D* h1 = (TH1D*)file -> Get((pseudo_num+"/"+name1).c_str());   // Gets the pseudo wiggle plot
    TH1D* h2 = (TH1D*)file -> Get((pseudo_num+"/"+name2).c_str());   // Gets the pseudo wiggle plot
    TH1D* h3 = (TH1D*)file -> Get((pseudo_num+"/"+name3).c_str());   // Gets the pseudo wiggle plot
    TH1D* h4 = (TH1D*)file -> Get((pseudo_num+"/"+name4).c_str());   // Gets the pseudo wiggle plot
    
    Double_t BinWidth = h1 -> GetBinWidth(1);      // Get bin width of pseudo wiggle plot
    
    Double_t BinNum = h1 -> GetNbinsX();           // Get bin number of pseudo wiggle plot
    
    Double_t Total_Time = BinWidth * BinNum;       // Find total time
    
    vector<TH1D*> tmp = GetRatioMethodHistFromPseudoExp(h1, h2, h3, h4, Total_Time , BinWidth , BinNum , 0);

    tmp.at(0) -> SetNameTitle(outname1.c_str(), outtitle1.c_str());
    tmp.at(0) -> GetXaxis() -> SetTitle("Time (ns)");
    tmp.at(0) -> GetYaxis() -> SetTitle("R(t)");
    
    tmp.at(1) -> SetNameTitle(outname2.c_str(), outtitle2.c_str());
    tmp.at(1) -> GetXaxis() -> SetTitle("Time (ns)");
    tmp.at(1) -> GetYaxis() -> SetTitle("D(t)");
    
    tmp.at(2) -> SetNameTitle(outname3.c_str(), outtitle3.c_str());
    tmp.at(2) -> GetXaxis() -> SetTitle("Time (ns)");
    tmp.at(2) -> GetYaxis() -> SetTitle("S(t)");

    tmp.at(0) -> Write();
    tmp.at(1) -> Write();
    tmp.at(2) -> Write();
    
  }

}
