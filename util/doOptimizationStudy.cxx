#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2.h"
#include "TLine.h"
#include "TGraph2D.h"
#include "TPad.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TPolyMarker.h"
#include "TGraph.h"
#include "TMath.h"

#include "SusyAnalysis/AtlasStyle.h"
#include "SusyAnalysis/AtlasLabels.h"
#include "SusyAnalysis/AtlasUtils.h"
#include "RooStats/NumberCountingUtils.h"

#define DEBUG false
#define NMODECPOINTS 19
#define NDLISLEPPOINTS 1

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Forward declerations
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void doSingleRegionOptimizationStudy ( 
                                       TString region, 
                                       TString chan, 
                                       std::map<TString,double>  &significanceMap, 
                                       std::map<TString,double>  &yieldMap, 
                                       bool doPlot, 
                                       std::string model
                                     );
void makePlot ( 
                std::map<TString,double> significanceMap, 
                std::map<TString,double> yieldMap, 
                TString region, 
                TString channel, 
                int plotMode,
                std::string model
              );

double getPrediction( 
                      TString source, 
                      TString region, 
                      TString channel, 
                      int sysMode
                    );

void getModeCDSIDs( 
                    TString dsids[] 
                  );

void getDLiSlepDSIDs (
                       TString dsids[] 
                     );

void   set_plot_style();

TGraph* ContourGraph(
                      TH2D* hist
                    );

TFile* m_inputFile = nullptr;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Function
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void doOptimizationStudy(TString region = "SRmT2b", TString chan = "EE", std::string model = "ModeC", bool multipleSR = true) 
int main(int argc, char** argv)
{
  // Inputs
  std::string regS  = "SRmT2,180";
  std::string chanS = "All";
  std::string model = "ModeC";
  std::string input = "/data/uclhc/uci/user/amete/analysis_n0222_run/hfts/mc15_13TeV_powheg.root";
  bool multipleSR   = true;

  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--region") == 0)
      regS       = argv[++i];
    else if (strcmp(argv[i], "--channel") == 0)
      chanS      = argv[++i];
    else if (strcmp(argv[i], "--model") == 0)
      model      = argv[++i];
    else if (strcmp(argv[i], "--multiSR") == 0)
      multipleSR = true; 
    else if (strcmp(argv[i], "--input") == 0)
      input      = argv[++i];
    else {
      return 0;
    }
  }

  // Print
  cout << "=======================================================" << endl;
  cout << "Running doOptimizationStudy w/ arguments:"               << endl;
  cout << "=======================================================" << endl;
  cout << "Input   = " << input      << endl;
  cout << "Region  = " << regS       << endl;
  cout << "Channel = " << chanS      << endl;
  cout << "Model   = " << model      << endl;
  cout << "multiSR = " << multipleSR << endl;
  cout << "=======================================================" << endl;

  // Convert 
  TString region; region.Form("%s",regS.c_str());
  TString chan;   chan.Form("%s",chanS.c_str());
      
  // Open file
  if(m_inputFile == nullptr) {
    m_inputFile = new TFile(input.c_str(),"READ");
  }

  if(!multipleSR) {
    std::map<TString,double> significanceMap, yieldMap;
    doSingleRegionOptimizationStudy(region, chan, significanceMap, yieldMap, true, model);
  } else {
    // Signal Regions
    bool orthogonal = false; // if false pwc otherwise stat. combination
    unsigned nSRs = 3;
    //TString regions[2] = { "SRmT2,90" ,"SRmT2,120" };
    //TString regions[2] = { "SRmT2,90,120" ,"SRmT2,120" };
    TString regions[3] = { "SRmT2,90" ,"SRmT2,120" ,"SRmT2,150" };
    //TString regions[3] = { "SRmT2,90,120" ,"SRmT2,120,150" ,"SRmT2,150" };

    // >> Don't touch anything below!!!!

    unsigned nSigPoints = 0.;
    if(model == "ModeC") nSigPoints = NMODECPOINTS;
    else if (model == "DLiSlep") nSigPoints = NDLISLEPPOINTS;
    else {
      std::cout << "Unknown mode " << model << ", quitting..." << std::endl;
      return 0;
    }

    // Calculate Significances for each SR
    std::vector< std::map<TString,double> > significances, yields; 
    for(unsigned int i=0; i<nSRs; i++) { 
      significances.push_back(std::map<TString,double>());
      yields.push_back(std::map<TString,double>());
      doSingleRegionOptimizationStudy(regions[i], "All", significances.at(i), yields.at(i),false,model); 
    }
    // Find the best SR per Signal Point
    std::map<TString,double> finalSignificances, finalYields;

    TString* SignalDSIDs = new TString[nSigPoints];
    if(model == "ModeC") getModeCDSIDs(SignalDSIDs);
    else if (model == "DLiSlep") { std::cout << "Model not implemented yet, quitting ..." << std::endl; return 0; getDLiSlepDSIDs(SignalDSIDs); }

    for(unsigned int i=0; i<nSigPoints; ++i) {
      if(DEBUG) std::cout << "Finding best SR for " << SignalDSIDs[i] << std::endl;
      double sign = 0., yield = 0.; int bestSR = -1;
      for(unsigned int j=0; j<nSRs; ++j) {
        if(!orthogonal) {
          if(significances.at(j)[SignalDSIDs[i]]>sign) {
            bestSR = j;
            sign   = significances.at(j)[SignalDSIDs[i]];
            yield  = yields.at(j)[SignalDSIDs[i]];
          }
        } else {
          std::cout << "Current SR " << j << " has significance " << significances.at(j)[SignalDSIDs[i]] << std::endl;
          sign   += pow(significances.at(j)[SignalDSIDs[i]],2);
        }
      }
      finalSignificances[SignalDSIDs[i]] = !orthogonal ? sign : sqrt(sign);
      finalYields[SignalDSIDs[i]]        = yield;
      if(DEBUG && !orthogonal) std::cout << "Best SR is " << bestSR << " with Significance " << finalSignificances[SignalDSIDs[i]] << std::endl;
      if(/*DEBUG &&*/ orthogonal) std::cout << "Final Significance is " << finalSignificances[SignalDSIDs[i]] << std::endl;
    }
    // Make Plot
    if(!orthogonal) makePlot(finalSignificances,finalYields,"PWC","All",0,model);
    else            makePlot(finalSignificances,finalYields,"Combined","All",0,model);
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Perform significance calculation for a single bin
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void doSingleRegionOptimizationStudy( 
                                      TString region, 
                                      TString chan, 
                                      std::map<TString,double>  &significanceMap, 
                                      std::map<TString,double>  &yieldMap, 
                                      bool doPlot, 
                                      std::string model
                                    )
{
  // Determine the number of signal points
  unsigned int nSigPoints = 0.;

  if(model == "ModeC") nSigPoints = NMODECPOINTS;
  else if (model == "DLiSlep") nSigPoints = NDLISLEPPOINTS;
  else {
    std::cout << "Unknown mode " << model << ", quitting..." << std::endl;
    return;
  }

  TString* SignalDSIDs = new TString[nSigPoints];
  if(model == "ModeC") getModeCDSIDs(SignalDSIDs);
  else if (model == "DLiSlep") getDLiSlepDSIDs(SignalDSIDs);

  double scaleFactor     = 10000.;
  double realtiveBGError = 0.30;
  if(region.Contains("SRmT2,90"))       realtiveBGError = 0.15;
  else if(region.Contains("SRmT2,120")) realtiveBGError = 0.25;
  else if(region.Contains("SRmT2,150")) realtiveBGError = 0.50;
  else if(region.Contains("SRmT2,180")) realtiveBGError = 0.75;
  int plotMode           = 0; // 0 : Significance - 1 : Yield

  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << " Making the calculation for channel " << chan <<  std::endl;
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << " Printing Region " << region << " for a luminosity of " << scaleFactor << " ipb"      << std::endl; 
  std::cout << " Total background uncertainty is " << realtiveBGError << std::endl; 
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
  std::cout << "==================================================================================================================================" << std::endl;
  std::cout << "|    Channel   |            ee            ||            mm             ||            em             ||           Total           |" << std::endl;
  std::cout << "==================================================================================================================================" << std::endl;

  // Background Part
  double background[6][4]     = {{0.}};
  double backgroundErr[6][4]  = {{0.}};
  unsigned int nProcess = 6;
  std::string processNames[6] = {"ttbar","singletop","VV","W","Z","Total"};
  unsigned int nChannel = 4;
  std::string channelNames[4] = {"EE","MM","EM","Total"};

  for(unsigned int process=0; process<nProcess-1; ++process) {
    std::cout << "|" << std::setw(10) << processNames[process] << std::setw(5) << "|";
    for(unsigned int channel=0; channel<nChannel-1; ++channel) {
      // Individual Channels
      background[process][channel] = scaleFactor*getPrediction(processNames[process],region,channelNames[channel],0);
      backgroundErr[process][channel] = scaleFactor*getPrediction(processNames[process],region,channelNames[channel],1);
      std::cout << std::setw(10) << std::setprecision(3) << background[process][channel]  << " +/- " 
                << std::setw(10) << std::setprecision(3) << backgroundErr[process][channel] << " || ";
      // Total Channel 
      background[process][nChannel-1] += background[process][channel];
      backgroundErr[process][nChannel-1] += pow(backgroundErr[process][channel],2);
    }
    backgroundErr[process][nChannel-1] = sqrt(backgroundErr[process][nChannel-1]); // Take sqrt(pow**2)
    std::cout << std::setw(10) << std::setprecision(3) << background[process][nChannel-1]  << " +/- " 
              << std::setw(10) << std::setprecision(3) << backgroundErr[process][nChannel-1] << " |" << std::endl;

  }

  std::cout << "==================================================================================================================================" << std::endl;
  std::cout << "|" << std::setw(10) << processNames[nProcess-1] << std::setw(5) << "|";
  for(unsigned int channel=0; channel<nChannel; ++channel) {
    for(unsigned int process=0; process<nProcess-1; ++process) {
      background[nProcess-1][channel]    += background[process][channel];
      backgroundErr[nProcess-1][channel] += pow(backgroundErr[process][channel],2);
    }
    backgroundErr[nProcess-1][channel] = sqrt(backgroundErr[nProcess-1][channel]);
    std::cout << std::setw(10) << std::setprecision(3) << background[nProcess-1][channel]  << " +/- " 
              << std::setw(10) << std::setprecision(3) << backgroundErr[nProcess-1][channel] << " || ";
  }
  std::cout << std::endl;
  std::cout << "==================================================================================================================================" << std::endl;

  // Signal Part
  for(unsigned int i=0; i<nSigPoints; i++) {

    double signal[4]    = {0.};
    double signalErr[4] = {0.};

    for(unsigned int channel= 0; channel<nChannel; ++channel) {
      if(channel<3) {
        signal[channel]    = scaleFactor*getPrediction(SignalDSIDs[i],region,channelNames[channel],0);
        signalErr[channel] = scaleFactor*getPrediction(SignalDSIDs[i],region,channelNames[channel],1);
      } else {
        for(unsigned int sum=0; sum<3; ++sum) {
          signal[3]    += signal[sum];
          signalErr[3] += pow(signalErr[sum],2);
        }
        signalErr[3] = sqrt(signalErr[3]);
      }
    }

    double signf = 0.;

    if(chan == "EE") {
      signf = RooStats::NumberCountingUtils::BinomialExpZ(signal[0],background[nProcess-1][0],realtiveBGError);
      yieldMap       [SignalDSIDs[i]] = signal[0];
    }
    else if(chan == "MM") {
      signf = RooStats::NumberCountingUtils::BinomialExpZ(signal[1],background[nProcess-1][1],realtiveBGError);
      yieldMap       [SignalDSIDs[i]] = signal[1];
    }
    else if(chan == "EM") {
      signf = RooStats::NumberCountingUtils::BinomialExpZ(signal[2],background[nProcess-1][2],realtiveBGError);
      yieldMap       [SignalDSIDs[i]] = signal[2];
    }
    else if(chan == "All") {
      signf  = pow(RooStats::NumberCountingUtils::BinomialExpZ(signal[0],background[nProcess-1][0],realtiveBGError),2);
      signf += pow(RooStats::NumberCountingUtils::BinomialExpZ(signal[1],background[nProcess-1][1],realtiveBGError),2);
      if(model == "ModeC") signf += pow(RooStats::NumberCountingUtils::BinomialExpZ(signal[2],background[nProcess-1][2],realtiveBGError),2); // Only for ModeC
      signf = sqrt(signf);
      yieldMap       [SignalDSIDs[i]] = signal[3];
    }

    std::cout << "|    "+SignalDSIDs[i]+"    |" << std::setw(10) << std::setprecision(3) << signal[0] << " +/- " << std::setw(10)  << std::setprecision(3) << signalErr[0]  << " || "
                                                << std::setw(10) << std::setprecision(3) << signal[1] << " +/- " << std::setw(10)  << std::setprecision(3) << signalErr[1]  << " || " 
                                                << std::setw(10) << std::setprecision(3) << signal[2] << " +/- " << std::setw(10)  << std::setprecision(3) << signalErr[2]  << " || "
                                                << std::setw(10) << std::setprecision(3) << signal[3] << " +/- " << std::setw(10)  << std::setprecision(3) << signalErr[3] << " || "
                                                << std::setw(10) << std::setprecision(3) << signf << " |" << std::endl;
    significanceMap[SignalDSIDs[i]] = signf; 
  }
  std::cout << "================================================================================================================================================" << std::endl;

  // Make the plot
  if(doPlot) makePlot(significanceMap,yieldMap,region,chan,plotMode,model);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Get Prediction
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getPrediction(TString source, TString region, TString channel, int mode)
{
 
  if(!m_inputFile->IsOpen()){
    std::cout << "Cannot find file, quitting..." << std::endl;
    return 0;
  }
  
  TTree* tree = (TTree*) m_inputFile->Get(source+"_CENTRAL");
  
  if(!tree){
    std::cout << "Cannot find tree, quitting..." << std::endl;
    return 0;
  }
  
  double integral = 0., error = 0.;
  TH1F* histo = new TH1F("histo","histo",1,0.5,1.5);
  histo->Sumw2();

  // Channel
  TString chCut = "";
  if(channel == "EE") chCut="l_flav[0]==0&&l_flav[1]==0";
  else if(channel == "MM") chCut="l_flav[0]==1&&l_flav[1]==1";
  else if(channel == "EM") chCut="l_flav[0]!=l_flav[1]";

  // OS
  TString osCut = "(l_q[0]*l_q[1])<0";

  // Lepton pT
  TString ptCut = "l_pt[0]>25.&&l_pt[1]>20.";

  // Jet-veto
  TString jetCut = "nCentralLJets==0&&nCentralBJets==0&&nForwardJets==0";

  // Mll
  TString mllCut = "mll>60.";
  if(channel=="EE"||channel=="MM") mllCut+="&&TMath::Abs(mll-90.2)>10."; // Z-veto

  // mT2
  TString mT2Cut = "";
  if(region.Contains("SRmT2")) {
    TObjArray* cutValues = region.Tokenize(","); 
    int cutValuesSize = cutValues->GetEntries();
    if(cutValuesSize == 2) { 
      TObjString *lower = (TObjString*) (*cutValues)[1]; 
      mT2Cut.Form("mT2lep>%s",lower->GetString().Data()); 
    }
    else if(cutValuesSize == 3) { 
      TObjString *lower = (TObjString*) (*cutValues)[1]; 
      TObjString *upper = (TObjString*) (*cutValues)[2]; 
      mT2Cut.Form("mT2lep>%s&&mT2lep<%s",lower->GetString().Data(),upper->GetString().Data());
    } 
  } else {
    std::cout << "Don't know region " << region << ", quitting ..." << std::endl;
    return 0;
  }

  // Full
  TString finalCut = "eventweight*(("     + chCut  + ")&&" 
                                    + "(" + osCut  + ")&&" 
                                    + "(" + ptCut  + ")&&" 
                                    + "(" + jetCut + ")&&" 
                                    + "(" + mllCut + ")&&" 
                                    + "(" + mT2Cut + "))";

  // Get integral and error
  tree->Draw("0.5>>histo",finalCut,"goff");
  integral = histo->IntegralAndError(0,-1,error);

  delete histo; histo=nullptr;
  delete tree; tree=nullptr;

  if(mode == 0)
    return integral;
  else if(mode == 1)
    return error;
  else{
    std::cout << "Unknown mode, returning zero..." << std::endl;
    return 0;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Plot Style
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void set_plot_style()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContourGraph
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraph* ContourGraph(TH2D* hist)
{
  TGraph* gr0  		= new TGraph();
  TH2D* histInternal  	= (TH2D*)  hist->Clone();
  //histInternal->Smooth(); // Smooth
  TGraph* gr1  		= (TGraph*) gr0->Clone();

  histInternal 		->SetContour( 1 );
  double pval  		= 0.05;
  double signif 	= TMath::NormQuantile(1-pval);
  histInternal		->SetContourLevel( 0, signif );
  histInternal		->Draw("CONT LIST");
  histInternal 		->SetDirectory(0);
  gPad 			->Update();

  TObjArray *contours 	= (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");

  TList *list 		= (TList*)contours->At(0);
  if( list->GetSize() == 0) return gr1;
  gr1 			= (TGraph*)list->First();

  return gr1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Get ModeC DSIDs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getModeCDSIDs(TString dsids[])
{
  TString SignalDSIDs[NMODECPOINTS] = { "392500","392501","392502","392504","392505","392506","392507","392508","392509", "392510",
                                        "392511","392512","392514","392515","392516","392517","392518","392519","392520" };

  std::copy(SignalDSIDs, SignalDSIDs+NMODECPOINTS,dsids);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Get DLiSlep DSIDs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getDLiSlepDSIDs(TString dsids[])
{
  TString SignalDSIDs[NDLISLEPPOINTS] = 
  {
   "111111" 
  };

  std::copy(SignalDSIDs, SignalDSIDs+NDLISLEPPOINTS,dsids);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Make Significance Plot
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void   makePlot( 
                 std::map<TString,double> significanceMap, 
                 std::map<TString,double> yieldMap, 
                 TString region, 
                 TString channel, 
                 int plotmode,
                 std::string model
               )
{
  // Set ATLAS Style
  SetAtlasStyle();

  // Print information if in debug
  for(std::map<TString,double>::iterator iter=significanceMap.begin(); iter!=significanceMap.end() && DEBUG; ++iter){
    std::cout << (*iter).first << " has a significance of " << std::setw(10)  << std::setprecision(3) << (*iter).second << std::endl;
  }

  // Read cross sections file
  std::string modelFile = "";
  if(model == "ModeC") modelFile = "/data/uclhc/uci/user/amete/grids/c1c1_slep.txt";
  else if(model == "DLiSlep") modelFile = "/export/home/amete/SimplifiedModelCrossSections/DLiSlep_Masses_Right.txt"; // Currently use RL

  std::ifstream inputFile(modelFile.c_str());

  // Variables
  int DSID; float mC1, /*mSl,*/ mN1/*, xsec_prod, uncert*/;
  std::string line;

  unsigned int nPoint = 0;
  double x[1000] = {0.}, y[1000] = {0.}, z[1000]= {0.};
  TString datasetId;

  if( inputFile.is_open() ){
    while( inputFile.good() ){
      getline (inputFile,line);
      try{
        sscanf(line.c_str(),"%i %f %f", &DSID,&mC1,&mN1);
      }
      catch(...){
        continue;
      }
      if( DSID > 100000 && DSID < 500000){
        x[nPoint] = mC1; 
        y[nPoint] = mN1;
        datasetId.Form("%i",DSID);
        if(plotmode == 0)
          z[nPoint] = significanceMap[datasetId];
        else if (plotmode == 1)
          z[nPoint] = yieldMap[datasetId];
        nPoint++;
      }
    }
  }

  // Hack Low Mass Issue
  //if( plotmode==0) { x[nPoint] = 112.5; y[nPoint] = 5; z[nPoint] = significanceMap["144902"]; nPoint++; }

  set_plot_style();
  gStyle->SetOptStat(kFALSE);

  // Create canvas
  TCanvas* c = new TCanvas( "c", "", 0, 0, 800, 800);  
  c->cd();

  // create and draw the frame
  int lowMassX  = 100.;
  int highMassX = 750.;
  int lowMassY  =   0.;
  int highMassY = 750.;

  if(model == "ModeC") {
    lowMassX  = 100.;
    highMassX = 700.;
    lowMassY  =   0.;
    highMassY = 500.;
  } else if(model == "DLiSlep") {
    lowMassX  = 100.;
    highMassX = 400.;
    lowMassY  =   0.;
    highMassY = 350.;
  }

  TH2F  *frame    = new TH2F("frame", "", 0.1*(highMassX-lowMassX),lowMassX, highMassX, 0.1*(highMassY-lowMassY), lowMassY, highMassY );
  TLine *diagonal = new TLine(lowMassX,lowMassX,200.,200.);
  diagonal->SetLineWidth(2);
  diagonal->SetLineStyle(2);

  TLatex *text     = new TLatex(0.36,0.53,"#Delta(m_{#tilde{#chi}_{1}^{#pm}} - m_{#tilde{#chi}_{1}^{0}}) = 0 ");
  text->SetNDC(kTRUE);
  text->SetTextSize(0.025);
  text->SetTextAngle(45);
 
  if(model == "ModeC") frame->SetXTitle( "m_{#tilde{#chi}_{1}^{#pm}} [GeV]" );
  else if(model == "DLiSlep") frame->SetXTitle( "m_{#tilde{l}_{L,R}^{#pm}} [GeV]" );
  frame->SetYTitle( "m_{#tilde{#chi}_{1}^{0}} [GeV]" );
  frame->GetYaxis()->SetTitleOffset(1.4);

  frame->GetXaxis()->SetTitleFont( 42 );
  frame->GetYaxis()->SetTitleFont( 42 );
  frame->GetZaxis()->SetTitleFont( 42 );
  frame->GetXaxis()->SetLabelFont( 42 );
  frame->GetYaxis()->SetLabelFont( 42 );
  frame->GetZaxis()->SetLabelFont( 42 );

  frame->GetXaxis()->SetTitleSize( 0.04 );
  frame->GetYaxis()->SetTitleSize( 0.04 );
  frame->GetZaxis()->SetTitleSize( 0.04 );
  frame->GetXaxis()->SetLabelSize( 0.04 );
  frame->GetYaxis()->SetLabelSize( 0.04 );
  frame->GetZaxis()->SetLabelSize( 0.04 );

  frame->GetXaxis()->SetNdivisions(505);

  // Markers
  TPolyMarker* gridpoints = new TPolyMarker(nPoint-1,x,y);
  //gridpoints            ->SetMarkerStyle(29);
  gridpoints            ->SetMarkerStyle(20);
  gridpoints            ->SetMarkerSize(1.5);
  gridpoints            ->SetMarkerColor(kBlack); 

  // Text
  if(plotmode==0)
    gStyle->SetPaintTextFormat(".1f");
  else if(plotmode==1)
    gStyle->SetPaintTextFormat(".0f");
  TH2D *numbers = new TH2D("numbers","numbers",(highMassX-lowMassX), lowMassX, highMassX, (highMassY-lowMassY), lowMassY, highMassY);
  for(unsigned int kk=0; kk < nPoint; ++kk) {
    //std::cout << "x " << x[kk] << " y " << y[kk] << " z " << z[kk] << std::endl;
    if(z[kk]>0)
      numbers->Fill(x[kk],y[kk],z[kk]);
    else
      numbers->Fill(x[kk],y[kk],0.001);
  }

  // TGraph2D
  TGraph2D* yield = new TGraph2D(nPoint-1,x,y,z);
  yield     ->GetZaxis()->SetTitleFont( 42 );
  yield     ->GetZaxis()->SetLabelFont( 42 );
  yield     ->GetZaxis()->SetTitleSize( 0.04 );
  yield     ->GetZaxis()->SetLabelSize( 0.04 );
  yield     ->GetZaxis()->SetTitleOffset( 1.5 );
  if(plotmode == 0)
    yield     ->GetZaxis()->SetTitle( "Zn" );
  else if(plotmode == 1)
    yield     ->GetZaxis()->SetTitle( "Signal Yield" );

  // Get the contour
  TH2D* yield2D     = (TH2D*) yield->GetHistogram()->Clone();
  TGraph* contGraph = (TGraph*) ContourGraph( yield2D )->Clone(); 
  contGraph->SetLineWidth(3);
  contGraph->SetLineColor(kRed+2);

  // Draw everything
  frame     ->Draw();
  diagonal  ->Draw();
  //text      ->Draw();
  yield     ->Draw("colz && same");
  //gridpoints->Draw("same && text");
  numbers   ->Draw("same && text");
  gPad      ->RedrawAxis();
  gPad      ->SetRightMargin(0.20);
  if(plotmode == 0) {
    contGraph->Draw("same");
    yield     ->GetZaxis()->SetRangeUser(0,3.5);
  }
  else if(plotmode == 1){
    yield     ->GetZaxis()->SetRangeUser(0,1.e2);
    gPad->SetLogz(1);
  }

  // Decoration
  TLatex *text1     = new TLatex(0.2,0.88,"#bf{#it{ATLAS}} Work In Progress");
  text1->SetNDC(kTRUE);
  text1->SetTextSize(0.04);
  text1->Draw();
  TLatex *text2     = new TLatex(0.2,0.82,"#scale[0.6]{#int} L dt = 10 fb^{-1}  #sqrt{s} = 13 TeV");
  text2->SetNDC(kTRUE);
  text2->SetTextSize(0.035);
  text2->Draw();
  TLatex *text3     = new TLatex(0.2,0.76,"Variable BG Uncertainty");
  text3->SetNDC(kTRUE);
  text3->SetTextSize(0.025);
  text3->Draw();
  TLatex *text4     = new TLatex(0.2,0.72,region+" - "+channel);
  text4->SetNDC(kTRUE);
  text4->SetTextSize(0.025);
  text4->Draw();
  TLatex *text5     = new TLatex(0.2,0.68,model.c_str());
  text5->SetNDC(kTRUE);
  text5->SetTextSize(0.025);
  text5->Draw();
  TLatex *text6     = new TLatex(0.2,0.64,"m_{#tilde{l}_{L},#tilde{#nu}} = (m_{#tilde{#chi}_{1}^{0}}+m_{#tilde{#chi}_{1}^{#pm}})/2");
  text6->SetNDC(kTRUE);
  text6->SetTextSize(0.025);
  if(model == "ModeC") text6->Draw();

  if(plotmode == 0)
    c->SaveAs("/data/uclhc/uci/user/amete/analysis_n0222_run/figures/powheg/10invfb_n0222_VariableUnc_"+region+"_"+channel+"_"+model+".eps");
  else if(plotmode == 1)
    c->SaveAs("/data/uclhc/uci/user/amete/analysis_n0222_run/figures/powheg/10invfb_n0222_VariableUnc_"+region+"_"+channel+"_"+model+"_Yield.eps");
}
