//////////////////////////////////////////////////////////////////////////////////////////
// \project:    ATLAS Experiment at CERN's LHC
// \package:    SusyAnalysis
// \class:      PlotMaker
// \file:       $Id$
// \author:     Alaettin.Serhan.Mete@cern.ch
// \history:    N/A 
// 
// Copyright (C) 2015 University of California, Irvine
//////////////////////////////////////////////////////////////////////////////////////////

#include "SusyAnalysis/AtlasStyle.h"
#include "SusyAnalysis/AtlasLabels.h"
#include "SusyAnalysis/AtlasUtils.h"
#include "SusyAnalysis/PlotMaker.h"
#include "SusyAnalysis/PlotMakerDefs.h"
#include "RooStats/NumberCountingUtils.h"

#include "TTree.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"

#include <stdlib.h> 

/// \brief Constructor
PlotMaker::PlotMaker()
{
  ////////////////////////////////////////////////////////////////////////////////////////
  // Set ATLAS Style
  SetAtlasStyle();
}

/// \brief Constructor
PlotMaker::~PlotMaker() 
{

}

/// \brief Main Function that plots and saves histograms
void PlotMaker::generatePlot(TString channel, TString region, TString variable)
{
  float luminosity = 3209.05; // in pb-1
  //float luminosity = 10000.0; // in pb-1
  int   drawRatio  = 1; // 0 : no - 1 : data/mc - 2 : zbi
  bool  countAbove = true;
  bool  blindData  = true;
  float blindThreshold = 90.;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Print information
  cout << "PlotMaker::INFO   Now plotting" << 
          " channel "    << channel  <<
          " region "     << region   <<
          " variable "   << variable << endl;

  cout << "PlotMaker::INFO   Included samples     : ";
  for(unsigned int i=0; i<m_sampleList.size(); ++i) {                          
    cout << m_sampleList.at(i) << " ";
  }
  cout << endl;

  cout << "PlotMaker::INFO   Included systematics : ";
  for(unsigned int i=0; i<m_systematicsList.size(); ++i) {                          
    cout << m_systematicsList.at(i) << " ";
  }
  cout << endl;

  ////////////////////////////////////////////////////////////////////////////////////////
  // Open the input ROOT file

  //if(m_inputROOTFile == nullptr) {
    m_inputROOTFile = new TFile(m_inputFile);
 // }
  if(!m_inputROOTFile->IsOpen()) {
    cerr << "PlotMaker::ERROR   Cannot read ROOT file " << m_inputFile << endl;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Find the Region Cut
  TString cut = RegionCuts[region];
  if(cut.Length()==0) {
    cerr << "PlotMaker::ERROR   Unknown region " << region << endl;
    return;
  }

  // Append Channel to the Cut
  if(channel.EqualTo("ee"))
    cut.Append("&&l_flav[0]==0&&l_flav[1]==0");
  else if(channel.EqualTo("em"))
    cut.Append("&&((l_flav[0]==0&&l_flav[1]==1)||(l_flav[0]==1&&l_flav[1]==0))");
  else if(channel.EqualTo("mm"))
    cut.Append("&&l_flav[0]==1&&l_flav[1]==1");
  else if(channel.EqualTo("sf"))
    cut.Append("&&((l_flav[0]==1&&l_flav[1]==1)||(l_flav[0]==0&&l_flav[1]==0))");
  else if(!channel.EqualTo("all")) {
    cerr << "PlotMaker::ERROR   Unknown channel " << channel << endl;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Find the variable in the Ntuple
  TString observable = VariableNames[variable];
  if(observable.Length()==0) {
    cerr << "PlotMaker::ERROR   Unknown variable " << variable << endl;
    return;
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Read the bin values set by the user and intialize the temporary histogram
  if(m_binValuesList.size()!=3) {
    cerr << "PlotMaker::ERROR   Wrong Bin Value List Size " << m_binValuesList.size() << endl;
    return;
  }
  double binWidth = ( atof(m_binValuesList.at(2).c_str()) - atof(m_binValuesList.at(1).c_str()) ) / atof(m_binValuesList.at(0).c_str()); 

  // Get the central histograms 
  TH1D     *histograms[100];
  getHistogramsSimple(m_inputROOTFile, observable, cut, histograms,"CENTRAL");

  // Massage the Data and add MC to SM stack 
  TGraphAsymmErrors *Data = new TGraphAsymmErrors();
  THStack *mcStack        = new THStack("mcStack","Standard Model");
  int      dataIndex      = -1; 
  std::vector<int> signalIndices;

  for(unsigned int i=0; i<m_sampleList.size(); ++i) {
    // Add to stack if background
    if(m_sampleList.at(i)!="Data") {
      if(m_sampleList.at(i)!="MM") histograms[i]->Scale(luminosity);
      if(m_sampleList.at(i).find("406")==std::string::npos &&
         m_sampleList.at(i).find("392")==std::string::npos ) // Don't add signal to stack
        mcStack->Add(histograms[i]);
      else
        signalIndices.push_back(i);
    }
    else {
      dataIndex = i;
      histograms[i]->SetLineWidth(0);
      histograms[i]->SetLineColor(0);
      // Convert data errors to Poisson
      convertErrorsToPoisson(histograms[i],Data,blindThreshold,blindData);
      Data->SetMarkerSize(1.2);
      Data->SetMarkerStyle(20);
      Data->SetLineWidth(2);
    }
  }

  // Build the legend
  TLegend* legend        = new TLegend(0.7,0.6,0.9,0.92);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  // Add Data
  if(dataIndex!=-1) legend->AddEntry(Data,SampleNames["Data"],"p");

  // Add SM Background and signal
  for(int i=m_sampleList.size()-1; i>-1 ; --i) {
    if(m_sampleList.at(i)!="Data")
    {
      if(m_sampleList.at(i).find("406")==std::string::npos ||
         m_sampleList.at(i).find("392")==std::string::npos)
        legend->AddEntry(histograms[i],SampleNames[m_sampleList.at(i)],"f");
      else
        legend->AddEntry(histograms[i],SampleNames[m_sampleList.at(i)],"l");
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // HERE SHOULD COME THE SYSTEMATICS LOOP

  // Get Total MC Histo
  TH1D* stackHisto = (TH1D*) mcStack->GetStack()->Last();
  double bkgError = 0., bkgTot = stackHisto->IntegralAndError(0,-1,bkgError);
  TH1D* dummyHisto = (TH1D*) histograms[0]->Clone();
  dummyHisto->Reset();

  TGraphAsymmErrors* nominalAsymErrors = TH1TOTGraph(stackHisto);
  nominalAsymErrors->SetMarkerSize(0);
  nominalAsymErrors->SetLineWidth(2);
  nominalAsymErrors->SetFillStyle(3004);
  nominalAsymErrors->SetFillColor(kBlack);

  // Add SM Bkg. to the legend
  legend->AddEntry(nominalAsymErrors,"Bkg. Uncert.","f"); 

  // Loop over and add the Tree based systematics
  TH1D* totalSysHisto = (TH1D*) stackHisto->Clone();
  totalSysHisto->Reset();
  TGraphAsymmErrors* transient; 
  TH1D    *sysHistograms[100];

  for(unsigned int i=0; i < m_systematicsList.size(); ++i) {

    // Retrieve the histograms
    getHistogramsSimple(m_inputROOTFile, observable, cut, sysHistograms, m_systematicsList.at(i));

    // Loop over samples and add to total systematics histo
    for(unsigned int j=0; j < m_sampleList.size(); ++j) {

      if(sysHistograms[j]!=NULL)
        totalSysHisto->Add(sysHistograms[j]);
      else if ( m_sampleList.at(j) != "Data" && m_sampleList.at(j) != "MM"/*"Fakes"*/ ) {
        cout << "PlotMaker::WARNING   Cannot find TTree for systematics " << 
                m_systematicsList.at(i) << " for sample " << m_sampleList.at(j)  << endl;
      }
    }

    // Add to the band
    transient = TH1TOTGraph(totalSysHisto);
    myAddtoBand(transient,nominalAsymErrors);
    totalSysHisto->Reset();
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // Draw the canvas
  TCanvas* canvas = new TCanvas("canvas","canvas",500,500);
  TPad*    topPad = new TPad("pTop","pTop",0,0.2,1,1);
  TPad*    botPad = new TPad("pBot","pBot",0,0.0,1,0.3);
  if(drawRatio) {
    topPad->Draw();
    botPad->Draw();
  }

  // Top Pad
  if(drawRatio) {
    topPad               ->cd();
    topPad               ->SetBottomMargin(0.15);
  }
  dummyHisto           ->Draw();
  mcStack              ->Draw("same && hists");
  nominalAsymErrors    ->Draw("same && E2");
  if(dataIndex!=-1)    Data->Draw("same && p");
  for(auto &isig : signalIndices) {
    histograms[isig]   ->Draw("same && hists");
  }
  legend               ->Draw();

  // Set a few fancy labels and set axis ranges
  TString ylabel = "Events";
  ylabel.Form("Events /%.2f",binWidth);
  TString xlabel = "";

  if( variable.EqualTo("ptL0") ) {
    xlabel = "p_{T,l0} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptL1") ) {
    xlabel = "p_{T,l1} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptJ0") ) {
    xlabel = "p_{T,j0} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptJ1") ) {
    xlabel = "p_{T,j1} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptvarcone20L0") ) {
    xlabel = "p_{T,l0}^{varcone20} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptvarcone20L1") ) {
    xlabel = "p_{T,l1}^{varcone20} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptvarcone30L0") ) {
    xlabel = "p_{T,l0}^{varcone30} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptvarcone30L1") ) {
    xlabel = "p_{T,l1}^{varcone30} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("etaL0") ) {
    xlabel = "#eta_{l0}";
  }
  else if( variable.EqualTo("etaL1") ) {
    xlabel = "#eta_{l1}";
  }
  else if( variable.EqualTo("etaJ0") ) {
    xlabel = "#eta_{j0}";
  }
  else if( variable.EqualTo("etaJ1") ) {
    xlabel = "#eta_{j1}";
  }
  else if( variable.EqualTo("phiL0") ) {
    xlabel = "#phi_{l0}";
  }
  else if( variable.EqualTo("phiL1") ) {
    xlabel = "#phi_{l1}";
  }
  else if( variable.EqualTo("phiJ0") ) {
    xlabel = "#phi_{j0}";
  }
  else if( variable.EqualTo("phiJ1") ) {
    xlabel = "#phi_{j1}";
  }
  else if( variable.EqualTo("mT2") ) {
    xlabel = "m_{T2} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("mll") ) {
    xlabel = "m_{ll} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("ptll") ) {
    xlabel = "p_{T,ll} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("drll") ) {
    xlabel = "#DeltaR_{ll}";
  }
  else if( variable.EqualTo("dphill") ) {
    xlabel = "#Delta#phi_{ll}";
  }
  else if( variable.EqualTo("metRel") ) {
    xlabel = "E_{T}^{miss,rel} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("met") ) {
    xlabel = "E_{T}^{miss} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("mT2lep") ) {
    xlabel = "m_{T2}(ll) [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("meff") ) {
    xlabel = "m_{eff} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("abs_deltaX") ) {
    xlabel = "|x_{1}-x_{2}|";
  }
  else if( variable.EqualTo("abs_cthllb") ) {
    xlabel = "|cos#theta_{b}|";
  }
  else if( variable.EqualTo("MDR_jigsaw") ) {
    xlabel = "m_{#Delta}^{R} [GeV]";
    ylabel.Append(" GeV");
  }
  else if( variable.EqualTo("RPT_jigsaw") ) {
    xlabel = "RPT";
  }
  else if( variable.EqualTo("gamInvRp1_jigsaw") ) {
    xlabel = "1/#gamma_{P}^{PP}";
  }
  else if( variable.EqualTo("DPB_vSS_jigsaw") ) {
    xlabel = "#Delta#phi_{#beta}^{R}";
  }
  else {
    xlabel = variable;
  }
  dummyHisto->GetXaxis()->SetTitle(xlabel); 
  dummyHisto->GetXaxis()->SetTitleSize(0.05);
  if(drawRatio) { 
    dummyHisto->GetXaxis()->SetLabelOffset(1.2); 
    dummyHisto->GetXaxis()->SetLabelSize(0.03);
  }
  else {
    dummyHisto->GetXaxis()->SetLabelSize(0.04);
  }
  dummyHisto->GetYaxis()->SetTitle(ylabel); 
  dummyHisto->GetYaxis()->SetTitleSize(0.05);
  dummyHisto->GetYaxis()->SetTitleOffset(1.2);
  dummyHisto->GetYaxis()->SetLabelSize(0.04);
  if(m_plotLog)
    dummyHisto->GetYaxis()->SetRangeUser(1.e-2,1000*pow(10,ceil(log(stackHisto->GetMaximum())/log(10))));
  else
    dummyHisto->GetYaxis()->SetRangeUser(stackHisto->GetMinimum()*0.8,stackHisto->GetMaximum()*1.20);

  gPad->RedrawAxis();
  if(m_plotLog)
    gPad->SetLogy(1);

  // Decoration
  char annoyingLabel1[100] = "#bf{#it{ATLAS}} Internal";
  TString luminosityInfo; luminosityInfo.Form("#sqrt{s} = 13 TeV, %.2f fb^{-1}",luminosity*0.001); 
  char* annoyingLabel2 = (char*) luminosityInfo.Data();
  char* annoyingLabel3 = (char*) region.Data();
  char* annoyingLabel4 = (char*) channel.Data();
  myText(0.20,0.88,kBlack,annoyingLabel1);
  myText(0.20,0.80,kBlack,annoyingLabel2);
  myText(0.20,0.72,kBlack,annoyingLabel3);
  myText(0.50,0.72,kBlack,annoyingLabel4);

  // Bottom Pad
  if(drawRatio) {
    botPad->cd();
    botPad->SetBottomMargin(0.3);

    // Dummy ratio histogram to set the scale and titles, etc.
    // I don't really like this part
    TH1D* ratio_original = (TH1D*) dummyHisto->Clone(); //(TH1D*) stackHisto->Clone();
    ratio_original->Reset();
    ratio_original->SetMarkerSize(1.2);
    ratio_original->SetMarkerStyle(20);
    ratio_original->SetLineColor(kBlack);
    ratio_original->SetLineWidth(2);
    ratio_original->GetXaxis()->SetTitle(xlabel);
    if(drawRatio == 1) ratio_original->GetYaxis()->SetTitle("Data/SM");
    else if(drawRatio == 2) ratio_original->GetYaxis()->SetTitle("Z_{Bi} (High-#Delta M)");
    ratio_original->GetXaxis()->SetLabelSize(0.1);
    ratio_original->GetXaxis()->SetLabelOffset(0.02);
    ratio_original->GetXaxis()->SetTitleSize(0.12);
    ratio_original->GetXaxis()->SetTitleOffset(1.);
    if(drawRatio == 1) ratio_original->GetYaxis()->SetRangeUser(0,2);
    else if(drawRatio == 2) ratio_original->GetYaxis()->SetRangeUser(0,6);
    ratio_original->GetYaxis()->SetLabelSize(0.1);
    ratio_original->GetYaxis()->SetTitleSize(0.12);
    ratio_original->GetYaxis()->SetTitleOffset(0.5);
    ratio_original->GetYaxis()->SetNdivisions(5);

    // Data/MC
    if(drawRatio == 1) {
      // Build the ratio's error band
      TGraphAsymmErrors* ratioBand   = new TGraphAsymmErrors( *nominalAsymErrors ); 
      buildRatioErrorBand(nominalAsymErrors,ratioBand);
    
      // Get the ratio
      // For Data/MC only use the statistical error for data
      // because we explicitly draw the MC error band
      TGraphAsymmErrors* nominalAsymErrorsNoError = new TGraphAsymmErrors( *nominalAsymErrors );
      for(int i=1; i<=nominalAsymErrorsNoError->GetN(); ++i) {
        nominalAsymErrorsNoError->SetPointError(i-1,0,0,0,0);
      }
      TGraphAsymmErrors* ratio_raw = myTGraphErrorsDivide(Data,nominalAsymErrorsNoError);
      TGraphAsymmErrors* ratio     = new TGraphAsymmErrors();
    
      double x1=0; double y1=0; unsigned int newIndex = 0.;
      for(int kk=0; kk<ratio_raw->GetN(); ++kk){
        ratio_raw->GetPoint(kk, x1,y1);
        if(y1 > 0.) {
          ratio->SetPoint(newIndex, x1, y1);
          ratio->SetPointError(newIndex, ratio_raw->GetErrorXlow(kk), ratio_raw->GetErrorXhigh(kk), ratio_raw->GetErrorYlow(kk), ratio_raw->GetErrorYhigh(kk));
          newIndex++;
        }
      }
      ratio->SetMarkerSize(1.2);
      ratio->SetMarkerStyle(20);
      ratio->SetLineColor(kBlack);
      ratio->SetLineWidth(2);
    
      // Make the line at 1 for the ratio plot
      TLine *line = new TLine(dummyHisto->GetXaxis()->GetXmin(),1,dummyHisto->GetXaxis()->GetXmax(),1);
      line->SetLineColor(kRed);
      line->SetLineStyle(7);
    
      // Draw
      ratio_original->Draw();
      ratioBand     ->Draw("same && E2");
      line          ->Draw();
      ratio         ->Draw("same && P && 0");
      gPad          ->SetGridy(1);
    } // drawRatio == 1
    else if (drawRatio == 2) {
      TGraphAsymmErrors* ratio     = new TGraphAsymmErrors();

      for(int kk=1; kk<stackHisto->GetXaxis()->GetNbins()+1; ++kk){
        float bgFromThreshold      = countAbove ? stackHisto->Integral(kk,-1) : stackHisto->Integral(0,kk);
        float xlowEdgeAtThreshold  = stackHisto->GetXaxis()->GetBinLowEdge(kk);
        float xhighEdgeAtThreshold = xlowEdgeAtThreshold + stackHisto->GetXaxis()->GetBinWidth(kk);
        float xcenterAtThreshold   = stackHisto->GetXaxis()->GetBinCenter(kk);

        // Currently only one signal
        float sigFromThrehold    = 0.;
        for(unsigned int isample=0; isample<m_sampleList.size(); ++isample) {
          if( m_sampleList.at(isample).find("392510"/*"392508"*/)!=std::string::npos ) {
            sigFromThrehold = countAbove ? histograms[isample]->Integral(kk,-1) : histograms[isample]->Integral(0,kk);
            break;
          }
        }

        // Calculate and print
        float bgunc = 0.3;
        //if(xlowEdgeAtThreshold > 149.)       { bgunc = 0.50; }
        //else if(xlowEdgeAtThreshold > 119.)  { bgunc = 0.25; }
        //else                                 { bgunc = 0.15; }
        float signfFromThreshold = RooStats::NumberCountingUtils::BinomialExpZ( sigFromThrehold, bgFromThreshold, bgunc);
        TString printline;
        if (countAbove) printline.Form("x [%*.2f,inf] : BG = %.2f (%% %.2f) \t SIG = %.2f \t ZBi %.2f", 3, xlowEdgeAtThreshold, bgFromThreshold, bgunc, sigFromThrehold, signfFromThreshold);
        else            printline.Form("x [0,%*.2f] : BG = %.2f (%% %.2f) \t SIG = %.2f \t ZBi %.2f", 3, xhighEdgeAtThreshold, bgFromThreshold, bgunc, sigFromThrehold, signfFromThreshold);
        std::cout << printline << std::endl;

        // Make the line at 1 for the ratio plot
        TLine *line1 = new TLine(dummyHisto->GetXaxis()->GetXmin(),1.64,dummyHisto->GetXaxis()->GetXmax(),1.64);
        line1->SetLineColor(kRed);
        line1->SetLineStyle(7);
        TLine *line2 = new TLine(dummyHisto->GetXaxis()->GetXmin(),3.00,dummyHisto->GetXaxis()->GetXmax(),3.00);
        line2->SetLineColor(kGreen+3);
        line2->SetLineStyle(7);
        TLine *line3 = new TLine(dummyHisto->GetXaxis()->GetXmin(),5.00,dummyHisto->GetXaxis()->GetXmax(),5.00);
        line3->SetLineColor(kOrange+3);
        line3->SetLineStyle(7);
    
        // Set and draw
        ratio->SetPoint(kk-1, xcenterAtThreshold, signfFromThreshold);
        ratio_original->Draw();
        ratio         ->Draw("same && P && 0");
        line1         ->Draw();
        line2         ->Draw();
        line3         ->Draw();
        gPad          ->SetGridy(1);
      } // drawRatio == 2
    }
  } // end of drawRatio

  TString plotName = channel + "_" + region + "_" + variable + ".eps" ;
  //plotName = dirOut + "/" + plotName;
  canvas->SaveAs(plotName);

  // Print yields
  for(unsigned int isample=0; isample<m_sampleList.size(); ++isample) {
    TString line;
    double totErr = 0., tot = histograms[isample]->IntegralAndError(0,-1,totErr);
    if( m_sampleList.at(isample).find("406")!=std::string::npos ||
        m_sampleList.at(isample).find("392")!=std::string::npos ) {
      double signf   = RooStats::NumberCountingUtils::BinomialExpZ( tot, bkgTot, 0.3 );
      double signf_a = RooStats::NumberCountingUtils::BinomialExpZ( tot*(10/3.21), bkgTot*(10/3.21), 0.3 );
      double signf_b = RooStats::NumberCountingUtils::BinomialExpZ( tot*(10/3.21), bkgTot*(10/3.21), 0.5 );
      line.Form("%9s \t %.2f +/- %.2f \t Zbi (actual) = %.2f - Zbi (*10/3.21) = %.2f - Zbi (*10/3.21 and 50%%) = %.2f", m_sampleList.at(isample).c_str(),tot,totErr,signf,signf_a,signf_b);
      std::cout << line << std::endl;
    } else {
      line.Form("%9s \t %.2f +/- %.2f", m_sampleList.at(isample).c_str(),tot,totErr);
      std::cout << line << std::endl;
    }
    if(isample==m_sampleList.size()-1) {
      line.Form("Total BG \t %.2f +/- %.2f",bkgTot,bkgError);
      std::cout << line << std::endl;
    }
  }
  

  // Delete unnecessary stuff to open up memory
  //delete[] histograms;
  //delete[] sysHistograms;
  delete Data;
  delete mcStack;
  delete nominalAsymErrors;
  //delete nominalAsymErrorsNoError;
  //delete ratioBand;
  delete dummyHisto;
  //delete ratio_raw;
  //delete ratio;
  delete totalSysHisto;
  delete legend;
  delete canvas;

  return;
}

/// \brief Function to add overflow to last bin
void PlotMaker::addOverFlowToLastBin(TH1* histo) {

  // Find last bin
  int lastBin          = histo->GetXaxis()->GetNbins();

  // Read values
  double lastBinValue  = histo->GetBinContent(lastBin);
  double lastBinError  = histo->GetBinError(lastBin);
  double overFlowValue = histo->GetBinContent(lastBin+1); 
  double overFlowError = histo->GetBinError(lastBin+1);

  // Set Values    
  histo->SetBinContent(lastBin+1,0);
  histo->SetBinError(lastBin+1,0);
  histo->SetBinContent(lastBin,lastBinValue+overFlowValue);
  histo->SetBinError(lastBin,sqrt(lastBinError*lastBinError+overFlowError*overFlowError));

}

/// \brief Function to add overflow to last bin
void PlotMaker::convertErrorsToPoisson(TH1* histo, TGraphAsymmErrors* graph, float blindThreshold = 90., bool blindData = false) 
{
  // Needed variables
  double value = 0;
  double error_poisson_up   = 0;
  double error_poisson_down = 0;
  double alpha = 0.158655, beta = 0.158655; // 68%

  // loop over bins and overwrite values
  for(int i=1; i<=histo->GetNbinsX(); i++){
    value = histo->GetBinContent(i);
    if(value!=0 && ((blindData && (histo->GetBinLowEdge(i)+histo->GetBinWidth(i))<=blindThreshold)||!blindData)) {
      error_poisson_up     = 0.5*TMath::ChisquareQuantile(1-beta,2*(value+1)) - value; 
      error_poisson_down   = value - 0.5*TMath::ChisquareQuantile(alpha,2*value);
      graph->SetPoint(i-1, histo->GetBinCenter(i), value);
      graph->SetPointError(i-1, 0., 0. , error_poisson_down, error_poisson_up);
    }
    else{
      graph->SetPoint(i-1, histo->GetBinCenter(i), 0.);
      graph->SetPointError(i-1, 0., 0., 0., 0.);
    }
  }
}

/// \brief Function to get higstograms
void PlotMaker::getHistogramsSimple(TFile* input, TString varToPlot, TString cutToApply, TH1D* histos[], TString variation)
{
  // Initialize the temp histogram
  TH1D *temp = new TH1D("temp","temp",atof(m_binValuesList.at(0).c_str()),
                                      atof(m_binValuesList.at(1).c_str()),
                                      atof(m_binValuesList.at(2).c_str())); 
  temp->Sumw2();

  // Here is the loop over samples
  for(unsigned int i=0; i<m_sampleList.size(); ++i) {

    // Return NULL if systematics and Data or Fakes
    if( (m_sampleList.at(i) == "Data" || m_sampleList.at(i) == "Fakes") && !variation.EqualTo("CENTRAL") ) {
      histos[i] = NULL;
      continue;
    }

    // Define variables internal to the loop
    TString treeName = ""; treeName.Form("%s_",m_sampleList.at(i).c_str()); treeName += variation;
    TString histoName = ""; histoName.Form("histo_%s",m_sampleList.at(i).c_str());

    // Get the TTree - continue the loop if not found 
    TTree* tree = (TTree*) input->Get(treeName);
    if(tree == NULL) {
      cout << "PlotMaker::WARNING   Cannot find TTree w/ name " << treeName << endl;
      histos[i] = NULL;
      continue;
    }

    // Fill the temp histogram
    if(m_convertToGeV)
      tree->Draw( varToPlot + "/1000.>>temp" , "eventweight*(" + cutToApply + ")" );
    else
      tree->Draw( varToPlot + ">>temp"       , "eventweight*(" + cutToApply + ")" );

    // Clone and beautify
    histos[i] = (TH1D*) temp->Clone();
    histos[i]->SetName(histoName);
    histos[i]->SetTitle(histoName);
    if(m_sampleList.at(i).find("406")==std::string::npos &&
       m_sampleList.at(i).find("392")==std::string::npos) {
      histos[i]->SetLineWidth(2);
      histos[i]->SetLineColor(kBlack);
      histos[i]->SetFillColor(SampleColors[m_sampleList.at(i)]);
    } else {
      histos[i]->SetLineWidth(2);
      histos[i]->SetLineColor(SampleColors[m_sampleList.at(i)]);
      histos[i]->SetFillColor(0);
      histos[i]->SetMarkerSize(0);
    }

    // Add the overflow to the last bin
    addOverFlowToLastBin(histos[i]);

    delete tree;
  } // end of loop over samples

  temp->Reset();
  delete temp;
}

/// \brief Function to build the ratio's error band
void PlotMaker::buildRatioErrorBand(TGraphAsymmErrors* input, TGraphAsymmErrors* output) 
{
  output->SetMarkerSize(0);
  for(int bin=0; bin < output->GetN(); bin++){ 
    output->GetY()[bin] = 1.; 
    if( input->GetY()[bin] > 0.0001 ) 
      output->GetEYhigh()[bin]=input->GetEYhigh()[bin]/input->GetY()[bin]; 
   else 
      output->GetEYhigh()[bin]= 0.; 
   if( input->GetY()[bin] > 0.0001 ) 
      output->GetEYlow()[bin]=input->GetEYlow()[bin]/input->GetY()[bin]; 
   else 
      output->GetEYlow()[bin]= 0.; 
   if( output->GetEYlow()[bin] > 1. ) 
     output->GetEYlow()[bin] = 1.; 
   if( output->GetEYhigh()[bin] > 1. ) 
     output->GetEYhigh()[bin] = 1.; 
  }
}
