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
  if(m_inputROOTFile == NULL) {
    m_inputROOTFile = new TFile(m_inputFile);
  }
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
  int binWidth = ( atoi(m_binValuesList.at(2).c_str()) - atoi(m_binValuesList.at(1).c_str()) ) / atoi(m_binValuesList.at(0).c_str()); 

  // Get the central histograms 
  TH1D     *histograms[100];
  getHistogramsSimple(m_inputROOTFile, observable, cut, histograms,"CENTRAL");

  // Massage the Data and add MC to SM stack 
  TGraphAsymmErrors *Data = new TGraphAsymmErrors();
  THStack *mcStack        = new THStack("mcStack","Standard Model");
  int      dataIndex      = -1; 

  // Ad-hoc fix, normalize Zmumu to data - currently!!!
  histograms[0]->Scale(histograms[1]->Integral()/histograms[0]->Integral());

  for(unsigned int i=0; i<m_sampleList.size(); ++i) {
    // Add to stack if background
    if(m_sampleList.at(i)!="Data") {
      mcStack->Add(histograms[i]);
    }
    else {
      dataIndex = i;
      histograms[i]->SetMarkerSize(0);
      histograms[i]->SetLineWidth(0);
      histograms[i]->SetLineColor(0);
      // Convert data errors to Poisson
      convertErrorsToPoisson(histograms[i],Data);
      Data->SetMarkerSize(1.2);
      Data->SetMarkerStyle(20);
      Data->SetLineWidth(2);
    }
  }

  // Build the legend
  TLegend* legend        = new TLegend(0.7,0.6,0.9,0.9);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);

  // Add Data
  legend->AddEntry(Data,SampleNames["Data"],"p");

  // Add SM Background
  for(int i=m_sampleList.size()-1; i>-1 ; --i) {
    if(m_sampleList.at(i)!="Data")
      legend->AddEntry(histograms[i],SampleNames[m_sampleList.at(i)],"f");
  }

  ////////////////////////////////////////////////////////////////////////////////////////
  // HERE SHOULD COME THE SYSTEMATICS LOOP

  // Get Total MC Histo
  TH1D* stackHisto = (TH1D*) mcStack->GetStack()->Last();

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
      else if ( m_sampleList.at(j) != "Data" && m_sampleList.at(j) != "Fakes" ) {
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
  topPad->Draw();
  botPad->Draw();

  // Top Pad
  topPad               ->cd();
  topPad               ->SetBottomMargin(0.15);
  histograms[dataIndex]->Draw("p");
  mcStack              ->Draw("same && hists");
  nominalAsymErrors    ->Draw("same && E2");
  Data                 ->Draw("same && p");
  legend               ->Draw();

  // Set a few fancy labels and set axis ranges
  TString ylabel = "";
  ylabel.Form("Events /%i",binWidth);
  if(m_convertToGeV)
    ylabel.Append(" GeV");
  TString xlabel = "";

  if( variable.EqualTo("mT2") ) {
    xlabel = "m_{T2} [GeV]";
  }
  else if( variable.EqualTo("mll") ) {
    xlabel = "m_{ll} [GeV]";
  }
  else if( variable.EqualTo("mjj") ) {
    xlabel = "m_{jj} [GeV]";
  }
  else if( variable.EqualTo("metRel") ) {
    xlabel = "E_{T}^{miss,rel} [GeV]";
  }
  else if( variable.EqualTo("met") ) {
    xlabel = "E_{T}^{miss} [GeV]";
  }
  else if( variable.EqualTo("ptL0") ) {
    xlabel = "p_{T,l0} [GeV]";
  }
  else if( variable.EqualTo("ptL1") ) {
    xlabel = "p_{T,l1} [GeV]";
  }
  else if( variable.EqualTo("ptvarcone20L0") ) {
    xlabel = "p_{T,l0}^{varcone20} [GeV]";
  }
  else if( variable.EqualTo("ptvarcone20L1") ) {
    xlabel = "p_{T,l1}^{varcone20} [GeV]";
  }
  else if( variable.EqualTo("ptvarcone30L0") ) {
    xlabel = "p_{T,l0}^{varcone30} [GeV]";
  }
  else if( variable.EqualTo("ptvarcone30L1") ) {
    xlabel = "p_{T,l1}^{varcone30} [GeV]";
  }
  else if( variable.EqualTo("etaL0") ) {
    xlabel = "#eta_{l0}";
  }
  else if( variable.EqualTo("etaL1") ) {
    xlabel = "#eta_{l1}";
  }
  else if( variable.EqualTo("phiL0") ) {
    xlabel = "#phi_{l0}";
  }
  else if( variable.EqualTo("phiL1") ) {
    xlabel = "#phi_{l1}";
  }
  else if( variable.EqualTo("ptll") ) {
    xlabel = "p_{T,ll} [GeV]";
  }
  else {
    xlabel = variable;
  }
  histograms[dataIndex]->GetXaxis()->SetTitle(xlabel); 
  histograms[dataIndex]->GetXaxis()->SetLabelOffset(1.2); 
  histograms[dataIndex]->GetXaxis()->SetLabelSize(0.03);
  histograms[dataIndex]->GetYaxis()->SetTitle(ylabel); 
  if(m_plotLog)
    histograms[dataIndex]->GetYaxis()->SetRangeUser(2.e-2,1000*pow(10,ceil(log(histograms[dataIndex]->GetMaximum())/log(10))));
  else
    histograms[dataIndex]->GetYaxis()->SetRangeUser(histograms[dataIndex]->GetMinimum()*0.8,histograms[dataIndex]->GetMaximum()*1.20);

  gPad->RedrawAxis();
  if(m_plotLog)
    gPad->SetLogy(1);

  // Decoration
  char annoyingLabel1[100] = "#bf{#it{ATLAS}} Internal", annoyingLabel2[100] = "#scale[0.6]{#int} L dt = 6.6 pb^{-1}  #sqrt{s} = 13 TeV";
  myText(0.20,0.88,kBlack,annoyingLabel1);
  myText(0.20,0.80,kBlack,annoyingLabel2);

  // Bottom Pad
  botPad->cd();
  botPad->SetBottomMargin(0.3);

  // Dummy ratio histogram to set the scale and titles, etc.
  // I don't really like this part
  TH1D* ratio_original = (TH1D*) histograms[dataIndex]->Clone();
  ratio_original->Reset();
  ratio_original->SetMarkerSize(1.2);
  ratio_original->SetMarkerStyle(20);
  ratio_original->SetLineColor(kBlack);
  ratio_original->SetLineWidth(2);
  ratio_original->GetXaxis()->SetTitle(xlabel);
  ratio_original->GetYaxis()->SetTitle("Data/SM");
  ratio_original->GetXaxis()->SetLabelSize(0.1);
  ratio_original->GetXaxis()->SetLabelOffset(0.02);
  ratio_original->GetXaxis()->SetTitleSize(0.12);
  ratio_original->GetXaxis()->SetTitleOffset(1.);
  ratio_original->GetYaxis()->SetRangeUser(0,2);
  ratio_original->GetYaxis()->SetLabelSize(0.1);
  ratio_original->GetYaxis()->SetTitleSize(0.12);
  ratio_original->GetYaxis()->SetTitleOffset(0.5);
  ratio_original->GetYaxis()->SetNdivisions(5);

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
  TLine *line = new TLine(histograms[dataIndex]->GetXaxis()->GetXmin(),1,histograms[dataIndex]->GetXaxis()->GetXmax(),1);
  line->SetLineColor(kRed);
  line->SetLineStyle(7);

  // Draw
  ratio_original->Draw();
  ratioBand     ->Draw("same && E2");
  line          ->Draw();
  ratio         ->Draw("same && P");
  gPad          ->SetGridy(1);

  TString plotName = channel + "_" + region + "_" + variable + ".eps" ;
  //plotName = dirOut + "/" + plotName;
  canvas->SaveAs(plotName);

  // Delete unnecessary stuff to open up memory
  //delete[] histograms;
  //delete[] sysHistograms;
  delete Data;
  delete mcStack;
  delete nominalAsymErrors;
  delete nominalAsymErrorsNoError;
  delete ratioBand;
  delete ratio_original;
  delete ratio_raw;
  delete ratio;
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
void PlotMaker::convertErrorsToPoisson(TH1* histo, TGraphAsymmErrors* graph) 
{
  // Needed variables
  double value = 0;
  double error_poisson_up   = 0;
  double error_poisson_down = 0;
  double alpha = 0.158655, beta = 0.158655; // 68%

  // loop over bins and overwrite values
  for(int i=1; i<=histo->GetNbinsX(); i++){
    value = histo->GetBinContent(i);
    if(value!=0) {
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
      //tree->Draw( varToPlot + "/1000.>>temp" , "eventweight*(" + cutToApply + ")" );
      tree->Draw( varToPlot + "/1000.>>temp" , "(" + cutToApply + ")" );
    else
      //tree->Draw( varToPlot + ">>temp"       , "eventweight*(" + cutToApply + ")" );
      tree->Draw( varToPlot + ">>temp"       , "(" + cutToApply + ")" );

    // Clone and beautify
    histos[i] = (TH1D*) temp->Clone();
    histos[i]->SetName(histoName);
    histos[i]->SetTitle(histoName);
    histos[i]->SetLineWidth(2);
    histos[i]->SetLineColor(kBlack);
    histos[i]->SetFillColor(SampleColors[m_sampleList.at(i)]);

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
