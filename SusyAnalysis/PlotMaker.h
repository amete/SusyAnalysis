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

#ifndef PLOTMAKER_H
#define PLOTMAKER_H

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"

using namespace std;

class PlotMaker {

  // Public
  public:

    /// \brief Constructor
    PlotMaker();
    /// \brief Constructor
    ~PlotMaker();

    /// \brief Set-Get Functions for input ROOT file
    void           setInputFile(TString arg) { m_inputFile = arg ; }; 
    TString        getInputFile(           ) { return m_inputFile; };

    /// \brief Set-Get Functions for sample list 
    void           setSampleList(vector<string> arg) { m_sampleList = arg;  };
    vector<string> getSampleList(                  ) { return m_sampleList; };
 
    /// \brief Set-Get Functions for systematics list
    void           setSystematicsList(vector<string> arg) { m_systematicsList = arg;  };
    vector<string> getSystematicsList(                  ) { return m_systematicsList; };

    /// \brief Set-Get Functions for bin values list
    void           setBinValuesList(vector<string> arg) { m_binValuesList = arg;  };
    vector<string> getBinValuesList(                  ) { return m_binValuesList; };

    /// \brief Set-Get Functions for GeV Flag
    void           setGeVFlag(bool arg) { m_convertToGeV = arg;  };
    bool           getGeVFlag(        ) { return m_convertToGeV; };         

    /// \brief Set-Get Functions for Log Flag
    void           setLogFlag(bool arg) { m_plotLog = arg;  };
    bool           getLogFlag(        ) { return m_plotLog; };         

    /// \brief Main Function that plots and saves histograms
    void           generatePlot(TString channel,TString region,TString variable);

    /// \brief Function to add overflow to last bin
    void           addOverFlowToLastBin(TH1* arg);

    /// \brief Function to get central histograms
    void           getHistogramsSimple(TFile* inputRootFile, TString variable, TString cut, TH1D* histos[], TString variation);

    /// \brief Function to add overflow to last bin
    void           convertErrorsToPoisson(TH1* inputHisto, TGraphAsymmErrors* outputGraph);

    /// \brief Function to build the ratio's error band
    void           buildRatioErrorBand(TGraphAsymmErrors* input, TGraphAsymmErrors* output);

  // Private
  private:

    TFile*         m_inputROOTFile;
    TString        m_inputFile;
    vector<string> m_sampleList;
    vector<string> m_systematicsList;
    vector<string> m_binValuesList;
    bool           m_convertToGeV;
    bool           m_plotLog;
};

#endif 
