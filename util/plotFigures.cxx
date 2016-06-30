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

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

#include "SusyAnalysis/PlotMaker.h"

using namespace std;

void tokenizeString(string input, vector<string> &result);
void help();

/// \brief Main Function
int main(int argc, char** argv)
{
  // Read user inputs
  string channels     = "mm"                                ; vector<string> channelList   ;
  string regions      = "CR2LOS"                            ; vector<string> regionList    ;
  string variables    = "l_pt"                              ; vector<string> variableList  ;
  //string samples      = "392508,392510,W,singletop,ttbar,VV,Z,Data"; vector<string> sampleList    ;
  //string samples      = "392508,392510,MM,singletop,ttbar,VV,Z,Data"; vector<string> sampleList    ;
  string samples      = "c1c1_slep_400.0_200.0,c1c1_slep_500.0_1.0,MM,singletop,ttbar,VV,Z,Data"; vector<string> sampleList    ;
  //string samples      = "392508,392510,W,singletop,ttbar,VV,Z"; vector<string> sampleList    ;
  //string samples      = "392508,392510,W,ttv,singletop,ttbar,VVV,VV,Z,Data"; vector<string> sampleList    ;
  string binValues    = "40,0,400"                          ; vector<string> binValueList  ;
  //string systematics  = "JER,JES_DN,JES_UP,EER_DN,EER_UP,EES_LOW_DN,EES_LOW_UP,EES_MAT_DN,EES_MAT_UP,EES_PS_DN,EES_PS_UP,EES_Z_DN,EES_Z_UP,ID_DN,ID_UP,MS_DN,MS_UP,RESOST,SCALEST_DN,SCALEST_UP"; 
  string systematics  = ""; 
  vector<string> systematicList;
  bool   convertToGeV = false;
  bool   plotLog = false;

  for(int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--channels") == 0)
      channels     = argv[++i];
    else if (strcmp(argv[i], "--regions") == 0)
      regions      = argv[++i];
    else if (strcmp(argv[i], "--variables") == 0)
      variables    = argv[++i];
    else if (strcmp(argv[i], "--systematics") == 0)
      systematics  = argv[++i];
    else if (strcmp(argv[i], "--samples") == 0)
      samples      = argv[++i];
    else if (strcmp(argv[i], "--binValues") == 0)
      binValues    = argv[++i];
    else if (strcmp(argv[i], "--isGeV") == 0)
      convertToGeV = true;
    else if (strcmp(argv[i], "--plotLog") == 0)
      plotLog = true;
    else {
      help();
      return 0; 
    }
  } // end of user inputs loop

  // Print
  cout << "=======================================================" << endl;
  cout << "Running PlotMaker w/ arguments:"                         << endl;
  cout << "=======================================================" << endl;
  cout << "Channels        = " << channels     << endl;
  cout << "Regions         = " << regions      << endl;
  cout << "Variables       = " << variables    << endl;
  cout << "Systematics     = " << systematics  << endl;
  cout << "Samples         = " << samples      << endl;
  cout << "Bin Values      = " << binValues    << endl;
  cout << "Convert to GeV  = " << convertToGeV << endl;
  cout << "=======================================================" << endl;

  // Tokenize inputs
  tokenizeString(channels   ,channelList   );
  tokenizeString(regions    ,regionList    );
  tokenizeString(variables  ,variableList  );
  tokenizeString(systematics,systematicList);
  tokenizeString(samples    ,sampleList    );
  tokenizeString(binValues  ,binValueList  );

  // Plot and Save
  PlotMaker* plots = new PlotMaker();
  plots->setInputFile("/data/uclhc/uci/user/amete/analysis_n0224_run/EWK2L/hfts_6_skimmed/HFT_COMBINED_13TeV.root"); // Rel 20.7
  //plots->setInputFile("/data/uclhc/uci/user/amete/analysis_n0224_run/EWK2L/hfts_4/mc15_13TeV.root"); // Rel 20.7
  //plots->setInputFile("/data/uclhc/uci/user/amete/analysis_n0222_run/hfts/mc15_13TeV_sherpa.root"); // Rel 20.1
  plots->setSampleList     (sampleList    );
  plots->setSystematicsList(systematicList);
  plots->setBinValuesList  (binValueList  );
  plots->setGeVFlag        (convertToGeV  );
  plots->setLogFlag        (plotLog       );

  // Loop over all combinations
  for(unsigned ii = 0; ii < channelList.size(); ++ii) {
    for(unsigned jj = 0; jj < regionList.size(); ++jj) {
      for(unsigned kk = 0; kk < variableList.size(); ++kk) {

        TString ch = ""; ch.Form("%s",channelList.at(ii).c_str() );
        TString rg = ""; rg.Form("%s",regionList.at(jj).c_str()  );
        TString vr = ""; vr.Form("%s",variableList.at(kk).c_str());

        plots->generatePlot(ch,rg,vr); 

      } // end of variable loop
    } // end of region loop
  } // end of channel loop

  return 0;
}

/// \brief Tokinize the string
void tokenizeString(string input, vector<string> &result) {

  for(/*unsigned int*/ size_t i=0,n; i<input.length(); i=n+1)
  {
    n = input.find_first_of(',',i);
    if (n == string::npos) {
        n = input.length();
    }
    string tmp = input.substr(i,n-i);
    result.push_back(tmp);
  }

}

/// \brief Print some help information to the user
void help()
{
   cout << "   Options: "                                                     << endl;
   cout << "   --channels    channels to be plotted"                          << endl;
   cout << "                 currently implemented ee/mm/em"                  << endl;
   cout << "                 default is mm"                                   << endl;

   cout << "   --regions     regions to be plotted "                           << endl;
   cout << "                 currently implemented --"                         << endl;
   cout << "                 default is CRzjets"                               << endl;

   cout << "   --variables   variables to be plotted "                         << endl;
   cout << "                 currently ptl0, ptl1, etal0, etal1, phil0, phil1" << endl;
   cout << "                 default is ptl0"                                  << endl;

   cout << "   --systematics systematics to be added to the error band"       << endl;
   cout << "                 currently only Tree based systematics"           << endl;

   cout << "   --samples     samples to be included in the plot"              << endl;
   cout << "                 currently implemented Zmumu, Data (default)"     << endl;
   cout << "                 order is preserved for the stack"                << endl;

   cout << "   --binValues   follows the logic: nbin,xmin,xmax"               << endl;
   cout << "                 default is 40,0,400"                             << endl;

   cout << "   --isGeV       converts variable into GeV"                      << endl;
   cout << "                 default is false"                                << endl;          

   cout << "   --plotLog     plots on a log scale"                            << endl;
   cout << "                 default is false"                                << endl;          
}
