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

#ifndef REGIONCUTS_H
#define REGIONCUTS_H

#include <map>
#include "boost/assign.hpp"
#include "TString.h"
#include "Rtypes.h"
#include "TMath.h"

using namespace std;

typedef map<TString, TString> TS2TSMap;
typedef map<TString, int>     ColorMap;

/// \brief Region Cut Definitions
TS2TSMap RegionCuts = boost::assign::map_list_of ("SRpremT2"  , "isOS&&(L2nCentralLightJets+L2nCentralBJets+L2nForwardJets)==0&&TMath::Abs(L2Mll/1000.-91.2)>10.&&L2METrel>40000.") 
                                                 ("SRmT2a"    , "isOS&&(L2nCentralLightJets+L2nCentralBJets+L2nForwardJets)==0&&TMath::Abs(L2Mll/1000.-91.2)>10.&&L2METrel>40000.&&MT2>90000.")
                                                 ("SRmT2b"    , "isOS&&(L2nCentralLightJets+L2nCentralBJets+L2nForwardJets)==0&&TMath::Abs(L2Mll/1000.-91.2)>10.&&L2METrel>40000.&&MT2>110000.")
                                                 ("TopCR"     , "isOS&&(L2nCentralLightJets+L2nCentralBJets)>1&&L2nCentralBJets>0&&TMath::Abs(L2Mll/1000.-91.2)>10.&&L2METrel>40000.")
                                                 ("CRZXSRmT2a", "isOS&&(L2nCentralLightJets+L2nCentralBJets+L2nForwardJets)==0&&TMath::Abs(L2Mll/1000.-91.2)<10.&&L2METrel>40000.&&MT2>90000.")
                                                 ("CRZXSRmT2b", "isOS&&(L2nCentralLightJets+L2nCentralBJets+L2nForwardJets)==0&&TMath::Abs(L2Mll/1000.-91.2)<10.&&L2METrel>40000.&&MT2>110000.")
                                                 ("CRWW"      , "isOS&&(L2nCentralLightJets+L2nCentralBJets+L2nForwardJets)==0&&TMath::Abs(L2Mll/1000.-91.2)>10.&&L2METrel>40000.&&MT2>50000.&&MT2<90000.")
                                                 ("CRzjets"   , "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&mll>60.");

/// \brief Sample Colors
ColorMap SampleColors = boost::assign::map_list_of ("Data"  , (int)kBlack    )
                                                   ("ttst"  , (int)kRed+1    )
                                                   ("WW"    , (int)kAzure+4  )
                                                   ("Zjets" , (int)kOrange-2 )
                                                   ("ZV"    , (int)kSpring+1 )
                                                   ("Higgs" , (int)kYellow-9 )
                                                   ("Fakes" , (int)kGray     )
                                                   ("Zmumu" , (int)kAzure-9  );

/// \brief Sample Names
TS2TSMap SampleNames = boost::assign::map_list_of ("Data"  , "Data 2015"          )
                                                  ("ttst"  , "t#bar{t} + Wt"      )
                                                  ("WW"    , "WW"                 )
                                                  ("Zjets" , "Z+jets"             )
                                                  ("ZV"    , "ZV"                 )
                                                  ("Higgs" , "Higgs"              )
                                                  ("Fakes" , "Fake leptons"       )
                                                  ("Zmumu" , "Z#rightarrow#mu#mu" );


/// \brief Define the variables
TS2TSMap VariableNames = boost::assign::map_list_of ("ptL0"          , "l_pt[0]"          )
                                                    ("ptL1"          , "l_pt[1]"          )
                                                    ("etaL0"         , "l_eta[0]"         )
                                                    ("etaL1"         , "l_eta[1]"         )
                                                    ("phiL0"         , "l_phi[0]"         )
                                                    ("phiL1"         , "l_phi[1]"         )
                                                    ("ptvarcone20L0" , "l_ptvarcone20[0]" )
                                                    ("ptvarcone20L1" , "l_ptvarcone20[1]" )
                                                    ("ptvarcone30L0" , "l_ptvarcone30[0]" )
                                                    ("ptvarcone30L1" , "l_ptvarcone30[1]" )
                                                    ("mll"           , "mll"              )
                                                    ("ptll"          , "ptll"             )
                                                    ("drll"          , "drll"             )
                                                    ("dphill"        , "dphill"           )
;

#endif
