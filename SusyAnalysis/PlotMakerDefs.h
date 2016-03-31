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
TS2TSMap RegionCuts = boost::assign::map_list_of ("CR2LOS"     , "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&mll>60.")
                                                 ("CRTTbarLike", "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&nCentralBJets>=1&&mll>20." )
                                                 ("CRWWLike"   , "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0" )
                                                 ("SR2LA-pre"  , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0" )
                                                 ("SR2LA-V1"   , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0&&R1>0.2&&mT2lep>30.&&TMath::Abs(cthllb)<0.8" )
                                                 ("SR2LA-V2"   , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0&&R1>0.2&&mT2lep>30.&&TMath::Abs(cthllb)<0.8&&DPB>1.5" )
                                                 ("SR2L-Danny" , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0&&mT2lep>80.&&R2>0.65&&DPB>(1.0*TMath::Abs(cthllb)+2)" )
                                                 ("SR2LEWK-pre", "l_pt[0]>25.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralLJets==0&&nCentralBJets==0&&nForwardJets==0&&!(mll<60.||(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<10.))" )
                                                 ("SR2LEWK-mT2-150", "l_pt[0]>25.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralLJets==0&&nCentralBJets==0&&nForwardJets==0&&!(mll<60.||(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<10.))&&mT2lep>150." )
;

/// \brief Sample Colors
ColorMap SampleColors = boost::assign::map_list_of ("Data"     , (int)kBlack     )
                                                   ("ttbar"    , (int)kAzure-9   )
                                                   ("singletop", (int)kRed+1     )
                                                   ("W"        , (int)kGray      )
                                                   ("Z"        , (int)kYellow-9  )
                                                   ("VV"       , (int)kGreen-3   )
                                                   //("WZ"       , (int)kOrange-2  )
                                                   //("ZZ"       , (int)kMagenta+1 )
                                                   ("406009"   , (int)kBlue+2    )
                                                   ("392510"   , (int)kBlue+2    )
;

/// \brief Sample Names
TS2TSMap SampleNames = boost::assign::map_list_of ("Data"     , "Data 2015"  )
                                                  ("ttbar"    , "t#bar{t}"   )
                                                  ("singletop", "Single Top" )
                                                  ("W"        , "W"          )
                                                  ("Z"        , "Z"          )
                                                  ("VV"       , "VV"         )
                                                  //("WZ"       , "WZ"         )
                                                  //("ZZ"       , "ZZ"         )
                                                  ("406009"   , "Stop 250, 160 GeV" )
                                                  //("406009"   , "(m_{#tilde{t}_1}, m_{#tilde{#chi}_{1}^{0}}) = 250, 160 GeV" )
                                                  ("392510"   , "C1C1 500, 1 GeV" )
;


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
                                                    ("ptJ0"          , "j_pt[0]"          )
                                                    ("ptJ1"          , "j_pt[1]"          )
                                                    ("etaJ0"         , "j_eta[0]"         )
                                                    ("etaJ1"         , "j_eta[1]"         )
                                                    ("phiJ0"         , "j_phi[0]"         )
                                                    ("phiJ1"         , "j_phi[1]"         )
                                                    ("flavJ0"        , "j_flav[0]"        )
                                                    ("flavJ1"        , "j_flav[1]"        )
                                                    ("met"           , "met"              )
                                                    ("meff"          , "meff"             )
                                                    ("mT2lep"        , "mT2lep"           )
                                                    ("R1"            , "R1"               )
                                                    ("R2"            , "R2"               )
                                                    ("abs_deltaX"    , "TMath::Abs(deltaX)" )
                                                    ("abs_cthllb"    , "TMath::Abs(cthllb)" )
                                                    ("nBaseJets"     , "nBaseJets"        )
                                                    ("nCentralLJets" , "nCentralLJets"    )
                                                    ("nCentralBJets" , "nCentralBJets"    )
                                                    ("nForwardJets"  , "nForwardJets"     )

;

#endif
