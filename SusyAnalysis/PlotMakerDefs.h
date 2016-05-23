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
TS2TSMap RegionCuts = boost::assign::map_list_of ("CR2LOS"         , "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&mll>60.")
                                                 ("CRTTbarLike"    , "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&nCentralBJets>=1&&mll>20." )
                                                 ("CRWWLike"       , "nBaseLeptons==2&&l_pt[0]>25.&&l_pt[1]>25.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0" )
                                                 ("SR2LA-pre"      , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0" )
                                                 ("SR2LA-V1"       , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0&&R1>0.2&&mT2lep>30.&&TMath::Abs(cthllb)<0.8" )
                                                 ("SR2LA-V2"       , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralBJets==0&&R1>0.2&&mT2lep>30.&&TMath::Abs(cthllb)<0.8&&DPB>1.5" )
                                                 ("SR2L-Danny"     , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nStop2lBJets==0&&MDR_jigsaw>95.&&RPT_jigsaw>0.5&&gamInvRp1_jigsaw>0.8&&DPB_vSS_jigsaw>(0.8*TMath::Abs(cthllb)+1.8)" )
                                                 ("SR2L-Serhan-lowDM"  , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&mT2lep>90.&&R2>0.7&&nStop2lBJets==0" )
                                                 ("SR2L-Serhan-highDM" , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&mT2lep>100.&&R2>0.7" )
                                                 ////("SR2L-Serhan"    , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&mT2lep>100.&&R1>0.35&&R2>(R1+0.2)" )
                                                 ////("SR2L-Serhan"    , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&mT2lep>100.&&R1>0.35&&R2>(R1+0.2)&&met>180." )
                                                 ////("SR2L-Serhan"    , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&mT2lep>100.&&R1>0.35&&R2>(R1+0.2)&&ptll<80." )
                                                 ////("SR2L-Serhan"    , "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&mT2lep>100.&&R1>0.35&&R2>(R1+0.2)&&RPT_jigsaw>0.6" )
                                                 ("SR2LpreMT2"    , "l_pt[0]>25.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralLJets==0&&nCentralBJets==0&&nForwardJets==0&&!(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<10.)&&mll>20." )
                                                 ("SR2LEWK-mT2-150", "l_pt[0]>25.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&nCentralLJets==0&&nCentralBJets==0&&nForwardJets==0&&!(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<10.)&&mll>20.&&mT2lep>150." )
;

/// \brief Sample Colors
ColorMap SampleColors = boost::assign::map_list_of ("Data"     , (int)kBlack     )
                                                   ("ttbar"    , (int)kAzure-9   )
                                                   ("singletop", (int)kOrange-3  )
                                                   ("W"        , (int)kGray      )
                                                   ("Z"        , (int)kYellow-9  )
                                                   ("VV"       , (int)kSpring+1  )
                                                   ("392505"   , (int)kGreen+1   )
                                                   ("392510"   , (int)kMagenta+1 )
                                                   ("406009"   , (int)kRed+1     )
                                                   ("406010"   , (int)kGreen+1   )
                                                   ("406011"   , (int)kMagenta   )
;

/// \brief Sample Names
TS2TSMap SampleNames = boost::assign::map_list_of ("Data"     , "Data 2015"        )
                                                  ("ttbar"    , "t#bar{t}"         )
                                                  ("singletop", "Single Top"       )
                                                  ("W"        , "W"                )
                                                  ("Z"        , "Z"                )
                                                  ("VV"       , "VV"               )
                                                  ("392505"   , "C1C1 300, 200 GeV"  )
                                                  ("392510"   , "C1C1 500, 1 GeV"  )
                                                  ("406009"   , "t1t1 250, 160 GeV")
                                                  ("406010"   , "t1t1 300, 150 GeV")
                                                  ("406011"   , "t1t1 300, 180 GeV")
;


/// \brief Define the variables
TS2TSMap VariableNames = boost::assign::map_list_of ("ptL0"            , "l_pt[0]"           )
                                                    ("ptL1"            , "l_pt[1]"           )
                                                    ("etaL0"           , "l_eta[0]"          )
                                                    ("etaL1"           , "l_eta[1]"          )
                                                    ("phiL0"           , "l_phi[0]"          )
                                                    ("phiL1"           , "l_phi[1]"          )
                                                    ("ptvarcone20L0"   , "l_ptvarcone20[0]"  )
                                                    ("ptvarcone20L1"   , "l_ptvarcone20[1]"  )
                                                    ("ptvarcone30L0"   , "l_ptvarcone30[0]"  )
                                                    ("ptvarcone30L1"   , "l_ptvarcone30[1]"  )
                                                    ("mll"             , "mll"               )
                                                    ("ptll"            , "ptll"              )
                                                    ("drll"            , "drll"              )
                                                    ("dphill"          , "dphill"            )
                                                    ("ptJ0"            , "j_pt[0]"           )
                                                    ("ptJ1"            , "j_pt[1]"           )
                                                    ("etaJ0"           , "j_eta[0]"          )
                                                    ("etaJ1"           , "j_eta[1]"          )
                                                    ("phiJ0"           , "j_phi[0]"          )
                                                    ("phiJ1"           , "j_phi[1]"          )
                                                    ("flavJ0"          , "j_flav[0]"         )
                                                    ("flavJ1"          , "j_flav[1]"         )
                                                    ("met"             , "met"               )
                                                    ("meff"            , "meff"              )
                                                    ("mT2lep"          , "mT2lep"            )
                                                    ("R1"              , "R1"                )
                                                    ("R2"              , "R2"                )
                                                    ("MDR_jigsaw"      , "MDR_jigsaw"        )
                                                    ("RPT_jigsaw"      , "RPT_jigsaw"        )
                                                    ("gamInvRp1_jigsaw", "gamInvRp1_jigsaw"  )
                                                    ("DPB_vSS_jigsaw"  , "DPB_vSS_jigsaw"    )
                                                    ("abs_deltaX"      , "TMath::Abs(deltaX)")
                                                    ("abs_cthllb"      , "TMath::Abs(cthllb)")
                                                    ("nBaseJets"       , "nBaseJets"         )
                                                    ("nCentralLJets"   , "nCentralLJets"     )
                                                    ("nCentralBJets"   , "nCentralBJets"     )
                                                    ("nForwardJets"    , "nForwardJets"      )
                                                    ("nStop2lLJets"    , "nStop2lLJets"      )
                                                    ("nStop2lBJets"    , "nStop2lBJets"      )
;

#endif
