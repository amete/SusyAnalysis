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

TString trigger   = "(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)";
TString ptCuts    = "l_pt[0]>25.&&l_pt[1]>20.&&mll>20.";
TString isOS      = "(l_q[0]*l_q[1])<0";
TString zVeto     = "!(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<10.)";
TString zSelect   = "TMath::Abs(mll-90.2)<10.";
TString cljVeto   = "((l_flav[0]==l_flav[1]&&nCentralLJets==0)||(l_flav[0]!=l_flav[1]&&nCentralLJets30==0))";
TString cljVeto20 = "nCentralLJets==0";
TString cljVeto30 = "nCentralLJets30==0";
TString cljVeto40 = "nCentralLJets40==0";
TString cljVeto50 = "nCentralLJets50==0";
TString cbjVeto   = "nCentralBJets==0";
TString cbjSelect = "nCentralBJets>0";
TString fjVeto    = "nForwardJets==0";

//(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)&&
//(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)&&
//(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)&&
//(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)&&
//(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)&&
//(pass_HLT_mu18_mu8noL1||pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14)&&

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
                                                 ("CR2L-VVDF-LJ"   , trigger+"&&"+ptCuts+"&&"+isOS                            +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]!=l_flav[1]&&mT2lep>50.&&mT2lep<75." )
                                                 ("VR2L-VVDF-LJ"   , trigger+"&&"+ptCuts+"&&"+isOS                            +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]!=l_flav[1]&&mT2lep>75.&&mT2lep<90." )
                                                 ("CR2L-VVSF-LJ"   , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zSelect               +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]==l_flav[1]&&mT2lep>90." )
                                                 ("CR2L-Top-LJ"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto                 +"&&"+cbjSelect+"&&"+fjVeto+"&&mT2lep>70.&&mT2lep<120." )
                                                 ("SR2L-preMT2-J"  , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  )
                                                 ("SR2L-preMT2-LJ" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +               "&&"+cbjVeto+"&&"+fjVeto )
                                                 ("SR2L-preMT2-20" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto20+"&&"+cbjVeto+"&&"+fjVeto )
                                                 ("SR2L-preMT2-30" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto30+"&&"+cbjVeto+"&&"+fjVeto )
                                                 ("SR2L-preMT2-40" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto40+"&&"+cbjVeto+"&&"+fjVeto )
                                                 ("SR2L-preMT2-50" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto50+"&&"+cbjVeto+"&&"+fjVeto )
                                                 ("SR2L-MT2120-LJ" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +               "&&"+cbjVeto+"&&"+fjVeto+"&&mT2lep>120." )
                                                 ("SR2L-MT2150-LJ" , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +               "&&"+cbjVeto+"&&"+fjVeto+"&&mT2lep>150." )
                                                 // Official EWK2L selection
                                                 ("CR2L-VVDF"      , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]!=l_flav[1]&&mT2lep>50.&&mT2lep<75." )
                                                 ("CR2L-VVSF"      , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zSelect+"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]==l_flav[1]&&mT2lep>90." )
                                                 ("VR2L-VVDF"      , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]!=l_flav[1]&&mT2lep>75.&&mT2lep<90." )
                                                 ("VR2L-VVSF"      , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&l_flav[0]==l_flav[1]&&mT2lep>75.&&mT2lep<90." )
                                                 ("CR2L-Top"       , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto  +"&&"+cbjSelect+"&&"+fjVeto+"&&l_flav[0]!=l_flav[1]&&mT2lep>70.&&mT2lep<120." )
                                                 ("SR2L-preMT2"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto+"&&"+fjVeto )
                                                 ("SR2L-MT290"     , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto+"&&"+fjVeto+"&&mT2lep>90." )
                                                 ("SR2L-MT2120"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto+"&&"+fjVeto+"&&mT2lep>120." )
                                                 ("SR2L-MT2150"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto+"&&"+fjVeto+"&&mT2lep>150." )
;

/// \brief Sample Colors
ColorMap SampleColors = boost::assign::map_list_of ("Data"     , (int)kBlack     )
                                                   ("ttbar"    , (int)kAzure-9   )
                                                   ("singletop", (int)kOrange-3  )
                                                   ("ttv"      , (int)kRed+1     )
                                                   ("W"        , (int)kGray      )
                                                   ("MM"       , (int)kGray      )
                                                   ("Z"        , (int)kYellow-9  )
                                                   ("VV"       , (int)kSpring+1  )
                                                   ("VVV"      , (int)kSpring+3  )
                                                   ("392501"   , (int)kGreen+1   )
                                                   ("392505"   , (int)kGreen+1   )
                                                   ("392508"   , (int)kGreen+1   )
                                                   ("392510"   , (int)kMagenta+1 )
                                                   ("406009"   , (int)kRed+1     )
                                                   ("406010"   , (int)kGreen+1   )
                                                   ("406011"   , (int)kMagenta   )
;

/// \brief Sample Names
TS2TSMap SampleNames = boost::assign::map_list_of ("Data"     , "Data 2015"        )
                                                  ("ttbar"    , "t#bar{t}"         )
                                                  ("singletop", "Single Top"       )
                                                  ("ttv"      , "t#bar{t}+V"       )
                                                  ("W"        , "W"                )
                                                  ("MM"       , "Fakes"            )
                                                  ("Z"        , "Z"                )
                                                  ("VV"       , "VV"               )
                                                  ("VVV"      , "VVV"              )
                                                  ("392501"   , "C1C1 200, 100 GeV")
                                                  ("392505"   , "C1C1 300, 200 GeV")
                                                  ("392508"   , "C1C1 400, 200 GeV")
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
                                                    ("nCentralLJets30" , "nCentralLJets30"   )
                                                    ("nCentralLJets40" , "nCentralLJets40"   )
                                                    ("nCentralLJets50" , "nCentralLJets50"   )
                                                    ("nCentralBJets"   , "nCentralBJets"     )
                                                    ("nForwardJets"    , "nForwardJets"      )
                                                    ("nStop2lLJets"    , "nStop2lLJets"      )
                                                    ("nStop2lBJets"    , "nStop2lBJets"      )
                                                    ("softestCentralLJetPt", "softestCentralLJetPt")
                                                    ("softestForwardJetPt" , "softestForwardJetPt" )
                                                    ("hardestCentralLJetPt", "hardestCentralLJetPt")
                                                    ("hardestForwardJetPt" , "hardestForwardJetPt" )
;

#endif
