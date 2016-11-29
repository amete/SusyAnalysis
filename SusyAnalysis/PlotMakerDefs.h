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

TString trigger   = "(((pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14||pass_HLT_mu18_mu8noL1)&&treatAsYear==2015)||((pass_HLT_2e17_lhvloose_nod0||pass_HLT_e17_lhloose_nod0_mu14||pass_HLT_mu22_mu8noL1)&&treatAsYear==2016))&&((!isMC&&runNumber<303560)||isMC)";
TString ptCuts    = "l_pt[0]>25.&&l_pt[1]>20.&&mll>40.";
TString isOS      = "(l_q[0]*l_q[1])<0";
TString zVeto     = "!(l_flav[0]==l_flav[1]&&TMath::Abs(mll-91.2)<10.)";
TString zSelect   = "TMath::Abs(mll-91.2)<10.";
TString cljVeto   = "((l_flav[0]==l_flav[1]&&nCentralLJets==0)||(l_flav[0]!=l_flav[1]&&nCentralLJets30==0))";
TString cljVeto20 = "nCentralLJets==0";
TString cljVeto30 = "nCentralLJets30==0";
TString cljVeto40 = "nCentralLJets40==0";
TString cljVeto50 = "nCentralLJets50==0";
TString cljVeto60 = "nCentralLJets60==0";
TString cbjVeto   = "nCentralBJets==0";
TString cbjSelect = "nCentralBJets>0";
//TString fjVeto    = "nForwardJets==0";
TString fjVeto    = "nForwardJets30==0"; // in n0228 the nForwardJets has 20 GeV threshold!!!

/// \brief Region Cut Definitions
TS2TSMap RegionCuts = boost::assign::map_list_of ("CR2L-inc"       , trigger                           ) // Skimming requires ==2 signal leptons w/ pT > 20 GeV and OS + mll > 40 GeV
                                                 ("CR2L-bselect"   , trigger+"&&"+cbjSelect            ) // Same skimming as above
                                                 ("CR2L-bveto"     , trigger+"&&"+cbjVeto              ) // Same skimming as above
                                                 ("CR2L-zbveto"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto +"&&nStop2lBJets==0" ) // Same skimming as above
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
                                                 ("SR2L-preMT2"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto )
                                                 ("SR2L-preMT250"  , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&mT2lep>50." )
                                                 ("SR2L-MT290"     , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&mT2lep>90." )
                                                 ("SR2L-MT2120"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&mT2lep>120." )
                                                 ("SR2L-MT2150"    , trigger+"&&"+ptCuts+"&&"+isOS+"&&"+zVeto  +"&&"+cljVeto  +"&&"+cbjVeto  +"&&"+fjVeto+"&&mT2lep>150." )
                                                 ("SR2L-preMT2J-new", trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cbjVeto )
                                                 ("SR2L-preMT2-new" , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto )
                                                 ("SR2L-NewSF"      , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto  +"&&l_flav[0]==l_flav[1]&&mll>100&&mT2lep>100"    )
                                                 ("CR2L-NewVVSF"    , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto  +"&&l_flav[0]==l_flav[1]&&mll<100&&mT2lep>100."   )
                                                 ("VR2L-NewVVSF"    , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto  +"&&l_flav[0]==l_flav[1]&&mll>100&&mT2lep>75.&&mT2lep<100." )
                                                 ("SR2L-NewDF"      , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto  +"&&l_flav[0]!=l_flav[1]&&mll>40&&mT2lep>100"     )
                                                 ("CR2L-NewVVDF"    , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto  +"&&l_flav[0]!=l_flav[1]&&mll>40&&mT2lep>50.&&mT2lep<75." )
                                                 ("VR2L-NewVVDF"    , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjVeto  +"&&l_flav[0]!=l_flav[1]&&mll>40&&mT2lep>75.&&mT2lep<100.")
                                                 ("CR2L-NewTop"     , trigger+"&&"+ptCuts+"&&"+isOS             +"&&"+cljVeto60+"&&"+cbjSelect+"&&l_flav[0]!=l_flav[1]&&mll>40&&mT2lep>70.&&mT2lep<120.")
;

/// \brief Sample Colors
ColorMap SampleColors = boost::assign::map_list_of ("Data"                 , (int)kBlack     )
                                                   ("ttbar"                , (int)kAzure-9   )
                                                   ("singletop"            , (int)kOrange-3  )
                                                   ("ttv"                  , (int)kRed+1     )
                                                   ("higgs"                , (int)kBlue+1    )
                                                   ("W"                    , (int)kGray      )
                                                   ("MM"                   , (int)kGray      )
                                                   ("Z"                    , (int)kYellow-9  )
                                                   ("VV"                   , (int)kSpring+1  )
                                                   ("VVV"                  , (int)kSpring+3  )
                                                   ("c1c1_slep_200.0_100.0", (int)kGreen+1   )
                                                   ("c1c1_slep_300.0_100.0", (int)kGreen+1   )
                                                   ("c1c1_slep_300.0_200.0", (int)kGreen+1   )
                                                   ("c1c1_slep_400.0_200.0", (int)kGreen+1   )
                                                   ("c1c1_slep_500.0_1.0"  , (int)kMagenta+1 )
                                                   ("c1c1_slep_700.0_1.0"  , (int)kMagenta+1 )
                                                   ("SlepSlep_350.5_1.0"   , (int)kGreen+1   )
                                                   ("SlepSlep_450.5_1.0"   , (int)kGreen+1   )
                                                   ("406009"               , (int)kRed+1     )
                                                   ("406010"               , (int)kGreen+1   )
                                                   ("406011"               , (int)kMagenta   )
;

/// \brief Sample Names
TS2TSMap SampleNames = boost::assign::map_list_of ("Data"                  , "Data"             )
                                                  ("ttbar"                 , "t#bar{t}"         )
                                                  ("singletop"             , "Single Top"       )
                                                  ("ttv"                   , "t#bar{t}+V"       )
                                                  ("W"                     , "W"                )
                                                  ("MM"                    , "Fakes"            )
                                                  ("higgs"                 , "Higgs"            )
                                                  ("Z"                     , "Z"                )
                                                  ("VV"                    , "VV"               )
                                                  ("VVV"                   , "VVV"              )
                                                  ("c1c1_slep_200.0_100.0" , "C1C1 200, 100 GeV")
                                                  ("c1c1_slep_300.0_100.0" , "C1C1 300, 100 GeV")
                                                  ("c1c1_slep_300.0_200.0" , "C1C1 300, 200 GeV")
                                                  ("c1c1_slep_400.0_200.0" , "C1C1 400, 200 GeV")
                                                  ("c1c1_slep_500.0_1.0"   , "C1C1 500, 1 GeV"  )
                                                  ("c1c1_slep_700.0_1.0"   , "C1C1 700, 1 GeV"  )
                                                  ("SlepSlep_350.5_1.0"    , "SlSl 350.5, 1 GeV"  )
                                                  ("SlepSlep_450.5_1.0"    , "SlSl 450.5, 1 GeV"  )
                                                  ("406009"                , "t1t1 250, 160 GeV")
                                                  ("406010"                , "t1t1 300, 150 GeV")
                                                  ("406011"                , "t1t1 300, 180 GeV")
;


/// \brief Define the variables
TS2TSMap VariableNames = boost::assign::map_list_of ("ptL0"                , "l_pt[0]"             )
                                                    ("ptL1"                , "l_pt[1]"             )
                                                    ("etaL0"               , "l_eta[0]"            )
                                                    ("etaL1"               , "l_eta[1]"            )
                                                    ("phiL0"               , "l_phi[0]"            )
                                                    ("phiL1"               , "l_phi[1]"            )
                                                    ("ptvarcone20L0"       , "l_ptvarcone20[0]"    )
                                                    ("ptvarcone20L1"       , "l_ptvarcone20[1]"    )
                                                    ("ptvarcone30L0"       , "l_ptvarcone30[0]"    )
                                                    ("ptvarcone30L1"       , "l_ptvarcone30[1]"    )
                                                    ("mll"                 , "mll"                 )
                                                    ("ptll"                , "ptll"                )
                                                    ("drll"                , "drll"                )
                                                    ("dphill"              , "dphill"              )
                                                    ("ptJ0"                , "j_pt[0]"             )
                                                    ("ptJ1"                , "j_pt[1]"             )
                                                    ("etaJ0"               , "j_eta[0]"            )
                                                    ("etaJ1"               , "j_eta[1]"            )
                                                    ("phiJ0"               , "j_phi[0]"            )
                                                    ("phiJ1"               , "j_phi[1]"            )
                                                    ("flavJ0"              , "j_flav[0]"           )
                                                    ("flavJ1"              , "j_flav[1]"           )
                                                    ("met"                 , "met"                 )
                                                    ("metRel"              , "metRel"              )
                                                    ("meff"                , "meff"                )
                                                    ("mT2lep"              , "mT2lep"              )
                                                    ("R1"                  , "R1"                  )
                                                    ("R2"                  , "R2"                  )
                                                    ("S1"                  , "S1"                  )
                                                    ("S2"                  , "S2"                  )
                                                    ("MDR_jigsaw"          , "MDR_jigsaw"          )
                                                    ("RPT_jigsaw"          , "RPT_jigsaw"          )
                                                    ("gamInvRp1_jigsaw"    , "gamInvRp1_jigsaw"    )
                                                    ("DPB_vSS_jigsaw"      , "DPB_vSS_jigsaw"      )
                                                    ("abs_deltaX"          , "TMath::Abs(deltaX)"  )
                                                    ("abs_cthllb"          , "TMath::Abs(cthllb)"  )
                                                    ("nBaseJets"           , "nBaseJets"           )
                                                    ("nCentralLJets"       , "nCentralLJets"       )
                                                    ("nCentralLJets30"     , "nCentralLJets30"     )
                                                    ("nCentralLJets40"     , "nCentralLJets40"     )
                                                    ("nCentralLJets50"     , "nCentralLJets50"     )
                                                    ("nCentralLJets60"     , "nCentralLJets60"     )
                                                    ("nCentralBJets"       , "nCentralBJets"       )
                                                    ("nForwardJets"        , "nForwardJets"        )
                                                    ("nForwardJets30"      , "nForwardJets30"      )
                                                    ("nForwardJets50"      , "nForwardJets50"      )
                                                    ("nStop2lLJets"        , "nStop2lLJets"        )
                                                    ("nStop2lBJets"        , "nStop2lBJets"        )
                                                    ("softestCentralLJetPt", "softestCentralLJetPt")
                                                    ("softestForwardJetPt" , "softestForwardJetPt" )
                                                    ("hardestCentralLJetPt", "hardestCentralLJetPt")
                                                    ("hardestForwardJetPt" , "hardestForwardJetPt" )
                                                    ("dphi_ll_met"         , "dphi_ll_met"         )
                                                    ("dphi_ljet_met"       , "dphi_ljet_met"       )
;

#endif
