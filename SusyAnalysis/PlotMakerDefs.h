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
;

/// \brief Sample Colors
ColorMap SampleColors = boost::assign::map_list_of ("Data"     , (int)kBlack     )
                                                   ("ttbar"    , (int)kAzure-9   )
                                                   ("singleTop", (int)kRed+1     )
                                                   ("W"        , (int)kGray      )
                                                   ("Z"        , (int)kYellow-9  )
                                                   ("WW"       , (int)kGreen-3   )
                                                   ("WZ"       , (int)kOrange-2  )
                                                   ("ZZ"       , (int)kMagenta+1 )
;

/// \brief Sample Names
TS2TSMap SampleNames = boost::assign::map_list_of ("Data"     , "Data 2015"  )
                                                  ("ttbar"    , "t#bar{t}"   )
                                                  ("singleTop", "Single Top" )
                                                  ("W"        , "W"          )
                                                  ("Z"        , "Z"          )
                                                  ("WW"       , "WW"         )
                                                  ("WZ"       , "WZ"         )
                                                  ("ZZ"       , "ZZ"         )
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
                                                    ("nBaseJets"     , "nBaseJets"        )
                                                    ("nCentralLJets" , "nCentralLJets"    )
                                                    ("nCentralBJets" , "nCentralBJets"    )
                                                    ("nForwardJets"  , "nForwardJets"     )

;

#endif