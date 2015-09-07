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

// std include(s)
//#include <cstdlib>  // << for atoi
//#include "unistd.h" // << for getopt

// analysis include(s)
#include "Superflow/Superflow.h"     
#include "Superflow/Superlink.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"
#include "SusyAnalysis/MT2_ROOT.h"

using namespace sflow;

///////////////////////////////////////////////////////////////////////
// Usage
///////////////////////////////////////////////////////////////////////
void usage(std::string progName)
{
  printf("=================================================================\n");
  printf("%s [options]\n",progName.c_str());
  printf("=================================================================\n");
  printf("Options:\n");
  printf("-h        Print this help\n");
  printf("-n        Number of events to be processed (default: -1)\n");
  printf("-f        Input file as *.root, list of *.root in a *.txt,\n"); 
  printf("          or a DIR/ containing *.root (default: none)\n");
  printf("=================================================================\n");
}

///////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  ///////////////////////////////////////////////////////////////////////
  // Read user inputs - NOT super safe so be careful :)
  unsigned int n_events  = -1;
  unsigned int n_skip_events  = 0;
  char *input_file = NULL;
  SuperflowRunMode run_mode = SuperflowRunMode::nominal;
  int c;

  opterr = 0;
  while ((c = getopt (argc, argv, "f:n:h")) != -1)
    switch (c)
      {
      case 'f':
        input_file = optarg;
        break;
      case 'n':
        n_events = atoi(optarg);
        break;
      case 'h':
        usage("makeMiniNtuple");
        return 1;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "makeMiniNtuple\t Option -%c requires an argument.\n", optopt);
        else if (optopt == 'n')
          fprintf (stderr, "makeMiniNtuple\t Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "makeMiniNtuple\t Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "makeMiniNtuple\t Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  // Catch problems or cast
  for (int index = optind; index < argc; index++)
    printf ("makeMiniNtuple\t Non-option argument %s\n", argv[index]);
  if (input_file==NULL) {
    printf("makeMiniNtuple\t An input file must be provided with option -f (a list, a DIR or single file)\n");
    return 0;  
  }
  

  // Print information
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t Running SusyAnalysis/makeMiniNtuple\n");
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t   Flags:\n");
  printf("makeMiniNtuple\t     Input file (-f)       : %s\n",input_file);
  printf("makeMiniNtuple\t     Number of events (-n) : %i\n",n_events );
  printf("makeMiniNtuple\t =================================================================\n");
  ///////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  // Setup Superflow 
  TChain* chain = new TChain("susyNt");
  chain->SetDirectory(0);

  bool inputIsFile = Susy::utils::endswith(input_file, ".root");
  bool inputIsList = Susy::utils::endswith(input_file, ".txt");
  bool inputIsDir  = Susy::utils::endswith(input_file, "/");

  if(inputIsFile) {
    ChainHelper::addFile(chain, input_file);
  } else if (inputIsList) {
    // If a list of ROOT files
    ChainHelper::addFileList(chain, input_file);
  } else if (inputIsDir) {
    ChainHelper::addFileDir(chain, input_file);
  } else {
    printf("makeMiniNtuple\t Cannot understand input %s",input_file);
    return 0;  
  }

  Superflow* cutflow = new Superflow(); // initialize the cutflow
  cutflow->setAnaName("SuperflowAna");                // arbitrary
  //cutflow->setAnaType(AnalysisType::Ana_2Lep);        // analysis type, passed to SusyNt ?
  cutflow->setAnaType(AnalysisType::Ana_Stop2L);        // analysis type, passed to SusyNt ?
  cutflow->setLumi(78.287);                           // set the MC normalized to X pb-1
  cutflow->setSampleName(input_file);                 // sample name, check to make sure it's set OK
  cutflow->setRunMode(run_mode);                      // make configurable via run_mode
  cutflow->setCountWeights(true);                     // print the weighted cutflows
  cutflow->setChain(chain);

  printf("makeMiniNtuple\t Total events available : %lli\n",chain->GetEntries());

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods begin here
  ///////////////////////////////////////////////////////////////////////

  *cutflow << CutName("read in") << [](Superlink* /*sl*/) -> bool { return true; };

  //  Cleaning Cuts
  int cutflags = 0;
  
  *cutflow << CutName("Pass GRL") << [&](Superlink* sl) -> bool {
      cutflags = sl->nt->evt()->cutFlags[sl->nt_sys];
      return (cutflags & ECut_GRL);
  };

  *cutflow << CutName("LAr error") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_LarErr);
  };
  
  *cutflow << CutName("Tile error") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_TileErr);
  };
  
  *cutflow << CutName("TTC veto") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_TTC);
  };

  *cutflow << CutName("pass good vertex") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_GoodVtx);
  };
  
  *cutflow << CutName("bad muon veto") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_BadMuon);
  };
  
  *cutflow << CutName("pass cosmic veto") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_Cosmic);
  };

  *cutflow << CutName("jet cleaning") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_BadJet);
  };
  
  //  Analysis Cuts
  *cutflow << CutName("exactly two baseline leptons") << [](Superlink* sl) -> bool {
      return (sl->baseLeptons->size() == 2);
  };

  *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
      return (sl->leptons->size() == 2);
  };

  //  Output Ntuple Setup
  //      > Ntuple variables

  // Event variables
  *cutflow << NewVar("event weight"); {
      *cutflow << HFTname("eventweight");
      *cutflow << [](Superlink* sl, var_double*) -> double { 
          return sl->weights->product();
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("Event run number"); {
    *cutflow << HFTname("runNumber");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->run; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("Event number"); {
    *cutflow << HFTname("eventNumber");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->nt->evt()->eventNumber; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("is Monte Carlo"); {
    *cutflow << HFTname("isMC");
    *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->nt->evt()->isMC ? true : false; };
    *cutflow << SaveVar();
  }

  // Lepton variables
  LeptonVector baseLeptons, signalLeptons;
  *cutflow << [&](Superlink* sl, var_void*) { 
    baseLeptons   = *sl->baseLeptons; 
    signalLeptons = *sl->leptons; 
  };

  *cutflow << NewVar("number of baseline leptons"); {
    *cutflow << HFTname("nBaseLeptons");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseLeptons.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of signal leptons"); {
    *cutflow << HFTname("nSignalLeptons");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return signalLeptons.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton flavor (0: e, 1: m)"); {
    *cutflow << HFTname("l_flav");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->isEle() ? 0 : 1);
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton pt"); {
    *cutflow << HFTname("l_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->Pt());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton eta"); {
    *cutflow << HFTname("l_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->Eta());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton phi"); {
    *cutflow << HFTname("l_phi");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->Phi());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton d0"); {
    *cutflow << HFTname("l_d0");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->d0);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton errD0"); {
    *cutflow << HFTname("l_errD0");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->errD0);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton d0sig"); {
    *cutflow << HFTname("l_d0sig");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->d0Sig());
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton z0"); {
    *cutflow << HFTname("l_z0");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->z0);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton errZ0"); {
    *cutflow << HFTname("l_errZ0");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->errZ0);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton z0sinTheta"); {
    *cutflow << HFTname("l_z0sinTheta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->z0SinTheta());
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton charge"); {
    *cutflow << HFTname("l_q");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->q);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton ptvarcone20"); {
    *cutflow << HFTname("l_ptvarcone20");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->ptvarcone20);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton ptvarcone30"); {
    *cutflow << HFTname("l_ptvarcone30");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
          out.push_back(lepton->ptvarcone30);
      }
      return out;
      };
    *cutflow << SaveVar();
  }

  // Leptonic event variables
  TLorentzVector lepton0 ; 
  TLorentzVector lepton1 ; 
  TLorentzVector dileptonP4 ;
  *cutflow << [&](Superlink* /*sl*/, var_void*) {
    lepton0 = *signalLeptons.at(0); 
    lepton1 = *signalLeptons.at(1); 
    dileptonP4 = lepton0 + lepton1;
  };

  *cutflow << NewVar("mass of di-lepton system, M_ll"); {
    *cutflow << HFTname("mll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return dileptonP4.M(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("Pt of di-lepton system, Pt_ll"); {
    *cutflow << HFTname("ptll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return dileptonP4.Pt(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("delta Phi of di-lepton system"); {
    *cutflow << HFTname("dphill");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.DeltaPhi(lepton1); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("delta R of di-lepton system"); {
    *cutflow << HFTname("drll");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return lepton0.DeltaR(lepton1); };
    *cutflow << SaveVar();
  }

  // Jet variables
  JetVector baseJets; //, centralLightJets, centralBJets, forwardJets;
  *cutflow << [&](Superlink* sl, var_void*) { 
    baseJets = *sl->baseJets; 
    //for(auto& jet : baseJets) {
    //  if(sl->tools->m_jetSelector.isCentralLightJet(jet))  { centralLightJets.push_back(jet); } 
    //  else if(sl->tools->m_jetSelector.isCentralBJet(jet)) { centralBJets.push_back(jet);     } 
    //  else if(sl->tools->m_jetSelector.isForwardJet(jet))  { forwardJets.push_back(jet);      } 
    //}
  };

  *cutflow << NewVar("number of baseline jets"); {
    *cutflow << HFTname("nBaseJets");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseJets.size(); };
    *cutflow << SaveVar();
  }

  //*cutflow << NewVar("number of central light jets"); {
  //  *cutflow << HFTname("nCentralLJets");
  //  *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCLJets(*sl->baseJets)/*(*baseJets)*/; };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("number of central b jets"); {
  //  *cutflow << HFTname("nCentralBJets");
  //  *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCBJets(*sl->baseJets)/*(*baseJets)*/; };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("number of forward jets"); {
  //  *cutflow << HFTname("nForwardJets");
  //  *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->baseJets)/*(*baseJets)*/; };
  //  *cutflow << SaveVar();
  //}

  *cutflow << NewVar("jet pt"); {
    *cutflow << HFTname("j_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& jet : baseJets) {
        out.push_back(jet->Pt());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("jet eta"); {
    *cutflow << HFTname("j_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& jet : baseJets) {
        out.push_back(jet->Eta());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("jet phi"); {
    *cutflow << HFTname("j_phi");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& jet : baseJets) {
        out.push_back(jet->Phi());
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("jet flavor (0: NA, 1: CL, 2: CB, 3: F)"); {
    *cutflow << HFTname("j_flav");
    *cutflow << [&](Superlink* sl, var_float_array*) -> vector<double> {
      vector<double> out; int flav = 0;
      for(auto& jet : baseJets) {
        if(sl->tools->m_jetSelector.isCentralLightJet(jet))  { flav = 1; } 
        else if(sl->tools->m_jetSelector.isCentralBJet(jet)) { flav = 2; } 
        else if(sl->tools->m_jetSelector.isForwardJet(jet))  { flav = 3; } 
        out.push_back(flav);
        flav=0;
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  // MET variables
  TLorentzVector met;
  *cutflow << NewVar("transverse missing energy (Et)"); {
    *cutflow << HFTname("met");
    *cutflow << [&](Superlink* sl, var_float*) -> double { 
      met.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
                     sl->met->Et*sin(sl->met->phi),
                     0.,
                     sl->met->Et);
      return met.Pt();
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("transverse missing energy (Phi)"); {
    *cutflow << HFTname("metPhi");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { return met.Phi(); };
    *cutflow << SaveVar();
  }

  // Other event variables

  // meff 
  double meff = 0.;
  *cutflow << NewVar("scalar sum pt of all leptons, jets and met"); {
    *cutflow << HFTname("meff");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      // Leptons
      for(auto& lepton : signalLeptons) {
        meff += lepton->Pt();
      }
      // Jets
      for(auto& jet : baseJets) {
        meff += jet->Pt();
      }
      // MET
      meff += met.Pt();
      return meff;
    };
    *cutflow << SaveVar();
  }

  // R1   
  double R1 = 0.;
  *cutflow << NewVar("met/meff"); {
    *cutflow << HFTname("R1");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      if(fabs(meff-1.e-5)>0.) R1 = met.Pt()/meff; 
      return R1; 
    };
    *cutflow << SaveVar();
  }

  // R2   
  double R2 = 0.;
  *cutflow << NewVar("met/(met + sum(lepton,pt))"); {
    *cutflow << HFTname("R2");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      double sumLepPt = 0.;
      // Leptons
      for(auto& lepton : signalLeptons) {
        sumLepPt += lepton->Pt();
      }
      if ( fabs((met.Pt()+sumLepPt)-1.e-5)>0. ) R2 = met.Pt()/(met.Pt()+sumLepPt);
      return R2; 
    };
    *cutflow << SaveVar();
  }

  // deltaX  
  double deltaX = 0.;
  *cutflow << NewVar("pz0 + pz1 / sqrt(s)"); {
    *cutflow << HFTname("deltaX");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      deltaX = (lepton0.Pz() + lepton1.Pz())/sqrt(13000.);
      return deltaX;
    };
    *cutflow << SaveVar();
  }

  // mT2  
  double mT2 = 0.;
  *cutflow << NewVar("stransvers mass of two leptons"); {
    *cutflow << HFTname("mT2lep");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      ComputeMT2 mycalc = ComputeMT2(lepton0,lepton1,met,0.,0.); // masses 0. 0.
      mT2 = mycalc.Compute();
      return mT2; 
    };
    *cutflow << SaveVar();
  }
 
  // cosTheta_b  
  double cthllb = 0.;
  *cutflow << NewVar("cosTheta_b"); {
    *cutflow << HFTname("cthllb");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      if(signalLeptons.at(0)->q*signalLeptons.at(1)->q > 0.) return -999.; // SS event
      TLorentzVector lepPos, lepNeg; 
      if(signalLeptons.at(0)->q > 0.) { 
        lepPos = lepton0; lepNeg = lepton1; 
      }      
      else { 
        lepPos = lepton1; lepNeg = lepton0; 
      }      
      TVector3 boost = dileptonP4.BoostVector();
      lepPos.Boost(-boost);
      lepNeg.Boost(-boost);
      cthllb = tanh((lepPos.Eta()-lepNeg.Eta())/2);
      return cthllb; 
    };
    *cutflow << SaveVar();
  }

  // Clear 
  *cutflow << [&](Superlink* /*sl*/, var_void*) { 
    signalLeptons.clear(); 
    baseJets.clear();
    meff=0.,R1=0.,R2=0.,deltaX=0.,mT2=0.,cthllb=0.;
    //centralLightJets.clear();
    //centralBJets.clear();
    //forwardJets.clear();
  };

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods end here
  ///////////////////////////////////////////////////////////////////////

  // Initialize the cutflow and start the event loop.
  chain->Process(cutflow, input_file, n_events, n_skip_events);

  delete cutflow;
  delete chain;

  // Print information
  printf("makeMiniNtuple\t =================================================================\n");
  printf("makeMiniNtuple\t All done!\n");
  printf("makeMiniNtuple\t =================================================================\n");
  return 0; 
}
