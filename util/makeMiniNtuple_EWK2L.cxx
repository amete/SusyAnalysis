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
#include "SusyNtuple/KinematicTools.h"

// RestFrames stuff
#include "RestFrames/RestFrames.hh"

using namespace sflow;
using namespace RestFrames;

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
// Print event information 
///////////////////////////////////////////////////////////////////////
void printEventInformation(Superlink* sl) {
  // event
  printf("==============================================\n");
  sl->nt->evt()->print();
  // pre-obj
  printf("==============================================\n");
  printf("Pre-objects: \n");
  for(auto &obj : *(sl->preElectrons))  { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->preMuons))      { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->preJets))       { obj->print(); }
  // base-obj
  printf("==============================================\n");
  printf("Base-objects: \n");
  for(auto &obj : *(sl->baseElectrons)) { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->baseMuons))     { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->baseJets))      { obj->print(); }
  // signal-obj
  printf("==============================================\n");
  printf("Signal-objects: \n");
  for(auto &obj : *(sl->electrons))     { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->muons))         { obj->print(); }
  printf("\n");
  for(auto &obj : *(sl->jets))          { obj->print(); }
  printf("==============================================\n");
  return;
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
        usage("makeMiniNtuple_EWK2L");
        return 1;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "makeMiniNtuple_EWK2L\t Option -%c requires an argument.\n", optopt);
        else if (optopt == 'n')
          fprintf (stderr, "makeMiniNtuple_EWK2L\t Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "makeMiniNtuple_EWK2L\t Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "makeMiniNtuple_EWK2L\t Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  // Catch problems or cast
  for (int index = optind; index < argc; index++)
    printf ("makeMiniNtuple_EWK2L\t Non-option argument %s\n", argv[index]);
  if (input_file==NULL) {
    printf("makeMiniNtuple_EWK2L\t An input file must be provided with option -f (a list, a DIR or single file)\n");
    return 0;  
  }
  

  // Print information
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
  printf("makeMiniNtuple_EWK2L\t Running SusyAnalysis/makeMiniNtuple_EWK2L\n");
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
  printf("makeMiniNtuple_EWK2L\t   Flags:\n");
  printf("makeMiniNtuple_EWK2L\t     Input file (-f)       : %s\n",input_file);
  printf("makeMiniNtuple_EWK2L\t     Number of events (-n) : %i\n",n_events );
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
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
    printf("makeMiniNtuple_EWK2L\t Cannot understand input %s",input_file);
    return 0;  
  }

  Superflow* cutflow = new Superflow(); // initialize the cutflow
  cutflow->setAnaName("SuperflowAna");                // arbitrary
  cutflow->setAnaType(AnalysisType::Ana_2Lep);        // analysis type, passed to SusyNt ?
  cutflow->setLumi(1.0);                              // set the MC normalized to X pb-1
  cutflow->setSampleName(input_file);                 // sample name, check to make sure it's set OK
  cutflow->setRunMode(run_mode);                      // make configurable via run_mode
  cutflow->setCountWeights(true);                     // print the weighted cutflows
  cutflow->setChain(chain);
  cutflow->nttools().initTriggerTool(ChainHelper::firstFile(input_file,0.)); 

  printf("makeMiniNtuple_EWK2L\t Total events available : %lli\n",chain->GetEntries());

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

 
  *cutflow << CutName("bad muon veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passBadMuon(sl->preMuons));
  };
  
  *cutflow << CutName("pass cosmic veto") << [&](Superlink* sl) -> bool {
      //if(sl->nt->evt()->eventNumber == 2031134 || 
      //   sl->nt->evt()->eventNumber == 2187215 ||
      //   sl->nt->evt()->eventNumber == 6813818 ||
      //   sl->nt->evt()->eventNumber == 7460893 /*||*/
      //  ) {
      //  printEventInformation(sl);
      //}
      return (sl->tools->passCosmicMuon(sl->baseMuons));
  };

  *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
      //std::cout << "CUTFLOW >> Event Number: " << sl->nt->evt()->eventNumber << std::endl;
      return (sl->tools->passJetCleaning(sl->baseJets));
  };

  *cutflow << CutName("pass good vertex") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_GoodVtx);
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
  //bool pass_HLT_mu18_mu8noL1           = false; // 
  //bool pass_HLT_2e12_lhloose_L12EM10VH = false; // ** 2015 dilepton triggers **
  //bool pass_HLT_e17_lhloose_mu14       = false; //
  //bool pass_HLT_mu20_mu8noL1           = false; // 
  //bool pass_HLT_2e15_lhloose_L12EM13VH = false; // ** 2016 dilepton triggers **
  ////bool pass_HLT_e17_lhloose_mu14       = false; // -> same in 2015

  *cutflow << NewVar("HLT_mu18_mu8noL1 trigger bit"); {
      *cutflow << HFTname("pass_HLT_mu18_mu8noL1");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu18_mu8noL1"); };
      *cutflow << SaveVar();
  }  

  *cutflow << NewVar("HLT_2e12_lhloose_L12EM10VH trigger bit"); {
      *cutflow << HFTname("pass_HLT_2e12_lhloose_L12EM10VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e12_lhloose_L12EM10VH"); };
      *cutflow << SaveVar();
  }  

  *cutflow << NewVar("HLT_e17_lhloose_mu14 trigger bit"); {
      *cutflow << HFTname("pass_HLT_e17_lhloose_mu14");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_mu14"); };
      *cutflow << SaveVar();
  }  

  *cutflow << NewVar("HLT_mu20_mu8noL1 trigger bit"); {
      *cutflow << HFTname("pass_HLT_mu20_mu8noL1");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu20_mu8noL1"); };
      *cutflow << SaveVar();
  }  

  *cutflow << NewVar("HLT_2e15_lhloose_L12EM13VH trigger bit"); {
      *cutflow << HFTname("pass_HLT_2e15_lhloose_L12EM13VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e15_lhloose_L12EM13VH"); };
      *cutflow << SaveVar();
  }  

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

  //*cutflow << NewVar("lepton d0"); {
  //  *cutflow << HFTname("l_d0");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->d0);
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("lepton errD0"); {
  //  *cutflow << HFTname("l_errD0");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->errD0);
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("lepton d0sig"); {
  //  *cutflow << HFTname("l_d0sig");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->d0Sig());
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("lepton z0"); {
  //  *cutflow << HFTname("l_z0");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->z0);
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("lepton errZ0"); {
  //  *cutflow << HFTname("l_errZ0");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->errZ0);
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("lepton z0sinTheta"); {
  //  *cutflow << HFTname("l_z0sinTheta");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->z0SinTheta());
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

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

  //*cutflow << NewVar("lepton ptvarcone20"); {
  //  *cutflow << HFTname("l_ptvarcone20");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->ptvarcone20);
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

  //*cutflow << NewVar("lepton ptvarcone30"); {
  //  *cutflow << HFTname("l_ptvarcone30");
  //  *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
  //    vector<double> out;
  //    for(auto& lepton : signalLeptons) {
  //        out.push_back(lepton->ptvarcone30);
  //    }
  //    return out;
  //    };
  //  *cutflow << SaveVar();
  //}

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
  JetVector baseJets, signalJets, centralLightJets, centralBJets, forwardJets, stop2lLightJets, stop2lBJets;
  *cutflow << [&](Superlink* sl, var_void*) { 
    baseJets = *sl->baseJets; 
    signalJets = *sl->jets; 
    for(auto& jet : baseJets) {
      if(sl->tools->m_jetSelector->isCentralLight(jet))  { centralLightJets.push_back(jet); } 
      else if(sl->tools->m_jetSelector->isCentralB(jet)) { centralBJets.push_back(jet);     } 
      else if(sl->tools->m_jetSelector->isForward(jet))  { forwardJets.push_back(jet);      } 
    }
    for(auto& jet : signalJets) {
      if(!sl->tools->jetSelector().isB(jet)) { stop2lLightJets.push_back(jet); }
      else                                   { stop2lBJets.push_back(jet);     }
    }
  };

  *cutflow << NewVar("number of baseline jets"); {
    *cutflow << HFTname("nBaseJets");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseJets.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central light jets"); {
    *cutflow << HFTname("nCentralLJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCLJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central b jets"); {
    *cutflow << HFTname("nCentralBJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfCBJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of forward jets"); {
    *cutflow << HFTname("nForwardJets");
    *cutflow << [](Superlink* sl, var_int*) -> int { return sl->tools->numberOfFJets(*sl->baseJets)/*(*baseJets)*/; };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central stop2l light jets"); {
    *cutflow << HFTname("nStop2lLJets");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return (int) stop2lLightJets.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central stop2l b jets"); {
    *cutflow << HFTname("nStop2lBJets");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return (int) stop2lBJets.size(); };
    *cutflow << SaveVar();
  }

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
        if(sl->tools->m_jetSelector->isCentralLight(jet))  { flav = 1; } 
        else if(sl->tools->m_jetSelector->isCentralB(jet)) { flav = 2; } 
        else if(sl->tools->m_jetSelector->isForward(jet))  { flav = 3; } 
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
      for(auto& jet : signalJets) {
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
      deltaX = (lepton0.Pz() + lepton1.Pz())/13000.;
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

  // Super-razor variables
  double MDR, shatr, cosThetaRp1, DPB, dphi_l1_l2, gamma_r;
  double dphi_vBeta_R_vBeta_T;
  TVector3 vBeta_z, pT_CM, vBeta_T_CMtoR, vBeta_r;
  *cutflow << NewVar("Super-razor variables -- shatr"); {
    *cutflow << HFTname("shatr");
    *cutflow << [&](Superlink* sl, var_float*) -> double {
      MDR = shatr = cosThetaRp1 = DPB = dphi_l1_l2 = gamma_r = -999.0;
      dphi_vBeta_R_vBeta_T = -999.0;
      kin::superRazor(signalLeptons, sl->met, vBeta_z, pT_CM,
                      vBeta_T_CMtoR, vBeta_r, shatr, DPB,
                      dphi_l1_l2, gamma_r, dphi_vBeta_R_vBeta_T,
                      MDR, cosThetaRp1);
      return shatr;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Super-razor variables -- MDR"); {
    *cutflow << HFTname("MDR");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return MDR;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Super-razor variables -- CosThetaR+1"); {
    *cutflow << HFTname("cosThetaRp1");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return cosThetaRp1;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Super-razor variables -- DPB"); {
    *cutflow << HFTname("DPB");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return DPB;
    };
    *cutflow << SaveVar();
  } 

  // Jigsaw variables
  double MDR_jigsaw       = -999.; // "m_Delta^R
  double DPB_vSS_jigsaw   = -999.; // "delta-phi-beta-R"
  double RPT_jigsaw       = -999.; // "R_p_{T}"
  double gamInvRp1_jigsaw = -999.;

  *cutflow << [&](Superlink* /*sl*/, var_void*) -> void {
    // declare the frames
    LabRecoFrame lab("lab", "lab");
    DecayRecoFrame ss("ss", "ss");
    DecayRecoFrame s1("s1", "s1");
    DecayRecoFrame s2("s2", "s2");
    VisibleRecoFrame v1("v1", "v1");
    VisibleRecoFrame v2("v2", "v2");
    InvisibleRecoFrame i1("i1", "i1");
    InvisibleRecoFrame i2("i2", "i2");

    // connect the frames
    lab.SetChildFrame(ss);
    ss.AddChildFrame(s1);
    ss.AddChildFrame(s2);
    s1.AddChildFrame(v1);
    s1.AddChildFrame(i1);
    s2.AddChildFrame(v2);
    s2.AddChildFrame(i2);

    // check that the decay tree is connected properly
    if(!lab.InitializeTree()) {
      printf("makeMiniNtuple_EWK2L\t RestFrames::InitializeTree ERROR (\"%i\")    Unable to initialize tree from lab frame. Exitting. ",__LINE__);
      exit(1); 
    }

    // define groups
    InvisibleGroup inv("inv", "invsible group jigsaws");
    inv.AddFrame(i1);
    inv.AddFrame(i2);

    CombinatoricGroup vis("vis", "visible object jigsaws");
    vis.AddFrame(v1);
    vis.SetNElementsForFrame(v1, 1, false);
    vis.AddFrame(v2);
    vis.SetNElementsForFrame(v2, 1, false);

    SetMassInvJigsaw MinMassJigsaw("MinMass", "Invisible system mass jigsaw");
    inv.AddJigsaw(MinMassJigsaw);

    SetRapidityInvJigsaw RapidityJigsaw("RapidityJigsaw", "invisible system rapidity jigsaw");
    inv.AddJigsaw(RapidityJigsaw);
    RapidityJigsaw.AddVisibleFrames(lab.GetListVisibleFrames());

    ContraBoostInvJigsaw ContraBoostJigsaw("ContraBoostJigsaw", "ContraBoost Invariant Jigsaw");
    inv.AddJigsaw(ContraBoostJigsaw);
    ContraBoostJigsaw.AddVisibleFrames((s1.GetListVisibleFrames()), 0);
    ContraBoostJigsaw.AddVisibleFrames((s2.GetListVisibleFrames()), 1);
    ContraBoostJigsaw.AddInvisibleFrame(i1, 0);
    ContraBoostJigsaw.AddInvisibleFrame(i2, 1);

    MinMassesCombJigsaw HemiJigsaw("hemi_jigsaw", "Minimize m_{v_{1,2}} jigsaw");
    vis.AddJigsaw(HemiJigsaw);
    HemiJigsaw.AddFrame(v1, 0);
    HemiJigsaw.AddFrame(v2, 1);

    // check that the jigsaws are in place
    if(!lab.InitializeAnalysis()) {
      printf("makeMiniNtuple_EWK2L\t RestFrames::InitializeAnalysis ERROR (\"%i\")    Unable to initialize analysis from lab frame. Exitting.",__LINE__);
      exit(1);
    }

    // clear the event for sho
    lab.ClearEvent();

    // set the met
    TVector3 met3vector(met.Px(), met.Py(), met.Pz());
    inv.SetLabFrameThreeVector(met3vector);

    // add leptons to the visible group
    // leptons holds TLorentzVectors
    vis.AddLabFrameFourVector(lepton0);
    vis.AddLabFrameFourVector(lepton1);

    // analayze that
    if(!lab.AnalyzeEvent()) {
      printf("makeMiniNtuple_EWK2L\t RestFrames::AnalyzeEvent ERROR. Exitting.");
      exit(1);
    }

    /// system mass
    double shat_jigsaw = ss.GetMass();

    // RATIO OF CM pT
    TVector3 vPTT = ss.GetFourVector(lab).Vect();
    RPT_jigsaw = vPTT.Pt() / (vPTT.Pt() + shat_jigsaw / 4.);

    // shapes
    gamInvRp1_jigsaw = ss.GetVisibleShape();

    // MDR_jigsaw
    MDR_jigsaw = 2.0 * v1.GetEnergy(s1);

    // BOOST ANGLES
    DPB_vSS_jigsaw = ss.GetDeltaPhiBoostVisible();
  };

  *cutflow << NewVar("Jigsaw variables -- MDR"); {
    *cutflow << HFTname("MDR_jigsaw");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return MDR_jigsaw;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Jigsaw variables -- RPT"); {
    *cutflow << HFTname("RPT_jigsaw");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return RPT_jigsaw;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Jigsaw variables -- gamInvRp1"); {
    *cutflow << HFTname("gamInvRp1_jigsaw");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return gamInvRp1_jigsaw;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("Jigsaw variables -- DPB_vSS"); {
    *cutflow << HFTname("DPB_vSS_jigsaw");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double {
      return DPB_vSS_jigsaw;
    };
    *cutflow << SaveVar();
  }

  // Clear 
  *cutflow << [&](Superlink* /*sl*/, var_void*) { 
    baseLeptons.clear();signalLeptons.clear(); 
    baseJets.clear(); signalJets.clear(); centralLightJets.clear(); centralBJets.clear(); forwardJets.clear(); stop2lLightJets.clear(); stop2lBJets.clear();
    meff=0.,R1=0.,R2=0.,deltaX=0.,mT2=0.,cthllb=0.,MDR_jigsaw=0.,RPT_jigsaw=0.,gamInvRp1_jigsaw=0.,DPB_vSS_jigsaw=0.;
  };

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods end here
  ///////////////////////////////////////////////////////////////////////

  // Initialize the cutflow and start the event loop.
  chain->Process(cutflow, input_file, n_events, n_skip_events);

  delete cutflow;
  delete chain;

  // Print information
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
  printf("makeMiniNtuple_EWK2L\t All done!\n");
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
  return 0; 
}
