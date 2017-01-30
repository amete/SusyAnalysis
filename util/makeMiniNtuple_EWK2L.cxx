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
  char *input_file  = nullptr;
  char *name_suffix = nullptr;
  SuperflowRunMode run_mode = SuperflowRunMode::all_syst; // SuperflowRunMode::all_syst; //SuperflowRunMode::nominal;
  int c;

  opterr = 0;
  while ((c = getopt (argc, argv, "f:s:n:h")) != -1)
    switch (c)
      {
      case 'f':
        input_file = optarg;
        break;
      case 's':
        name_suffix = optarg;
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
  if (input_file==nullptr) {
    printf("makeMiniNtuple_EWK2L\t An input file must be provided with option -f (a list, a DIR or single file)\n");
    return 0;  
  }
  

  // Print information
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
  printf("makeMiniNtuple_EWK2L\t Running SusyAnalysis/makeMiniNtuple_EWK2L\n");
  printf("makeMiniNtuple_EWK2L\t =================================================================\n");
  printf("makeMiniNtuple_EWK2L\t   Flags:\n");
  printf("makeMiniNtuple_EWK2L\t     Input file (-f)         : %s\n",input_file);
  printf("makeMiniNtuple_EWK2L\t     Output file suffix (-s) : %s\n",name_suffix);
  printf("makeMiniNtuple_EWK2L\t     Number of events (-n)   : %i\n",n_events );
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
  if(name_suffix != nullptr)
    cutflow->setFileSuffix(name_suffix);              // setting a suffix to the file name    
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
      cutflags = sl->nt->evt()->cutFlags[NtSys::NOM];
      return (sl->tools->passGRL(cutflags));
  };

  *cutflow << CutName("LAr error") << [&](Superlink* sl) -> bool {
      return (sl->tools->passLarErr(cutflags));
  };
  
  *cutflow << CutName("Tile error") << [&](Superlink* sl) -> bool {
      return (sl->tools->passTileErr(cutflags));
  };
  
  *cutflow << CutName("TTC veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passTTC(cutflags));
  };

  *cutflow << CutName("SCT seu") << [&](Superlink* sl) -> bool {
      return (sl->tools->passSCTErr(cutflags));
  };
 
  *cutflow << CutName("pass good vertex") << [&](Superlink* sl) -> bool {
      return (sl->tools->passGoodVtx(cutflags));
  };

  *cutflow << CutName("bad muon veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passBadMuon(sl->preMuons));
  };
  
  *cutflow << CutName("pass cosmic veto") << [&](Superlink* sl) -> bool {
      return (sl->tools->passCosmicMuon(sl->baseMuons));
  };

  *cutflow << CutName("jet cleaning") << [&](Superlink* sl) -> bool {
      return (sl->tools->passJetCleaning(sl->baseJets));
  };
   
  //// MET Cleaning - temporary
  //#warning "MET CLEANING IS TURNED ON!!!"
  //*cutflow << CutName("met cleaning") << [&](Superlink* sl) -> bool {
  //    return (sl->tools->passMetCleaning(sl->met));
  //};

  //  Analysis Cuts
  //*cutflow << CutName("at least one baseline jet") << [](Superlink* sl) -> bool {
  //  return (sl->baseJets->size()>=1);
  //};

  //*cutflow << CutName("at least one signal jet") << [](Superlink* sl) -> bool {
  //  int nSignalJets = 0;
  //  for(auto& jet : *sl->baseJets) {
  //    if(fabs(jet->Eta())<4.5)  { nSignalJets++; } 
  //  }
  //  return (nSignalJets>=1);
  //};

  //*cutflow << CutName("at least one b jet") << [](Superlink* sl) -> bool {
  //  int nSignalJets = 0;
  //  for(auto& jet : *sl->baseJets) {
  //    if(sl->tools->m_jetSelector->isB(jet))  { nSignalJets++; } 
  //  }
  //  return (nSignalJets>=1);
  //};

  *cutflow << CutName("exactly two baseline leptons") << [](Superlink* sl) -> bool {
      //if(sl->nt->evt()->eventNumber == 150148 || 
      //   sl->nt->evt()->eventNumber == 154664 //||
      //   sl->nt->evt()->eventNumber == 150393 ||
      //   sl->nt->evt()->eventNumber == 150396 //||
      //  ) {
      //  printEventInformation(sl);
      //}
      return (sl->baseLeptons->size() == 2);
  };

  *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
      //std::cout << "SERHAN " << sl->nt->evt()->eventNumber; // << std::endl;
      //printf(" event %*.4f pile-up %*.4f ",10,sl->nt->evt()->w,10,sl->nt->evt()->wPileup);
      //sl->weights->product(); 
      return (sl->leptons->size() == 2);
  };

  *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
      return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
  };

  *cutflow << CutName("lepton pTs > 20 GeV") << [](Superlink* sl) -> bool {
      return ((sl->leptons->at(0)->Pt()>20) && (sl->leptons->at(1)->Pt()>20));
  };

  *cutflow << CutName("mll > 40 GeV") << [](Superlink* sl) -> bool {
      return ((*sl->leptons->at(0) + *sl->leptons->at(1)).M() > 40.);
  };

  //*cutflow << CutName("loose jet-veto") << [](Superlink* sl) -> bool {
  //  bool veto = false;
  //  for(auto& jet : *sl->baseJets) {
  //    if(sl->tools->m_jetSelector->isCentralLight(jet))  { veto = jet->Pt()>50. ? true : false; if(veto) break; } 
  //    else if(sl->tools->m_jetSelector->isForward(jet))  { veto = true; break; } 
  //  }
  //  return !veto;
  //};
  //*cutflow << CutName("b-veto") << [](Superlink* sl) -> bool {
  //  int nSignalJets = 0;
  //  for(auto& jet : *sl->baseJets) {
  //    if(sl->tools->m_jetSelector->isB(jet))  { nSignalJets++; } 
  //  }
  //  return (nSignalJets==0);
  //};

  //*cutflow << CutName("MET > 50 GeV") << [](Superlink* sl) -> bool {
  //    return (sl->met->Et > 50.0);
  //};

  //*cutflow << CutName("z-veto") << [](Superlink* sl) -> bool {
  //    TLorentzVector lepton0 = *(*sl->leptons).at(0);
  //    TLorentzVector lepton1 = *(*sl->leptons).at(1);
  //    bool isSF =  (*sl->leptons).at(0)->isEle() == (*sl->leptons).at(1)->isEle() ? true : false;
  //    return ((isSF && fabs((lepton0+lepton1).M() - 91.2) > 10.) || !isSF);
  //};

  //*cutflow << CutName("mT2 > 90 GeV") << [](Superlink* sl) -> bool {
  //    TLorentzVector lepton0 = *(*sl->leptons).at(0);
  //    TLorentzVector lepton1 = *(*sl->leptons).at(1);
  //    TLorentzVector met;
  //    met.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
  //                   sl->met->Et*sin(sl->met->phi),
  //                   0.,
  //                   sl->met->Et);
  //    ComputeMT2 mycalc = ComputeMT2(lepton0,lepton1,met,0.,0.); // masses 0. 0.
  //    double mT2 = mycalc.Compute();
  //    return (mT2 > 90.0);
  //};

  //  Output Ntuple Setup
  //      > Ntuple variables

  // Event variables
  //bool pass_HLT_mu18_mu8noL1            = false; // 
  //bool pass_HLT_2e12_lhloose_L12EM10VH  = false; // ** 2015 dilepton triggers **
  //bool pass_HLT_e17_lhloose_mu14        = false; //
  //bool pass_HLT_mu20_mu8noL1            = false; // 
  //bool pass_HLT_2e15_lhvloose_L12EM13VH = false; // ** 2016 dilepton triggers **
  //bool pass_HLT_e17_lhloose_mu14        = false; // -> same in 2015
  //bool pass_HLT_2e17_lhvloose_nod0      = false;
  //bool pass_HLT_mu22_mu8noL1            = false;
  //bool pass_HLT_e17_lhloose_nod0_mu14   = false;

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

  *cutflow << NewVar("HLT_2e15_lhvloose_L12EM13VH trigger bit"); {
      *cutflow << HFTname("pass_HLT_2e15_lhvloose_L12EM13VH");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e15_lhvloose_L12EM13VH"); };
      *cutflow << SaveVar();
  } 

  *cutflow << NewVar("HLT_2e17_lhvloose_nod0 trigger bit"); {
      *cutflow << HFTname("pass_HLT_2e17_lhvloose_nod0");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0"); };
      *cutflow << SaveVar();
  } 

  *cutflow << NewVar("HLT_mu22_mu8noL1 trigger bit"); {
      *cutflow << HFTname("pass_HLT_mu22_mu8noL1");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1"); };
      *cutflow << SaveVar();
  } 

  *cutflow << NewVar("HLT_e17_lhloose_nod0_mu14 trigger bit"); {
      *cutflow << HFTname("pass_HLT_e17_lhloose_nod0_mu14");
      *cutflow << [](Superlink* sl, var_bool*) -> bool { return sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14"); };
      *cutflow << SaveVar();
  } 

  *cutflow << NewVar("treatAsYear"); { 
      *cutflow << HFTname("treatAsYear");
      *cutflow << [](Superlink* sl, var_double*) -> int { 
          return sl->nt->evt()->treatAsYear;
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("event weight"); {
      *cutflow << HFTname("eventweight");
      *cutflow << [](Superlink* sl, var_double*) -> double { 
          return sl->weights->product() * sl->nt->evt()->wPileup;
          //return sl->weights->product();
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("event weight no pileup"); {
      *cutflow << HFTname("eventweight_nopileup");
      *cutflow << [](Superlink* sl, var_double*) -> double { 
          return sl->weights->product();
      };
      *cutflow << SaveVar();
  }

  *cutflow << NewVar("event weight pileup"); {
      *cutflow << HFTname("eventweight_onlypileup");
      *cutflow << [](Superlink* sl, var_double*) -> double { 
          return sl->nt->evt()->wPileup;
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

  *cutflow << NewVar("lepton type"); {
    *cutflow << HFTname("l_type");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->mcType);
      }
      return out;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("lepton origin"); {
    *cutflow << HFTname("l_origin");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> out;
      for(auto& lepton : signalLeptons) {
        out.push_back(lepton->mcOrigin);
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
  JetVector baseJets, signalJets, 
            centralLightJets, centralLightJets30, centralLightJets40, centralLightJets50, centralLightJets60, centralLightJets70,
            centralBJets,  
            forwardJets, forwardJets30, forwardJets40, forwardJets50,
            stop2lLightJets, stop2lBJets;
  *cutflow << [&](Superlink* sl, var_void*) { 
    baseJets = *sl->baseJets; 
    signalJets = *sl->jets; 
    for(auto& jet : baseJets) {
      if(sl->tools->m_jetSelector->isCentralLight(jet))  { centralLightJets.push_back(jet); 
                                                           if(jet->Pt()>30.) { centralLightJets30.push_back(jet); } 
                                                           if(jet->Pt()>40.) { centralLightJets40.push_back(jet); } 
                                                           if(jet->Pt()>50.) { centralLightJets50.push_back(jet); } 
                                                           if(jet->Pt()>60.) { centralLightJets60.push_back(jet); } 
                                                           if(jet->Pt()>70.) { centralLightJets70.push_back(jet); } 
                                                         } 
      else if(sl->tools->m_jetSelector->isCentralB(jet)) { centralBJets.push_back(jet);     } 
      else if(sl->tools->m_jetSelector->isForward(jet))  { forwardJets.push_back(jet);      
                                                           if(jet->Pt()>30.) { forwardJets30.push_back(jet); } 
                                                           if(jet->Pt()>40.) { forwardJets40.push_back(jet); } 
                                                           if(jet->Pt()>50.) { forwardJets50.push_back(jet); } 
                                                         } 
    }
    std::sort(centralLightJets.begin()  , centralLightJets.end()  , comparePt);
    std::sort(centralLightJets30.begin(), centralLightJets30.end(), comparePt);
    std::sort(centralLightJets40.begin(), centralLightJets40.end(), comparePt);
    std::sort(centralLightJets50.begin(), centralLightJets50.end(), comparePt);
    std::sort(centralLightJets60.begin(), centralLightJets60.end(), comparePt);
    std::sort(centralLightJets70.begin(), centralLightJets70.end(), comparePt);
    std::sort(centralBJets.begin()      , centralBJets.end()      , comparePt);
    std::sort(forwardJets.begin()       , forwardJets.end()       , comparePt);
    std::sort(forwardJets30.begin()     , forwardJets30.end()     , comparePt);
    std::sort(forwardJets40.begin()     , forwardJets40.end()     , comparePt);
    std::sort(forwardJets50.begin()     , forwardJets50.end()     , comparePt);
    for(auto& jet : signalJets) {
      if(!sl->tools->jetSelector().isB(jet)) { stop2lLightJets.push_back(jet); }
      else                                   { stop2lBJets.push_back(jet);     }
    }
  };

  // Lowest pT jets
  *cutflow << NewVar("pt of softest centralLightJet"); {
    *cutflow << HFTname("softestCentralLJetPt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> float { 
        float value = centralLightJets.size() > 0 ? (centralLightJets.at(centralLightJets.size()-1))->Pt() : 0. ; 
        return value;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("pt of softest ForwardJet"); {
    *cutflow << HFTname("softestForwardJetPt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> float { 
        float value = forwardJets.size() > 0 ? (forwardJets.at(forwardJets.size()-1))->Pt() : 0.;
        return value;
    };
    *cutflow << SaveVar();
  }

  // Highest pT jets
  *cutflow << NewVar("pt of hardest centralLightJet"); {
    *cutflow << HFTname("hardestCentralLJetPt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> float { 
        float value = centralLightJets.size() > 0 ? (centralLightJets.at(0))->Pt() : 0. ; 
        return value;
    };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("pt of hardest ForwardJet"); {
    *cutflow << HFTname("hardestForwardJetPt");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> float { 
        float value = forwardJets.size() > 0 ? (forwardJets.at(0))->Pt() : 0.;
        return value;
    };
    *cutflow << SaveVar();
  }

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

  *cutflow << NewVar("number of central light jets pt 30"); {
    *cutflow << HFTname("nCentralLJets30");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return centralLightJets30.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central light jets pt 40"); {
    *cutflow << HFTname("nCentralLJets40");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return centralLightJets40.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central light jets pt 50"); {
    *cutflow << HFTname("nCentralLJets50");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return centralLightJets50.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central light jets pt 60"); {
    *cutflow << HFTname("nCentralLJets60");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return centralLightJets60.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of central light jets pt 70"); {
    *cutflow << HFTname("nCentralLJets70");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return centralLightJets70.size(); };
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

  *cutflow << NewVar("number of forward jets pt 30 GeV"); {
    *cutflow << HFTname("nForwardJets30");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return forwardJets30.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of forward jets pt 40 GeV"); {
    *cutflow << HFTname("nForwardJets40");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return forwardJets40.size(); };
    *cutflow << SaveVar();
  }

  *cutflow << NewVar("number of forward jets 50 GeV"); {
    *cutflow << HFTname("nForwardJets50");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return forwardJets50.size(); };
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

  // DeltaPhi(ll,met)
  *cutflow << NewVar("deltaPhi(ll,met)"); {
    *cutflow << HFTname("dphi_ll_met");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
        return dileptonP4.DeltaPhi(met);
    };
    *cutflow << SaveVar();
  }

  // DeltaPhi(leading light jet,met)
  *cutflow << NewVar("deltaPhi(leading light jet,met)"); {
    *cutflow << HFTname("dphi_ljet_met");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
        double result = -999.;
        if(centralLightJets.size()>0) {
          TLorentzVector leadingLightJet = *centralLightJets.at(0);
          result = leadingLightJet.DeltaPhi(met);
        }
        return result;
    };
    *cutflow << SaveVar();
  }

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

  // S1
  *cutflow << NewVar("lep_pt/(jet_pt + met)"); {
    *cutflow << HFTname("S1");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      // Leptons
      TLorentzVector leps(0.,0.,0.,0.); TLorentzVector jets(0.,0.,0.,0.);
      for(auto& lepton : signalLeptons) {
        leps += *lepton;
      }
      // Jets
      for(auto& jet : signalJets) {
        jets += *jet;
      }
      double num   = leps.Pt();
      double denom = jets.Pt()+met.Pt();
      return denom > 0. ? num/denom : -999.;
    };
    *cutflow << SaveVar();
  }

  // S2
  *cutflow << NewVar("met/jet_pt"); {
    *cutflow << HFTname("S2");
    *cutflow << [&](Superlink* /*sl*/, var_float*) -> double { 
      TLorentzVector jets(0.,0.,0.,0.);
      // Jets
      for(auto& jet : signalJets) {
        jets += *jet;
      }
      double num   = met.Pt();
      double denom = jets.Pt();
      return denom > 0. ? num/denom : -999.;
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

  // metRel
  double metRel = 0.;
  *cutflow << NewVar("met rel."); {
    *cutflow << HFTname("metRel");
    *cutflow << [&](Superlink* sl, var_float*) -> double { 
      metRel = kin::getMetRel(sl->met,signalLeptons,signalJets);
      return metRel;
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

  // CUTFLOWWWW
  *cutflow << [&](Superlink* sl, var_void*) {
      std::cout <<  "EVENT :: " <<  sl->nt->evt()->eventNumber  << " " << 
                    centralLightJets.size()     << " " <<
                    centralBJets.size()         << " " <<
                    forwardJets30.size()        << " " <<
                    mT2                         << " "
                    << std::endl;
  };

  // Clear 
  *cutflow << [&](Superlink* /*sl*/, var_void*) { 
    baseLeptons.clear();signalLeptons.clear(); 
    baseJets.clear(); signalJets.clear(); centralLightJets.clear(); centralLightJets30.clear(); centralLightJets40.clear(); centralLightJets50.clear(); centralLightJets60.clear(); centralLightJets70.clear();  
    centralBJets.clear(); forwardJets.clear(); forwardJets30.clear(); forwardJets40.clear(); forwardJets50.clear(); stop2lLightJets.clear(); stop2lBJets.clear();
    meff=0.,R1=0.,R2=0.,deltaX=0.,mT2=0.,cthllb=0.,MDR_jigsaw=0.,RPT_jigsaw=0.,gamInvRp1_jigsaw=0.,DPB_vSS_jigsaw=0.,metRel=0.;
  };

  ///////////////////////////////////////////////////////////////////////
  // Systematics
  ///////////////////////////////////////////////////////////////////////

  ////////////////////////////////////
  // weight systematics
  ////////////////////////////////////

  // electron eff
  *cutflow << NewSystematic("shift in electron ID efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ID_UP, SupersysWeight::EL_EFF_ID_DN);
      *cutflow << TreeName("EL_EFF_ID");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in electron ISO efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::EL_EFF_ISO_UP, SupersysWeight::EL_EFF_ISO_DN);
      *cutflow << TreeName("EL_EFF_Iso");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in electron RECO efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::EL_EFF_RECO_UP, SupersysWeight::EL_EFF_RECO_DN);
      *cutflow << TreeName("EL_EFF_Reco");
      *cutflow << SaveSystematic();
  }

  // muon eff
  *cutflow << NewSystematic("muon eff stat uncertainty"); {
      *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_UP, SupersysWeight::MUON_EFF_STAT_DN);
      *cutflow << TreeName("MUON_EFF_STAT");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon eff stat uncertainty (low pt)"); {
      *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_STAT_LOWPT_UP, SupersysWeight::MUON_EFF_STAT_LOWPT_DN);
      *cutflow << TreeName("MUON_EFF_STAT_LOWPT");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon eff syst uncertainty"); {
      *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYS_UP, SupersysWeight::MUON_EFF_SYS_DN);
      *cutflow << TreeName("MUON_EFF_SYS");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon eff syst uncertainty (low pt"); {
      *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_SYS_LOWPT_UP, SupersysWeight::MUON_EFF_SYS_LOWPT_DN);
      *cutflow << TreeName("MUON_EFF_SYS_LOWPT");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon eff iso stat uncertainty"); {
      *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_STAT_UP, SupersysWeight::MUON_EFF_ISO_STAT_DN);
      *cutflow << TreeName("MUON_ISO_STAT");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon eff iso syst uncertainty"); {
      *cutflow << WeightSystematic(SupersysWeight::MUON_EFF_ISO_SYS_UP, SupersysWeight::MUON_EFF_ISO_SYS_DN);
      *cutflow << TreeName("MUON_ISO_SYS");
      *cutflow << SaveSystematic();
  }

  // flavor tagging eff
  *cutflow << NewSystematic("shift in b-tag efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::FT_EFF_B_UP, SupersysWeight::FT_EFF_B_DN);
      *cutflow << TreeName("FT_EFF_B");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in c-tag efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::FT_EFF_C_UP, SupersysWeight::FT_EFF_C_DN);
      *cutflow << TreeName("FT_EFF_C");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in light tag (i.e. mis-tag) efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::FT_EFF_LT_UP, SupersysWeight::FT_EFF_LT_DN);
      *cutflow << TreeName("FT_EFF_Light");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift flavor tagging extrapolation ?"); {
      *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAP_UP, SupersysWeight::FT_EFF_EXTRAP_DN);
      *cutflow << TreeName("FT_EFF_extrapolation");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift flavor tagging extrapolation - charm ?"); {
      *cutflow << WeightSystematic(SupersysWeight::FT_EFF_EXTRAPC_UP, SupersysWeight::FT_EFF_EXTRAPC_DN);
      *cutflow << TreeName("FT_EFF_extrapolation_charm");
      *cutflow << SaveSystematic();
  }

  // jvt eff
  *cutflow << NewSystematic("shift in JVT efficiency"); {
      *cutflow << WeightSystematic(SupersysWeight::JVT_EFF_UP, SupersysWeight::JVT_EFF_DN);
      *cutflow << TreeName("JET_JVTEff");
      *cutflow << SaveSystematic();
  }

  // pileup
  *cutflow << NewSystematic("shift in data mu (pile-up)"); {
      *cutflow << WeightSystematic(SupersysWeight::PILEUP_UP, SupersysWeight::PILEUP_DN);
      *cutflow << TreeName("PILEUP");
      *cutflow << SaveSystematic();
  }


  ////////////////////////////////////
  // shape systematics
  ////////////////////////////////////

  // egamma
  *cutflow << NewSystematic("shift in e-gamma resolution (UP)"); {
      *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_UP);
      *cutflow << TreeName("EG_RESOLUTION_ALL_UP");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in e-gamma resolution (DOWN)"); {
      *cutflow << EventSystematic(NtSys::EG_RESOLUTION_ALL_DN);
      *cutflow << TreeName("EG_RESOLUTION_ALL_DN");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in e-gamma scale (UP)"); {
      *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_UP);
      *cutflow << TreeName("EG_SCALE_ALL_UP");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("shift in e-gamma scale (DOWN)"); {
      *cutflow << EventSystematic(NtSys::EG_SCALE_ALL_DN);
      *cutflow << TreeName("EG_SCALE_ALL_DN");
      *cutflow << SaveSystematic();
  }
  // muon
  *cutflow << NewSystematic("muon ID (UP)"); {
      *cutflow << EventSystematic(NtSys::MUONS_ID_UP);
      *cutflow << TreeName("MUONS_ID_UP");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon ID (DOWN)"); {
      *cutflow << EventSystematic(NtSys::MUONS_ID_DN);
      *cutflow << TreeName("MUONS_ID_DN");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon MS (UP)"); {
      *cutflow << EventSystematic(NtSys::MUONS_MS_UP);
      *cutflow << TreeName("MUONS_MS_UP");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon MS (DOWN)"); {
      *cutflow << EventSystematic(NtSys::MUONS_MS_DN);
      *cutflow << TreeName("MUONS_MS_DN");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon scale shift (UP)"); {
      *cutflow << EventSystematic(NtSys::MUONS_SCALE_UP);
      *cutflow << TreeName("MUONS_SCALE_UP");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("muon scale shift (DN)"); {
      *cutflow << EventSystematic(NtSys::MUONS_SCALE_DN);
      *cutflow << TreeName("MUONS_SCALE_DN");
      *cutflow << SaveSystematic();
  }

  // jet
  *cutflow << NewSystematic("JER"); {
      *cutflow << EventSystematic(NtSys::JER);
      *cutflow << TreeName("JER");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("JES NP set 1 (up)"); {
      *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_UP);
      *cutflow << TreeName("JET_GroupedNP_1_UP");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("JES NP set 1 (down)"); {
      *cutflow << EventSystematic(NtSys::JET_GroupedNP_1_DN);
      *cutflow << TreeName("JET_GroupedNP_1_DN");
      *cutflow << SaveSystematic();
  }

  // met
  *cutflow << NewSystematic("MET TST Soft-Term resolution (parallel)"); {
      *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPara);
      *cutflow << TreeName("MET_SoftTrk_ResoPara");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("MET TST Soft-Term resolution (perpendicular)"); {
      *cutflow << EventSystematic(NtSys::MET_SoftTrk_ResoPerp);
      *cutflow << TreeName("MET_SoftTrk_ResoPerp");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("MET TST Soft-Term shift in scale (UP)"); {
      *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleUp);
      *cutflow << TreeName("MET_SoftTrk_ScaleUp");
      *cutflow << SaveSystematic();
  }
  *cutflow << NewSystematic("MET TST Soft-Term shift in scale (DOWN)"); {
      *cutflow << EventSystematic(NtSys::MET_SoftTrk_ScaleDown);
      *cutflow << TreeName("MET_SoftTrk_ScaleDown");
      *cutflow << SaveSystematic();
  }
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
