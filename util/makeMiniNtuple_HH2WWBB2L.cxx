//////////////////////////////////////////////////////////////////////////////////////////
// \project:    ATLAS Experiment at CERN's LHC
// \package:    SusyAnalysis
// \file:       $Id$
// \author:     Alaettin.Serhan.Mete@cern.ch
// \history:    N/A 
// 
// Copyright (C) 2018 University of California, Irvine
//////////////////////////////////////////////////////////////////////////////////////////

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
// Main function
///////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

  ///////////////////////////////////////////////////////////////////////
  // Read user inputs - NOT super safe so be careful :)
  unsigned int n_events   = -1;
  unsigned int dilep_flav = 0;
  unsigned int n_skip_events  = 0;
  char *input_file  = nullptr;
  char *name_suffix = nullptr;
  SuperflowRunMode run_mode = SuperflowRunMode::nominal; // SuperflowRunMode::all_syst; //SuperflowRunMode::nominal;
  int c;

  opterr = 0;
  while ((c = getopt (argc, argv, "f:s:n:l:h")) != -1)
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
      case 'l':
        dilep_flav = atoi(optarg);
        break;
      case 'h':
        usage("makeMiniNtuple_HH2WWBB2L");
        return 1;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "makeMiniNtuple_HH2WWBB2L\t Option -%c requires an argument.\n", optopt);
        else if (optopt == 'n')
          fprintf (stderr, "makeMiniNtuple_HH2WWBB2L\t Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "makeMiniNtuple_HH2WWBB2L\t Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "makeMiniNtuple_HH2WWBB2L\t Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  // Catch problems or cast
  for (int index = optind; index < argc; index++)
    printf ("makeMiniNtuple_HH2WWBB2L\t Non-option argument %s\n", argv[index]);
  if (input_file==nullptr) {
    printf("makeMiniNtuple_HH2WWBB2L\t An input file must be provided with option -f (a list, a DIR or single file)\n");
    return 0;  
  }
  

  // Print information
  printf("makeMiniNtuple_HH2WWBB2L\t =================================================================\n");
  printf("makeMiniNtuple_HH2WWBB2L\t Running SusyAnalysis/makeMiniNtuple_HH2WWBB2L\n");
  printf("makeMiniNtuple_HH2WWBB2L\t =================================================================\n");
  printf("makeMiniNtuple_HH2WWBB2L\t   Flags:\n");
  printf("makeMiniNtuple_HH2WWBB2L\t     Input file (-f)         : %s\n",input_file);
  printf("makeMiniNtuple_HH2WWBB2L\t     Output file suffix (-s) : %s\n",name_suffix);
  printf("makeMiniNtuple_HH2WWBB2L\t     Number of events (-n)   : %i\n",n_events);
  printf("makeMiniNtuple_HH2WWBB2L\t     Dilepton flavor  (-l)   : %i\n",dilep_flav); // 0: ee - 1: mm - 2: em+me
  printf("makeMiniNtuple_HH2WWBB2L\t =================================================================\n");
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
    printf("makeMiniNtuple_HH2WWBB2L\t Cannot understand input %s",input_file);
    return 0;  
  }

  Superflow* cutflow = new Superflow(); // initialize the cutflow
  cutflow->setAnaName("SuperflowAna");                // arbitrary
  cutflow->setAnaType(AnalysisType::Ana_WWBB);        // analysis type, passed to SusyNt ?
  cutflow->setLumi(1.0);                              // set the MC normalized to X pb-1
  cutflow->setSampleName(input_file);                 // sample name, check to make sure it's set OK
  cutflow->setRunMode(run_mode);                      // make configurable via run_mode
  cutflow->setCountWeights(true);                     // print the weighted cutflows
  if(name_suffix != nullptr)
    cutflow->setFileSuffix(name_suffix);              // setting a suffix to the file name    
  cutflow->setChain(chain);
  cutflow->nttools().initTriggerTool(ChainHelper::firstFile(input_file,0.)); 

  printf("makeMiniNtuple_HH2WWBB2L\t Total events available : %lli\n",chain->GetEntries());

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

  *cutflow << CutName("SCT error") << [&](Superlink* sl) -> bool {
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

  // Triggers
  *cutflow << CutName("trigger") << [&](Superlink* sl) -> bool {
      if (dilep_flav == 0)
        return (sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_2e17_lhvloose_nod0"));
      else if (dilep_flav == 1)
        return (sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_mu22_mu8noL1"));
      else if (dilep_flav == 2)
        return (sl->tools->triggerTool().passTrigger(sl->nt->evt()->trigBits, "HLT_e17_lhloose_nod0_mu14"));
  };

  // Leptons
  *cutflow << CutName("at least two baseline leptons") << [&](Superlink* sl) -> bool {
      if (dilep_flav == 0)
        return (sl->baseLeptons->size() >= 2 && sl->baseLeptons->at(0)->isEle() && sl->baseLeptons->at(1)->isEle());
      else if (dilep_flav == 1)
        return (sl->baseLeptons->size() >= 2 && !sl->baseLeptons->at(0)->isEle() && !sl->baseLeptons->at(1)->isEle());
      else if (dilep_flav == 2)
        return (sl->baseLeptons->size() >= 2 && (sl->baseLeptons->at(0)->isEle() != sl->baseLeptons->at(1)->isEle()));
  };

  *cutflow << CutName("exactly two signal leptons") << [](Superlink* sl) -> bool {
      return (sl->leptons->size() == 2);
  };

  *cutflow << CutName("correct dilepton flavor") << [&](Superlink* sl) -> bool {
      if (dilep_flav == 0)
        return (sl->leptons->at(0)->isEle() && sl->leptons->at(1)->isEle());
      else if (dilep_flav == 1)
        return (!sl->leptons->at(0)->isEle() && !sl->leptons->at(1)->isEle());
      else if (dilep_flav == 2)
        return ((sl->leptons->at(0)->isEle() != sl->leptons->at(1)->isEle()));
  };

  *cutflow << CutName("opposite sign") << [](Superlink* sl) -> bool {
      return ((sl->leptons->at(0)->q * sl->leptons->at(1)->q) < 0);
  };

  *cutflow << CutName("mll > 20 GeV") << [](Superlink* sl) -> bool {
      return ((*sl->leptons->at(0) + *sl->leptons->at(1)).M() > 20.);
  };

  *cutflow << CutName("z-veto") << [](Superlink* sl) -> bool {
      TLorentzVector lepton0 = *(*sl->leptons).at(0);
      TLorentzVector lepton1 = *(*sl->leptons).at(1);
      bool isSF =  (*sl->leptons).at(0)->isEle() == (*sl->leptons).at(1)->isEle() ? true : false;
      return ((isSF && fabs((lepton0+lepton1).M() - 91.2) > 20.) || !isSF);
  };

  *cutflow << CutName("subleading lepton pt > 20 GeV") << [](Superlink* sl) -> bool {
      return (sl->leptons->at(1)->Pt() > 20.);
  };

  JetVector bjets;
  *cutflow << CutName("exactly two b-jets") << [&](Superlink* sl) -> bool {
      bjets.clear();
      for(auto& jet : *sl->jets) {
        if(sl->tools->jetSelector().isB(jet))  { bjets.push_back(jet); } 
      }
      return (bjets.size()==2);
  };
   
  *cutflow << CutName("mbb in (100, 140) GeV") << [&](Superlink* /*sl*/) -> bool {
      return ((*bjets.at(0) + *bjets.at(1)).M() > 100. && (*bjets.at(0) + *bjets.at(1)).M() < 140.);
  };

  *cutflow << CutName("mT2 in (100, 140) GeV") << [&](Superlink* sl) -> bool {
      TLorentzVector lepton0 = *(*sl->leptons).at(0);
      TLorentzVector lepton1 = *(*sl->leptons).at(1);
      TLorentzVector bjet0   = *bjets.at(0);
      TLorentzVector bjet1   = *bjets.at(1);
      TLorentzVector met;
      met.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
                     sl->met->Et*sin(sl->met->phi),
                     0.,
                     sl->met->Et);
      ComputeMT2 mycalc = ComputeMT2(lepton0+lepton1,bjet0+bjet1,met,0.,0.); // masses 0. 0.
      double mT2 = mycalc.Compute();
      return (mT2 > 100. && mT2 < 140.);
  };

  *cutflow << CutName("dRll < 0.9") << [](Superlink* sl) -> bool {
      TLorentzVector lepton0 = *(*sl->leptons).at(0);
      TLorentzVector lepton1 = *(*sl->leptons).at(1);
      return (lepton0.DeltaR(lepton1) < 0.9);
  };

  *cutflow << CutName("HT2R > 0.8") << [&](Superlink* sl) -> bool {
      TLorentzVector lepton0 = *(*sl->leptons).at(0);
      TLorentzVector lepton1 = *(*sl->leptons).at(1);
      TLorentzVector bjet0   = *bjets.at(0);
      TLorentzVector bjet1   = *bjets.at(1);
      TLorentzVector met;
      met.SetPxPyPzE(sl->met->Et*cos(sl->met->phi),
                     sl->met->Et*sin(sl->met->phi),
                     0.,
                     sl->met->Et);
      double HT2R = ((met+lepton0+lepton1).Pt()+(bjet0+bjet1).Pt())/(met.Pt()+lepton0.Pt()+lepton1.Pt()+bjet0.Pt()+bjet1.Pt());
      return (HT2R > 0.8);
  };

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods end here
  ///////////////////////////////////////////////////////////////////////

  // Initialize the cutflow and start the event loop.
  chain->Process(cutflow, input_file, n_events, n_skip_events);

  delete cutflow;
  delete chain;

  // Print information
  printf("makeMiniNtuple_HH2WWBB2L\t =================================================================\n");
  printf("makeMiniNtuple_HH2WWBB2L\t All done!\n");
  printf("makeMiniNtuple_HH2WWBB2L\t =================================================================\n");
  return 0; 
}
