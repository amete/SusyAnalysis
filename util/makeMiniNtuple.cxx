// std include(s)
//#include <cstdlib>  // << for atoi
//#include "unistd.h" // << for getopt

// analysis include(s)
#include "Superflow/Superflow.h"     
#include "Superflow/Superlink.h"
#include "SusyNtuple/ChainHelper.h"
#include "SusyNtuple/string_utils.h"

using namespace sflow;

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
  while ((c = getopt (argc, argv, "f:n:")) != -1)
    switch (c)
      {
      case 'f':
        input_file = optarg;
        break;
      case 'n':
        n_events = atoi(optarg);
        break;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "makeMiniNtuple \t Option -%c requires an argument.\n", optopt);
        else if (optopt == 'n')
          fprintf (stderr, "makeMiniNtuple \t Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "makeMiniNtuple \t Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                   "makeMiniNtuple \t Unknown option character `\\x%x'.\n",
                   optopt);
        return 1;
      default:
        abort ();
      }

  // Catch problems or cast
  for (int index = optind; index < argc; index++)
    printf ("makeMiniNtuple \t Non-option argument %s\n", argv[index]);
  if (input_file==NULL) {
    printf("makeMiniNtuple \t An input file must be provided with option -f (a list, a DIR or single file)\n");
    return 0;  
  }
  

  // Print information
  printf("makeMiniNtuple \t =================================================================\n");
  printf("makeMiniNtuple \t Running SusyAnalysis/makeMiniNtuple\n");
  printf("makeMiniNtuple \t =================================================================\n");
  printf("makeMiniNtuple \t   Flags:\n");
  printf("makeMiniNtuple \t     Input file (-f)       : %s\n",input_file);
  printf("makeMiniNtuple \t     Number of events (-n) : %i\n",n_events );
  printf("makeMiniNtuple \t =================================================================\n");
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
    ChainHelper::addFileList(chain, input_file);
  } else if (inputIsDir) {
    ChainHelper::addFileDir(chain, input_file);
  } else {
    printf("makeMiniNtuple \t Cannot understand input %s",input_file);
    return 0;  
  }

  Superflow* cutflow = new Superflow(); // initialize the cutflow
  cutflow->setAnaName("SuperflowAna");                // arbitrary
  cutflow->setAnaType(AnalysisType::Ana_2Lep);        // analysis type, passed to SusyNt ?
  cutflow->setLumi(LUMI_A_A4);                        // set the MC normalized to lumi periods A1-A4
  cutflow->setSampleName(input_file);                 // sample name, check to make sure it's set OK
  cutflow->setRunMode(run_mode);                      // make configurable via run_mode
  cutflow->setCountWeights(true);                     // print the weighted cutflows
  cutflow->setChain(chain);

  printf("makeMiniNtuple \t Total events available : %lli\n",chain->GetEntries());

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

  *cutflow << CutName("bad muon veto") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_BadMuon);
  };
  
  *cutflow << CutName("jet cleaning") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_BadJet);
  };
  
  *cutflow << CutName("pass good vertex") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_GoodVtx);
  };
  
  *cutflow << CutName("pass cosmic veto") << [&](Superlink* /*sl*/) -> bool {
      return (cutflags & ECut_Cosmic);
  };

  //  Analysis Cuts
  *cutflow << CutName("at least 2 baseline leptons") << [](Superlink* sl) -> bool {
      return sl->baseLeptons->size() >= 2;
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

  // Lepton variables
  LeptonVector baseLeptons;
  *cutflow << [&](Superlink* sl, var_void*) { baseLeptons = *sl->baseLeptons; };

  *cutflow << NewVar("number of baseline leptons"); {
    *cutflow << HFTname("nBaseLeptons");
    *cutflow << [&](Superlink* /*sl*/, var_int*) -> int { return baseLeptons.size(); };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton flavor (0: e, 1: m)"); {
    *cutflow << HFTname("l_flav");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> lep_flav;
      for(unsigned int i = 0; i < baseLeptons.size(); i++) {
        lep_flav.push_back(baseLeptons.at(i)->isEle() ? 0 : 1);
      }
      return lep_flav;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton pt"); {
    *cutflow << HFTname("l_pt");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> lep_pt;
      for(unsigned int i = 0; i< baseLeptons.size(); i++) {
        lep_pt.push_back(baseLeptons.at(i)->Pt());
      }
      return lep_pt;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton eta"); {
    *cutflow << HFTname("l_eta");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> lep_eta;
      for(unsigned int i = 0; i < baseLeptons.size(); i++) {
        lep_eta.push_back(baseLeptons.at(i)->Eta());
      }
      return lep_eta;
    };
    *cutflow << SaveVar();
  }
  *cutflow << NewVar("lepton phi"); {
    *cutflow << HFTname("l_phi");
    *cutflow << [&](Superlink* /*sl*/, var_float_array*) -> vector<double> {
      vector<double> lep_phi;
      for(unsigned int i = 0; i < baseLeptons.size(); i++) {
        lep_phi.push_back(baseLeptons.at(i)->Phi());
      }
      return lep_phi;
    };
    *cutflow << SaveVar();
  }
  // Jet variables
  // MET variables

  ///////////////////////////////////////////////////////////////////////
  // Superflow methods end here
  ///////////////////////////////////////////////////////////////////////

  // Initialize the cutflow and start the event loop.
  chain->Process(cutflow, input_file, n_events, n_skip_events);

  delete cutflow;
  delete chain;

  // Print information
  printf("makeMiniNtuple \t =================================================================\n");
  printf("makeMiniNtuple \t All done!\n");
  printf("makeMiniNtuple \t =================================================================\n");
  return 0; 
}
