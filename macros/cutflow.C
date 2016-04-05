// \brief : Quick script to check cutflow for stop 2L
//void cutflow(TString filename = "/scratch/amete/Stop2LTest/CENTRAL_410000.root")
//void cutflow(TString filename = "../../../cutflow/CENTRAL_361068.root")
void cutflow(TString filename = "../../../cutflow/CENTRAL_406009.root")
{
  TCut channels[2]     = { "l_flav[0]==l_flav[1]", "(l_flav[0]!=l_flav[1])" }; // SF - DF
  TCut signs[2]        = { "(l_q[0]*l_q[1])<0"   , "(l_q[0]*l_q[1])>0" };      // OS - SS
  TCut selection[7]    = { "mll>20.", 
                           "mll>20&&met>40", 
                           "mll>20&&met>40&&R1>0.2", 
                           "mll>20&&met>40&&R1>0.2&&MDR_jigsaw>95.",
                           "mll>20&&met>40&&R1>0.2&&MDR_jigsaw>95.&&RPT_jigsaw>0.5",
                           "mll>20&&met>40&&R1>0.2&&MDR_jigsaw>95.&&RPT_jigsaw>0.5&&gamInvRp1_jigsaw>0.8",
                           "mll>20&&met>40&&R1>0.2&&MDR_jigsaw>95.&&RPT_jigsaw>0.5&&gamInvRp1_jigsaw>0.8&&DPB_vSS_jigsaw>(0.8*TMath::Abs(cthllb)+1.8)",
                           //"mll>20&&met>40&&R1>0.2&&mT2lep>20.",
                           //"mll>20&&met>40&&R1>0.2&&mT2lep>20.&&TMath::Abs(deltaX)<0.05",
                           //"mll>20&&met>40&&R1>0.2&&mT2lep>20.&&TMath::Abs(deltaX)<0.05&&R2>0.5",
                           //"mll>20&&met>40&&R1>0.2&&mT2lep>20.&&TMath::Abs(deltaX)<0.05&&R2>0.5&&TMath::Abs(cthllb)<0.8"
  };
  TString selHRName[7] = { "mll>20 GeV\t\t", 
                           "met>40 GeV\t\t", 
                           "R1>0.2\t\t\t", 
                           "MDR>95. GeV\t\t",
                           "RPT>0.5\t\t\t",
                           "gamInvRp1>0.8\t\t",
                           "DPB_vSS>(0.8*|cthllb|+1.8)"
                           //"mT2lep>20.",
                           //"|deltaX|<0.05",
                           //"R2>0.5\t",
                           //"|cthllb|<0.8"
  };

  TFile* file = new TFile(filename);
  if(!file->IsOpen()) return;
  TTree* tree = (TTree*) file->Get("superNt");
  if(!tree) return;
  double yields[2][2][7] = {{0.}};
  TH1F* histo = new TH1F("histo","histo",1,0.5,1.5); histo->Sumw2();

  // Read
  for(unsigned int i=0; i<2; ++i) {
    for(unsigned int j=0; j<2; ++j) {
      for(unsigned int k=0; k<7; ++k) {
        tree->Draw("1>>histo",channels[i]+signs[j]+selection[k],"goff");
        yields[i][j][k]=histo->Integral();
      } // selection
    } // signs
  } // channels

  // Print
  std::cout << "\n\n\t\t\t\tSFOS\t\tSFSS\t\tDFOS\t\tDFSS" << std::endl;
  for(unsigned int i=0; i<7; ++i) {
    std::cout << selHRName[i] << "\t";
    for(unsigned int j=0; j<2; ++j) {
      for(unsigned int k=0; k<2; ++k) {
        std::cout << yields[j][k][i] << "\t\t";
      } // signs
    } // channels
    std::cout << std::endl;
  } // selection
  std::cout << "\n\n" << std::endl;
  
  //tree->SetScanField(0);
  //tree->Scan("eventNumber:l_pt[0]:l_pt[1]:j_pt[0]:j_pt[1]:j_pt[2]:j_pt[3]:j_pt[4]:j_pt[5]:j_pt[6]:met:R1","eventNumber==7713609");
  //tree->Scan("eventNumber:mll:met:R1:mT2lep:deltaX:(TMath::Abs(deltaX)<0.05)",channels[1]+signs[1]+selection[3]);
}
