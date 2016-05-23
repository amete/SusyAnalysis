void make2DPlots(TString proc = "ttbar_CENTRAL", TString var1 = "R1", TString var2 = "R2" , TString cut = "l_pt[0]>20.&&l_pt[1]>20.&&(l_q[0]*l_q[1])<0&&l_flav[0]!=l_flav[1]&&mT2lep>100")
{
  TFile* file = new TFile("../mc15_13TeV_sherpa.root");
  if(!file->IsOpen()) return;
  TTree* tree = (TTree*) file->Get(proc);
  if(tree==nullptr) return;
  TH2D* hist = new TH2D("hist","hist",20,0,1,20,0,1);
  hist->Sumw2();
  hist->GetXaxis()->SetTitle(var1);
  hist->GetYaxis()->SetTitle(var2);
  hist->SetTitle("");
  TString drawthis = var2 + ":" + var1 + ">>hist"; //drawthis.Form("%s:%s>>hist",var1,var2);
  tree->Draw(drawthis,cut,"goff");
  hist->Scale(1/hist->Integral());
  TCanvas* can = new TCanvas("can","can",500,500);
  can->SetFillColor(0);
  can->cd();
  gStyle->SetOptStat(0);
  hist->Draw("COLZ");
  can->SaveAs(proc+"_"+var1+"_vs_"+var2+"_OSDF.eps");
  return;
}
