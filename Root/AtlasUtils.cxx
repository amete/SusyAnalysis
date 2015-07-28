#include <iostream>
#include <cmath>

#include "SusyAnalysis/AtlasUtils.h"

#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPave.h"

#include "TObject.h"
#include "TList.h"

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color,string additional) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  string text = "ATLAS"; 
  if( additional != "" ) text = text + " " + additional;
  l.DrawLatex(x,y,text.c_str());
}

TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2) {
 
  const Int_t debug=0; 

  if (!g1) printf("**myTGraphErrorsDivide: g1 does not exist !  \n"); 
  if (!g2) printf("**myTGraphErrorsDivide: g2 does not exist !  \n"); 


  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
   printf("**myTGraphErrorsDivide: vector do not have same number of entries !  \n"); 
  }

  TGraphErrors* g3= new TGraphErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0.;
  Double_t dx1=0.,dy1=0.,       dy2=0.;

  Int_t iv=0;
  for (Int_t i1=0; i1<n1; i1++) {
   for (Int_t i2=0; i2<n2; i2++) {
     //if (debug) printf("**myTGraphErrorsDivide: %d  %d !  \n",i1,i2);

    g1->GetPoint(i1,x1,y1);
    g2->GetPoint(i2,x2,y2);
    if (x1!=x2) {
      //printf("**myTGraphErrorsDivide: %d x1!=x2  %f %f  !  \n",iv,x1,x2);
    }else{
      //if (debug) printf("**myTGraphErrorsDivide: %d x1=x2  %f %f  !  \n",iv,x1,x2);
     dx1  = g1->GetErrorX(i1);
     if (y1!=0) dy1  = g1->GetErrorY(i1)/y1;
     if (y2!=0) dy2  = g2->GetErrorY(i2)/y2;
   
     if (debug)
      printf("**myTGraphErrorsDivide: %d x1=%f x2=%f y1=%f y2=%f  \n",iv,x1,x2,y1,y2);

     if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
     else        g3->SetPoint(iv, x1,y2);
   
     Double_t e=0.;
     if (y1!=0 && y2!=0) e=std::sqrt(dy1*dy1+dy2*dy2)*(y1/y2); 
     g3->SetPointError(iv,dx1,e);


     if (debug) {
       //Double_t g3y, g3x,g3e;
       //g3->GetPoint(iv, g3y,g3x);
       //g3e=g3->GetErrorY(iv);
       //printf("%d g3y= %f g3e=%f  \n",iv,g3y,g3e);
     }
     iv++;
    }
    //    printf("**myTGraphErrorsDivide: ...next  \n");
   }
  }  
  return g3;

}


TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

  const Int_t debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  Double_t   x1=0.,   y1=0., x2=0., y2=0.;
  Double_t dx1h=0., dx1l=0.;
  Double_t dy1h=0., dy1l=0.;
  Double_t dy2h=0., dy2l=0.;

  //Double_t* X1 = g1->GetX();
  //Double_t* Y1 = g1->GetY();
  Double_t* EXhigh1 = g1->GetEXhigh();
  Double_t* EXlow1 =  g1->GetEXlow();
  Double_t* EYhigh1 = g1->GetEYhigh();
  Double_t* EYlow1 =  g1->GetEYlow();

  //Double_t* X2 = g2->GetX();
  //Double_t* Y2 = g2->GetY();
  //Double_t* EXhigh2 = g2->GetEXhigh();
  //Double_t* EXlow2 =  g2->GetEXlow();
  Double_t* EYhigh2 = g2->GetEYhigh();
  Double_t* EYlow2 =  g2->GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
   
    //if (debug)
    //printf("%d x1=%f x2=%f y1=%f y2=%f  \n",i,x1,x2,y1,y2);
    if (debug)
      printf("%d dy1=%f %f dy2=%f %f sqrt= %f %f \n",i,dy1l,dy1h,dy2l,dy2h,
	     std::sqrt(dy1l*dy1l+dy2l*dy2l), std::sqrt(dy1h*dy1h+dy2h*dy2h));

    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    Double_t el=0.; Double_t eh=0.;

    if (y1!=0. && y2!=0.) el=std::sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
    if (y1!=0. && y2!=0.) eh=std::sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);

    if (debug) printf("dx1h=%f  dx1l=%f  el=%f  eh=%f \n",dx1h,dx1l,el,eh);
    g3->SetPointError(i,dx1h,dx1l,el,eh);

  }  
  return g3;

}


void myTGraphErrorsAdd(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2){

  if( g1->GetN() != g2->GetN() ){
    cout<<"TGraphErrorsAdd: NOT same number of points!"<<endl;
    return;
  }

  for(int i=0; i<g1->GetN(); i++){

    float y1   = g1->GetY()[i]; 
    float ehy1 = g1->GetEYhigh()[i]; 
    float ely1 = g1->GetEYlow()[i]; 

    float y2   = g2->GetY()[i]; 
    float ehy2 = g2->GetEYhigh()[i]; 
    float ely2 = g2->GetEYlow()[i]; 

    g1->GetY()[i] = y1+y2; 
    g1->GetEYhigh()[i] = sqrt(ehy1*ehy1 + ehy2*ehy2); 
    g1->GetEYlow()[i]  = sqrt(ely1*ehy1 + ely2*ely2); 
    
  }//End looping through points

}//End add function


TGraphAsymmErrors* myMakeBand(TGraphAsymmErrors* g0, 
			      TGraphAsymmErrors* g1,
			      TGraphAsymmErrors* g2) {
  // default is g0
    //const Int_t debug=0;

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;
  //Double_t dx1=0.;
  Double_t dum;
  for (Int_t i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    // if (y1==0) y1=1;
    //if (y2==0) y2=1;

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    Double_t tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    //Double_t y3=1.;
    Double_t y3=y0;
    g3->SetPoint(i,x1,y3);

    Double_t binwl=(x1-x3)/2.;
    Double_t binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,binwl,binwh,(y3-y2),(y1-y3));

  }
  return g3;

}

void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2) {

  Double_t  x1=0., y1=0.,  y2=0., y0=0;
  //Double_t dx1=0.;
  //Double_t dum;

  if (g1->GetN()!=g2->GetN())
    std::cout << " graphs have not the same # of elements " << std::endl;

  Double_t* EYhigh = g2-> GetEYhigh();
  Double_t* EYlow  = g2-> GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    if (y1==0) y1=1;
    if (y2==0) y2=1;

    //    if (i==g1->GetN()-1) x2=x1;
    //    else                 g2->GetPoint(i+1,x2,dum);
    //    if (i==0)            x3=x1;
    //    else                 g2->GetPoint(i-1,x3,dum);

    Double_t eyh=0., eyl=0.;
    //if (y1<y2) {y2=y1; y1=tmp;}
    //Double_t y3=1.;

    //printf("%d: y1=%f y2=%f Eyhigh= %f Eylow= %f \n",i,y1,y2,EYhigh[i],EYlow[i]);

    y0=y1-y2;
    if (y0!=0) {
     if (y0>0){
      eyh=EYhigh[i];
      eyh=std::sqrt(eyh*eyh+y0*y0);
      //printf("high: %d: y0=%f eyh=%f  \n",i,y0,eyh);
      g2->SetPointEYhigh(i,eyh);
     } else {
      eyl=EYlow[i];
      eyl=std::sqrt(eyl*eyl+y0*y0);
      // printf("low: %d: y0=%f eyl=%f  \n",i,y0,eyl);
      g2->SetPointEYlow (i,eyl);
     }
    }
  }
  return;

}

TGraphAsymmErrors* TH1TOTGraph(TH1 *h1){


  if (!h1) std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

 TGraphAsymmErrors* g1= new TGraphAsymmErrors();

 Double_t x, y, ex, ey;
/*
 for (Int_t i=1; i<=h1->GetNbinsX()+1; i++) {
   y=h1->GetBinContent(i-1);
  ey=h1->GetBinError(i-1);
   x=h1->GetBinCenter(i-1);
  ex=h1->GetBinWidth(i-1)/2.;
  //   cout << " x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;
   g1->SetPoint(i,x,y);
   g1->SetPointError(i,ex,ex,ey,ey);
 }
*/

 // Don't care about the underflow/overflow
 for (Int_t i=1; i<=h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i);
  ey=h1->GetBinError(i);
   x=h1->GetBinCenter(i);
  ex=h1->GetBinWidth(i)/2.;
   g1->SetPoint(i-1,x,y);
   g1->SetPointError(i-1,ex,ex,ey,ey);
 }

 //g1->Print();

 return g1;
}

void TGraphToTH2(TH2* h2, TGraph* g){

  if( h2 == 0x0 ){
    std::cout<<"h2 == NULL!"<<std::endl;
    return;
  }
  if( g  == 0x0){
    std::cout<<"g == NULL!"<<std::endl;
    return;
  }

  double x = 0; 
  double y = 0;

  for(int i=0; i<g->GetN(); i++){ //Loop over points
    g->GetPoint(i, x, y); 
    h2->Fill(x, y); 
  }//End looping over points


}//End TGraphToTH2

/*
void myText(Double_t x,Double_t y,Color_t color,char *text) {

  //Double_t tsize=0.05;
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
 

void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,char *text) 
{

  Double_t tsize=0.06;

  TLatex l; l.SetTextAlign(12); //l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");

  mbox->SetFillColor(mcolor);
  mbox->SetFillStyle(1001);
  mbox->Draw();

  TLine mline;
  mline.SetLineWidth(4);
  mline.SetLineColor(1);
  mline.SetLineStyle(1);
  Double_t y_new=(y1+y2)/2.;
  mline.DrawLineNDC(x1,y_new,x2,y_new);

}
*/

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,char *text,Float_t msize) 
{
  Double_t tsize=0.06;
  TMarker *marker = new TMarker(x-(0.4*tsize),y,8);
  marker->SetMarkerColor(color);  marker->SetNDC();
  marker->SetMarkerStyle(mstyle);
  marker->SetMarkerSize(msize);
  marker->Draw();

  TLatex l; l.SetTextAlign(12); //l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

// ADDED BY BKT --------------------------------------------------------------

TCanvas* sqAtlasCanvas(const char* name, bool ratio){
  TCanvas* atlas_square = new TCanvas(name,name,0.,0.,600,600);
  atlas_square->ToggleEventStatus();
  if( !ratio )   return atlas_square; 
  
  string pad1name = name; 
  pad1name = pad1name + " p";
  TPad *p1 = new TPad(pad1name.c_str(),pad1name.c_str(),0.05,0.40,0.95,1.0);
  p1->SetBottomMargin(0.00);
  p1->SetNumber(1);

  string pad2name = name; 
  pad2name = pad1name + " r";
  TPad *p2 = new TPad(pad2name.c_str(),pad2name.c_str(),0.05,0.1,0.95,0.4);
  p2->SetNumber(2);
  p2->SetTopMargin(0.00); 

  atlas_square-> cd();
  p1->Draw();
  p1->cd(); 

  atlas_square-> cd();
  p2->Draw();
  p2->cd();

  atlas_square->cd(0);   

  return atlas_square; 

}

TCanvas* recAtlasCanvas(const char* name, bool ratio){
  TCanvas* atlas_rectangular = new TCanvas(name,name,0.,0.,800,600);
  atlas_rectangular->ToggleEventStatus();
  if( !ratio )   return atlas_rectangular; 
  
  string pad1name = name; 
  pad1name = pad1name + " p";
  TPad *p1 = new TPad(pad1name.c_str(),pad1name.c_str(),0.05,0.40,0.95,1.0);
  p1->SetBottomMargin(0.00);
  p1->SetNumber(1);

  string pad2name = name; 
  pad2name = pad1name + " r";
  TPad *p2 = new TPad(pad2name.c_str(),pad2name.c_str(),0.05,0.1,0.95,0.4);
  p2->SetNumber(2);
  p2->SetTopMargin(0.00); 
  
  atlas_rectangular->cd();
  p1->Draw();
  p1->cd(); 

  atlas_rectangular->cd();
  p2->Draw();
  p2->cd();

  atlas_rectangular->cd(0);   

  return atlas_rectangular; 
}

TCanvas* colzAtlasCanvas(const char* name){

  TCanvas* c_colz = recAtlasCanvas(name);
  //c_colz->SetLeftMargin(0.25); 
  c_colz->SetTopMargin(0.1);
  c_colz->SetRightMargin(0.25);

  return c_colz; 

}

void ratioCanvas(TCanvas* c){

   TVirtualPad* _tv = c->cd();
   TPad* _pTop = new TPad("pTop","pTop",0,0.3,1,1);
   _pTop->SetTopMargin(0.075);
   _pTop->SetBottomMargin(0.02);
   _pTop->SetRightMargin(0.05);
   _pTop->SetLeftMargin(0.15);
   _pTop->SetNumber(1);
   TPad* _pBot = new TPad("pBot","pBot",0,0,1,0.3);
   _pBot->SetTopMargin(0.04);
   _pBot->SetBottomMargin(0.4);
   _pBot->SetRightMargin(0.044);
   _pBot->SetLeftMargin(0.15);
   _pBot->SetNumber(2);
   _tv->cd();
 
   _pBot->Draw(); 
   _pTop->Draw();
   _pTop->cd();

}

void myPrint(TCanvas* c, string name){
  
  string eps = name + ".eps";
  string png = name + ".png";
  string pdf = name + ".pdf"; 

  c->Print(eps.c_str());
  c->Print(png.c_str());
  c->Print(pdf.c_str());
  
  return;
}

void SetAltDots(TH1* h, int color){
  
  h->SetMarkerStyle(24);
  h->SetMarkerColor(color);
  h->SetLineColor(color);
  h->Paint(); 
  TCanvas cDummy; 
  TPaveStats* p = (TPaveStats*) h->FindObject("stats");
  if( p ){
    p->SetLineColor(color);
    p->SetTextColor(color);
  }

  return; 

}

void SetUnderFlow(TH1F* h){
   
  int nX = h->GetXaxis()->GetFirst(); 

  float underflow = 0.; 
  float underflowerr = 0.; 
  for(int i=0; i < nX; i++){
    underflow = underflow+h->GetBinContent(i); 
    underflowerr = underflowerr+h->GetBinError(i)*h->GetBinError(i); 
    h->SetBinContent(i, 0.); 
    h->SetBinError(i, 0.); 
  }

  h->SetBinContent( nX, underflow );
  h->SetBinError(nX, sqrt(underflowerr)); 
   
  return; 
}

void SetOverFlow(TH1F* h){
   
  int nX = h->GetXaxis()->GetLast(); 

  float overflow = 0.; 
  float overflowerr = 0.; 
  for(int i=nX+1; i<=h->GetNbinsX()+1; i++){
    overflow = overflow+h->GetBinContent(i); 
    overflowerr = overflowerr+h->GetBinError(i)*h->GetBinError(i); 
    h->SetBinContent(i, 0.); 
    h->SetBinError(i, 0.); 
  }

  h->SetBinContent( nX, overflow );
  h->SetBinError(nX, sqrt(overflowerr)); 
   
  return; 
}

void SetEdges(TH1F* h){
  SetUnderFlow(h); 
  SetOverFlow(h);
  return;
}

void CleanLegend(TLegend* leg){

  leg->SetFillColor(kWhite); 
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);

  return;
}

void TwoD2Latex(TH2* h, string filename, int /*prec*/){

   ofstream out;
 
   filename = filename + ".tex";
   out.open(filename.c_str()); 
 
   //Print header info
   out << "%% This was made with TwoD2Latex " << endl;
   out << "\\documentclass[12pt]{article}  " << endl 
       << "\\usepackage{geometry} " << endl 
       << "\\usepackage{graphicx} " << endl 
       << "\\usepackage{latexsym} " << endl 
       << "\\usepackage{atlasphysics} "<<endl
       << "\\usepackage{pdflscape} "<<endl
       << endl
       << "\\geometry{a4paper} " << endl
       << endl
       << "\\title{" << h->GetTitle() << "} " << endl   
       << "\\author{Brokk Toggerson} " << endl
       << "%%"<<"\\date{} % delete this line to display the current date " 
       << endl
       << "%%"<<"%%"<<"%%"<<" BEGIN DOCUMENT " << endl
       << "\\begin{document} " << endl
       << "\\begin{landscape} "<< endl 
       << "\\begin{table} " << endl
       << "\\begin{center} " << endl
       << "\\begin{tabular}{|l|" ;
 
   //Set the properties for the columns
   for(int c=1; c<=h->GetXaxis()->GetNbins(); c++) out << "r|";
   out<<"} "<<"\\"<<"hline " <<endl;
 
   //Column title
   out<<h->GetTitle()<<" &"<<"\\"<<"multicolumn{"<<h->GetXaxis()->GetNbins()<<"}{c|}"
      <<"{"<<h->GetXaxis()->GetTitle()<<"} "<<"\\"<<"\\"<<" \\hline"<<endl;
   //Column headings
   out<<h->GetYaxis()->GetTitle()<<" ";
   for(int c=1; c<=h->GetXaxis()->GetNbins(); c++){
     string columnLabel = h->GetXaxis()->GetBinLabel(c);  
     if( !columnLabel.empty() ) out << "& " << h->GetXaxis()->GetBinLabel(c);
     else out << "&" << h->GetXaxis()->GetBinCenter(c); 

   }
   out<<" \\"<<"\\"<<" \\"<<"hline " <<endl;
 
   //Fill the rows
/*
   for(int r=1; r<=h->GetYaxis()->GetNbins(); r++){
     string rowLabel = h->GetYaxis()->GetBinLabel(c);  
     if( !rowLabel.empty() ) out<<h->GetYaxis()->GetBinLabel(r)<<" ";
     //else out << h->GetYaxis()->GetBinCenter(r); 
     for(int c=1; c<=h->GetXaxis()->GetNbins(); c++) {
       float val = h->GetBinContent(c, r);
       float err = h->GetBinError(c, r); 
       out<<"& "<<setprecision(prec)<<val<<" $ \\pm $ "<<err;
     }//End looping over columns 
     out<<"\\"<<"\\"<<endl;
   }//End looping over rows
   out<<"\\"<<"hline "<<endl;
*/
   //The tail stuff
   out << "\\end{tabular} " << endl
       << "\\end{center} " << endl
       << "\\end{table} " << endl
       << "\\end{landscape}" <<endl
       << endl
       << "\\end{document} "<<endl;
 
   out.close();
 
   return; 
 
}//End TwoD2Latex()


void myDraw1D(TCanvas* c, TH1* hist, string option){

  float yMax=hist->GetMaximum();

  //find if the pad has already a histogram and rescale the Y-axis accordingly.
  TList* dummy = 0x0;
  dummy = c->GetListOfPrimitives(); 
  TIter next(dummy); 
  TObject *obj; 
  TH1* h1=NULL; //The first histogram drawn
  int nHisto=0;
  while ((obj=next())) { 
    if (obj->InheritsFrom("TH1")) { 
      TH1 *h = (TH1*)obj;
      nHisto++;
      if(h1==NULL) h1=h;
      if(h->GetMaximum()>yMax)	yMax=h->GetMaximum();
    }//If obj is TH1 
  }//End looping through list of primitives

  if(h1){
    float sf = 1.2; 
    if( c->GetLogy() > 0 ) sf = 10.; 
    h1->SetMaximum(yMax*sf);
  }
  
  hist->Draw(option.c_str());

}//End myDraw1D


void PrintTLorentzVector(TLorentzVector& lv){
  cout<<"(E, px, py, pz)  = ("<<lv.E()<<", "<<lv.Px()<<", "<<lv.Py()<<", "<<lv.Pz()<<")"<<endl;
}


void fixRatioAxes(TAxis* xt, TAxis* /*yt*/, TAxis* xb, TAxis* yb){
  
  xt->SetLabelSize(0.);
  xt->SetTitleSize(0.);

  xb->SetLabelSize(0.12);
  yb->SetLabelSize(0.12);

  xb->SetTitleSize(0.15);
  xb->SetTitleOffset(1.2);

  yb->SetTitleSize(0.11);
  yb->SetTitleOffset(0.65);
  yb->SetNdivisions(507);

}

void cdfDataErr(TGraphAsymmErrors*& gout){

  /** From http://www-cdf.fnal.gov/physics/statistics/notes/pois_eb.txt */

  for(int n=0; n < gout->GetN(); n++){

    float val = gout->GetY()[n]; 
    float errP = 0.;
    float errN = 0.;
    if( val > 0.1 ){
      errP =  0.5 + sqrt(val+0.25); 
      errN = -0.5 + sqrt(val+0.25);
    }

    gout->SetPointError(n, 0., 0., errN, errP); 
    
  }//End looping through points

}//End cdfDataErr

void skyDataErr(TGraphAsymmErrors*& gout){

  for(int n=0; n < gout->GetN(); n++){

    float y = gout->GetY()[n]; 

    float errP = 0.;
    if( y > 0 ){
      float y1 = y + 1; 
      float d = 1. - 1./(9.*y1) + 1./(3.*sqrt(y1)); 
      errP = y1*d*d*d - y; 
    }
    float errN = 0; 
    if( y > 0 ){
      float d = 1. - 1./(9.*y) - 1./(3.*sqrt(y)); 
      errN = y - y*d*d*d;
    }

    gout->SetPointError(n, 0., 0., errN, errP); 

  }//End looping through points 

}

void removeZero(TGraphAsymmErrors*& gout){

  for(int n=gout->GetN()-1; n>=0; n--)
    if( fabs( gout->GetY()[n] ) < 0.00001 ) gout->RemovePoint(n);
    

}

TLine makeLOne(TH1F* h_ratio){

  TLine ans(h_ratio->GetXaxis()->GetXmin(), 1., 
	    h_ratio->GetXaxis()->GetXmax(), 1.);
  ans.SetLineWidth(1); 
  ans.SetLineStyle(2); 
  ans.SetLineColor(kBlack); 
  return ans;

}
