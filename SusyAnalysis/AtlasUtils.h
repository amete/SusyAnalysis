//
//   @file    AtlasUtils.h         
//   
//
//   @author M.Sutton
// 
//   Copyright (C) 2010 Atlas Collaboration
//
//   $Id: AtlasUtils.h 90845 2013-05-20 22:32:06Z amete $


#ifndef __ATLASUTILS_H
#define __ATLASUTILS_H

#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLorentzVector.h"
#include "TH2.h"

#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

void ATLAS_LABEL(Double_t x,Double_t y,Color_t color=1,string additional=""); 

TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2);

TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

void myTGraphErrorsAdd(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2); 

TGraphAsymmErrors* myMakeBand(TGraphAsymmErrors* g0, 
			      TGraphAsymmErrors* g1,
			      TGraphAsymmErrors* g2);

void myAddtoBand(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2);

TGraphAsymmErrors* TH1TOTGraph(TH1 *h1);

void TGraphToTH2(TH2* h2, TGraph* g); 

void myText(Double_t x,Double_t y,Color_t color,char *text);

void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,char *text);

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle,char *text,Float_t msize=2.); 

// ADDED by BKT ---------------------------------------------------------------

TCanvas* sqAtlasCanvas(const char* name, bool ratio=false);

TCanvas* recAtlasCanvas(const char* name, bool ratio=false); 

TCanvas* colzAtlasCanvas(const char* name);

void ratioCanvas(TCanvas* c); 

void myPrint(TCanvas* c, string name); 

void SetAltDots(TH1* h, int color); 

void SetUnderFlow(TH1F* h); 

void SetOverFlow(TH1F* h); 

void SetEdges(TH1F* h); 

void CleanLegend(TLegend* leg); 

void TwoD2Latex(TH2* h, string filename, int prec = 3); 

void myDraw1D(TCanvas* c, TH1* hist, string option = ""); 

void PrintTLorentzVector(TLorentzVector& lv); 

void fixRatioAxes(TAxis* xt, TAxis* yt, TAxis* xb, TAxis* yb);

void cdfDataErr(TGraphAsymmErrors*& gout);

void skyDataErr(TGraphAsymmErrors*& gout); 

void removeZero(TGraphAsymmErrors*& gout); 

TLine makeLOne(TH1F* h_ratio); 

#endif // __ATLASUTILS_H
