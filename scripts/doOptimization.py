#!/usr/bin/python
import ROOT,sys,math

LUMI    = 35000.
VERBOSE = False
CHANNEL = "df"
BINNED  = True 
SIGNAL  = "LOWDM" # "LOWDM" "HIGHDM"

# Define the cuts
def trigger():
    return "(((pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14||pass_HLT_mu18_mu8noL1)&&treatAsYear==2015)||((pass_HLT_2e17_lhvloose_nod0||pass_HLT_e17_lhloose_nod0_mu14||pass_HLT_mu22_mu8noL1)&&treatAsYear==2016))"

def osptcuts(cut1,cut2):
    return "(l_q[0]*l_q[1])<0&&l_pt[0]>%s.&&l_pt[1]>%s."%(cut1,cut2)

def issf():
    return "l_flav[0]==l_flav[1]"

def isdf():
    return "l_flav[0]!=l_flav[1]"

def mllcut(cut1,cut2=1.e10):
    return "mll>%s.&&mll<%s"%(cut1,cut2)

def mt2cut(cut1,cut2=1.e10):
    return "mT2lep>%s.&&mT2lep<%s"%(cut1,cut2)

def zveto(cut):
    return "TMath::Abs(mll-91.2)>%s."%(cut)

def zselect(cut):
    return "TMath::Abs(mll-91.2)<%s"%(cut)

def cljveto(cut):
    if "None" in cut: return "1"
    if "20" in cut: return "nCentralLJets==0"
    return "nCentralLJets%s==0"%(cut)

def cbjveto(cut):
    if "None" in cut: return "1"
    if "20" in cut: return "nCentralBJets==0"
    return "nCentralBJets==0"

def fjveto(cut):
    if "None" in cut: return "1"
    if "20" in cut: return "nForwardJets==0"
    return "nForwardJets%s==0"%(cut)

# Calculate ZBi
def calculateZBi(selection,lumi,signal,backgroundList):
    if VERBOSE:
        print "="*35
    # Get Total BG
    total_bg = 0.
    for bkg in backgroundList:
        htemp = ROOT.TH1F('htemp','htemp',1,0.5,1.5)
        htemp.Sumw2()
        backgroundList[bkg].Draw('1>>htemp','%s*eventweight*(%s)'%(lumi,selection),'goff') 
        error    = ROOT.Double(0.)
        integral = htemp.IntegralAndError(0,-1,error);
        total_bg = total_bg + integral 
        if VERBOSE:
            print '%*s = %*.2f +/- %*.2f'%(10,bkg,5,integral,5,error)       
        htemp.Clear()
    if VERBOSE:
        print "="*35
        print '%*s = %.2f'%(10,"Total BG",total_bg)       
        print "="*35
    # Get the signal and calculate
    htemp = ROOT.TH1F('htemp','htemp',1,0.5,1.5)
    htemp.Sumw2()
    signal.Draw('1>>htemp','%s*eventweight*(%s)'%(lumi,selection),'goff') 
    error    = ROOT.Double(0.)
    total_sig = htemp.IntegralAndError(0,-1,error);
    significance = ROOT.RooStats.NumberCountingUtils.BinomialExpZ(total_sig,total_bg,0.3)
    if VERBOSE:
        print '%*.2f +/- %*.2f with ZBi = %.2f'%(5,total_sig,5,error,significance)       
    htemp.Clear()
    if VERBOSE:
        print "="*35
    return significance,total_sig,total_bg

# Main function
def main():
    ROOT.RooStats.NumberCountingUtils

    # Open the file
    inputFile = ROOT.TFile('/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/hfts/HFT_COMBINED_13TeV.root','READ')

    # Open the output file
    scantype = "unbinned"
    if BINNED : scantype = "binned"
    outputFile = open('optimization_result_%s_%s_%s.txt'%(scantype,CHANNEL,SIGNAL),'w')

    # Get Background Trees
    backgroundList = {'ttbar' : 0. ,'singletop' : 0. , 'ttv' : 0. ,'VV' : 0. ,'VVV' : 0. ,'higgs' : 0. ,'W' : 0. ,'Z' : 0. }
    for bkg in backgroundList:
       backgroundList[bkg] = inputFile.Get('%s_CENTRAL'%(bkg)) 
       if VERBOSE:
           print backgroundList[bkg]

    # Get Signal Trees
    if SIGNAL == "LOWDM":
        signalList     = {'c1c1_slep_300.0_100.0' : 0.}
    elif SIGNAL == "HIGHDM":
        signalList     = {'c1c1_slep_700.0_1.0'   : 0.}
    else:
        signalList     = {'c1c1_slep_150.0_50.0'  : 0.,
                          'c1c1_slep_200.0_100.0' : 0.,
                          'c1c1_slep_300.0_100.0' : 0.,
                          'c1c1_slep_300.0_200.0' : 0.,
                          'c1c1_slep_400.0_100.0' : 0.,
                          'c1c1_slep_400.0_200.0' : 0.,
                          'c1c1_slep_400.0_300.0' : 0.,
                          'c1c1_slep_500.0_1.0'   : 0.,
                          'c1c1_slep_500.0_100.0' : 0.,
                          'c1c1_slep_500.0_200.0' : 0.,
                          'c1c1_slep_500.0_300.0' : 0.,
                          'c1c1_slep_600.0_1.0'   : 0.,
                          'c1c1_slep_600.0_100.0' : 0.,
                          'c1c1_slep_600.0_200.0' : 0.,
                          'c1c1_slep_700.0_1.0'   : 0.,
                          'c1c1_slep_700.0_100.0' : 0.,
                          'c1c1_slep_700.0_300.0' : 0.}

    for sig in signalList:
       signalList[sig] = inputFile.Get('%s_CENTRAL'%(sig)) 
       if VERBOSE:
           print signalList[sig]

    # Baseline Selection
    if "sf" in CHANNEL:
        baseSelection = trigger()           + "&&" +\
                        issf()              + "&&" +\
                        osptcuts("25","20") + "&&" +\
                        mllcut("40")
    elif "df" in CHANNEL:
        baseSelection = trigger()           + "&&" +\
                        isdf()              + "&&" +\
                        osptcuts("25","20") + "&&" +\
                        mllcut("40")
    else:
        print "Unknown channel %s, quitting..."%(CHANNEL)
        return

    # zcut   : |Mll-91.2| > X GeV
    # mllcut : X_1 < Mll < X_2 (BINNED-only)
    if "sf" in CHANNEL:
        zcuts   = ["10"]
        mllcuts = ["100","150","200","300"]
    else:
        zcuts   = ["-1"]
        mllcuts = ["-1"]
    cljcuts = ["None","20","30","40","50","60"]
    cbjcuts = ["20"]
    fjcuts  = ["None"]
    mt2cuts = ["100","150","200","300"]

    # Single inclusive bin cut-and-count
    if not BINNED:
        outputFile.write( "="*150+"\n" )
        outputFile.write( "%*s \t %*s \t %*s \t %*s \t %*s \t %*s \t %*s \t %*s\n"%( 10,"Z-veto [GeV]",
                                                                                     10,"CLJ-veto [GeV]",
                                                                                     10,"CBJ-veto [GeV]",
                                                                                     10,"FJ-veto [GeV]",
                                                                                     10,"mT2 [GeV]",
                                                                                     10,"Total SIG", 
                                                                                     10,"Total BG", 
                                                                                     10,"ZBi") )
        outputFile.write( "="*150+"\n" )
        for sig in signalList:
            outputFile.write( "~"*150+"\n" )
            outputFile.write( "%*s Analyzing Signal Point %s \n"%(50,"",sig) )
            outputFile.write( "~"*150+"\n" )
            for zcut in zcuts:
                for cljcut in cljcuts:
                    for cbjcut in cbjcuts:
                        for fjcut in fjcuts:
                            for mtcut in mt2cuts:
                                selection = baseSelection   + "&&" +\
                                            zveto(zcut)     + "&&" +\
                                            cljveto(cljcut) + "&&" +\
                                            cbjveto(cbjcut) + "&&" +\
                                            fjveto(fjcut)   + "&&" +\
                                            mt2cut(mtcut)
                                sign,syield,bgyield = calculateZBi(selection,LUMI,signalList[sig],backgroundList)  
                                outputFile.write( "%*s \t %*s \t\t %*s \t\t %*s \t %*s \t %*.2f \t %*.2f \t %*.2f\n"%(10,zcut,
                                                                                                                      10,cljcut,
                                                                                                                      10,cbjcut,
                                                                                                                      10,fjcut,
                                                                                                                      10,mtcut,
                                                                                                                      10,syield,
                                                                                                                      10,bgyield,
                                                                                                                      10,sign) )
                                outputFile.write( "%*s -  %*s    - COMBINED ZBi = %*.2f\n"%(60,sig,5,CHANNEL,10,sign) )
        outputFile.write( "="*150+"\n" )
    # Orthogonal mll-mT2 bins combined in quadrature
    else:
        outputFile.write( "="*150+"\n" )
        outputFile.write( "%*s \t %*s \t %*s \t %*s \t %*s \t %*s \t %*s \t %*s\n"%( 10,"CLJ-veto [GeV]",
                                                                                     10,"CBJ-veto [GeV]",
                                                                                     10,"FJ-veto [GeV]",
                                                                                     10,"mll-window [GeV]",
                                                                                     10,"mT2-window [GeV]",
                                                                                     10,"Total SIG", 
                                                                                     10,"Total BG", 
                                                                                     10,"ZBi") )
        outputFile.write( "="*150+"\n" )
        for sig in signalList:
            outputFile.write( "~"*150+"\n" )
            outputFile.write( "%*s Analyzing Signal Point %s \n"%(50,"",sig) )
            outputFile.write( "~"*150+"\n" )
            for cljcut in cljcuts:
                for cbjcut in cbjcuts:
                    for fjcut in fjcuts:
                        outputFile.write( "\n"+"="*150+"\n" )
                        combinedSign = 0.
                        for jj,mtcut in enumerate(mt2cuts):
                            lowercut_mt2 = "%s"%(mt2cuts[jj])
                            uppercut_mt2 = "%s"%("1.e10")
                            if jj < len(mt2cuts)-1: uppercut_mt2 = "%s"%(mt2cuts[jj+1])
                            skipMll = False
                            for ii,mcut in enumerate(mllcuts):
                                lowercut_mll = "%s"%(mllcuts[ii])
                                uppercut_mll = "%s"%("1.e10")
                                if skipMll == True: break
                                if lowercut_mt2 == "300": 
                                    skipMll = True
                                    if CHANNEL == "sf" : lowercut_mll = "150"
                                elif ii < len(mllcuts)-1: uppercut_mll = "%s"%(mllcuts[ii+1])
                                selection = baseSelection                     + "&&" +\
                                            cljveto(cljcut)                   + "&&" +\
                                            cbjveto(cbjcut)                   + "&&" +\
                                            fjveto(fjcut)                     + "&&" +\
                                            mllcut(lowercut_mll,uppercut_mll) + "&&" +\
                                            mt2cut(lowercut_mt2,uppercut_mt2)
                                sign,syield,bgyield = calculateZBi(selection,LUMI,signalList[sig],backgroundList)  
                                outputFile.write( "%*s \t %*s \t\t %*s \t\t %*s%*s-%*s \t %*s%*s-%*s \t %*.2f \t %*.2f \t %*.2f\n"%(10,cljcut,
                                                                                                                                    10,cbjcut,
                                                                                                                                    10,fjcut,
                                                                                                                                    10,lowercut_mll,
                                                                                                                                    3,"", 7,uppercut_mll,
                                                                                                                                    10,lowercut_mt2,
                                                                                                                                    3,"", 7,uppercut_mt2,
                                                                                                                                    10,syield,
                                                                                                                                    10,bgyield,
                                                                                                                                    10,sign) )
                                if sign > 0: combinedSign = combinedSign + sign*sign
                        outputFile.write( "="*150+"\n" )
                        outputFile.write( "%*s -  %*s    - COMBINED ZBi = %*.2f"%(60,sig,5,CHANNEL,10,math.sqrt(combinedSign)) )
                        outputFile.write( "\n"+"="*150+"\n\n" )
        outputFile.write( "="*150+"\n" )

    # Once done close the files
    inputFile.Close()      
    outputFile.close()
 
if '__main__' in __name__:
    main()
