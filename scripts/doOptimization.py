#!/usr/bin/python
import ROOT,sys,math

VERBOSE = False
CHANNEL = "sf"
BINNED  = False
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
def calculateZBi(selection,lumi,signalList,backgroundList):
    if VERBOSE:
        print "="*100
    # Get Total BG
    significances = []
    signalyields  = []
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
        print "="*100
        print '%*s = %.2f'%(10,"Total BG",total_bg)       
        print "="*100
    # Get the signal and calculate
    for sig in signalList:
        htemp = ROOT.TH1F('htemp','htemp',1,0.5,1.5)
        htemp.Sumw2()
        signalList[sig].Draw('1>>htemp','%s*eventweight*(%s)'%(lumi,selection),'goff') 
        error    = ROOT.Double(0.)
        integral = htemp.IntegralAndError(0,-1,error);
        significance = ROOT.RooStats.NumberCountingUtils.BinomialExpZ(integral,total_bg,0.3)
        significances.append(significance)
        signalyields.append(integral)
        if VERBOSE:
            print '%*s = %*.2f +/- %*.2f with ZBi = %.2f'%(10,sig,5,integral,5,error,significance)       
        htemp.Clear()
    if VERBOSE:
        print "="*100
    return significances,signalyields,total_bg

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
        signalList     = {'c1c1_slep_300.0_100.0' : 0., 
                          'c1c1_slep_500.0_200.0' : 0.,
                          'c1c1_slep_700.0_1.0'   : 0.}
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
        zcuts   = ["10"  ,"20", "30", "40", "50", "60", "70", "80", "90"]
        mllcuts = ["100","150","200","300"]
    else:
        zcuts   = ["-1"]
        mllcuts = ["-1"]
    cljcuts = ["60"]
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
                            signs,syields,bgyield = calculateZBi(selection,35000.,signalList,backgroundList)  
                            outputFile.write( "%*s \t %*s \t\t %*s \t\t %*s \t %*s \t %*.2f \t %*.2f \t %*.2f\n"%(10,zcut,
                                                                                                                10,cljcut,
                                                                                                                10,cbjcut,
                                                                                                                10,fjcut,
                                                                                                                10,mtcut,
                                                                                                                10,syields[0],
                                                                                                                10,bgyield,
                                                                                                                10,signs[0]) )
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
                            signs,syields,bgyield = calculateZBi(selection,35000.,signalList,backgroundList)  
                            outputFile.write( "%*s \t %*s \t\t %*s \t\t %*s%*s-%*s \t %*s%*s-%*s \t %*.2f \t %*.2f \t %*.2f\n"%(10,cljcut,
                                                                                                                                10,cbjcut,
                                                                                                                                10,fjcut,
                                                                                                                                10,lowercut_mll,
                                                                                                                                3,"", 7,uppercut_mll,
                                                                                                                                10,lowercut_mt2,
                                                                                                                                3,"", 7,uppercut_mt2,
                                                                                                                                10,syields[0],
                                                                                                                                10,bgyield,
                                                                                                                                10,signs[0]) )
                            if signs[0] > 0: combinedSign = combinedSign + signs[0]*signs[0]
                    outputFile.write( "="*150+"\n" )
                    outputFile.write( "%*sCOMBINED ZBi = %*.2f"%(50,"",10,math.sqrt(combinedSign)) )
                    outputFile.write( "\n"+"="*150+"\n\n" )
        outputFile.write( "="*150+"\n" )

    # Once done close the files
    inputFile.Close()      
    outputFile.close()
 
if '__main__' in __name__:
    main()
