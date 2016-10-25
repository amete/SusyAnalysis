import ROOT

VERBOSE=False

# Define the cuts
def trigger():
    return "(((pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14||pass_HLT_mu18_mu8noL1)&&treatAsYear==2015)||((pass_HLT_2e17_lhvloose_nod0||pass_HLT_e17_lhloose_nod0_mu14||pass_HLT_mu22_mu8noL1)&&treatAsYear==2016))"

def osptcuts(cut1,cut2):
    return "(l_q[0]*l_q[1])<0&&l_pt[0]>%s.&&l_pt[1]>%s."%(cut1,cut2)

def issf():
    return "l_flav[0]==l_flav[1]"

def isdf():
    return "l_flav[0]!=l_flav[1]"

def mllcut(cut):
    return "mll>%s."%(cut)

def mt2cut(cut):
    return "mT2lep>%s."%(cut)

def zveto(cut):
    return "!(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<%s.)"%(cut)

def zselect(cut):
    return "TMath::Abs(mll-90.2)<%s"%(cut)

def cljveto(cut):
    if "None" in cut: return ""
    return "nCentralLJets%s==0"%(cut)

def cbjveto():
    if "None" in cut: return ""
    return "nCentralBJets==0"

def fjveto(cut):
    if "None" in cut: return ""
    return "nForwardJets%s==0"%(cut)

# Calculate ZBi
def calculateZBi(selection,lumi,signalList,backgroundList):
    if VERBOSE:
        print "="*100
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
        if VERBOSE:
            print '%*s = %*.2f +/- %*.2f with ZBi = %.2f'%(10,sig,5,integral,5,error,significance)       
        htemp.Clear()
    if VERBOSE:
        print "="*100

# Main function
def main():
    ROOT.RooStats.NumberCountingUtils
    # Open the file
    inputFile = ROOT.TFile('HFT_COMBINED_13TeV.root','READ')

    # Get Background Trees
    backgroundList = {'ttbar' : 0. ,'singletop' : 0. , 'ttv' : 0. ,'VV' : 0. ,'VVV' : 0. ,'higgs' : 0. ,'W' : 0. ,'Z' : 0. }
    for bkg in backgroundList:
       backgroundList[bkg] = inputFile.Get('%s_CENTRAL'%(bkg)) 
       if VERBOSE:
           print backgroundList[bkg]

    # Get Signal Trees
    signalList     = {'c1c1_slep_700.0_1.0' : 0.}
    for sig in signalList:
       signalList[sig] = inputFile.Get('%s_CENTRAL'%(sig)) 
       if VERBOSE:
           print signalList[sig]

    # Bseline Selection
    baseSelection = trigger()           + "&&" +\
                    osptcuts("25","20") + "&&" +\
                    mllcut("40") 
    #calculateZBi(selection,35000.,signalList,backgroundList)

    # z-veto  10,15,20
    # cljveto none,20,30,40,50
    # fjveto  none,20,30,40,50
    # mt2     90, 120, 150, 180, 210, 250
    for zwindow in ("10","15","20"):
        for cljcut in ("None","20","30","40","50"):
            for fjcut in ("None","20","30","40","50"):
                for mtcut in ("90","120","150","180","210","250"):
                selection = baseSelection + "&&" +\
                            zveto("10")         + "&&" +\
                            cljveto("")         + "&&" +\
                            cbjveto()           + "&&" +\
                            fjveto("")          + "&&" +\
                            mt2cut("150")
 

    # Once done close the file
    inputFile.Close()      
 
if '__main__' in __name__:
    main()
