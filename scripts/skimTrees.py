#!/bin/python
import ROOT as r
from combineHFTs import Background 

def main():
    files=[]
    filelist_dir     = "/data/uclhc/uci/user/amete/analysis_n0228/inputs_EWK2L/"
    data_sample_dir  = "/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/outputs/" 
    mc_sample_dir    = "/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/outputs/" 
    output_dir       = "/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/outputs_skimmed/" 
    trigger          = "(((pass_HLT_2e12_lhloose_L12EM10VH||pass_HLT_e17_lhloose_mu14||pass_HLT_mu18_mu8noL1)&&treatAsYear==2015)||((pass_HLT_2e17_lhvloose_nod0||pass_HLT_e17_lhloose_nod0_mu14||pass_HLT_mu22_mu8noL1)&&treatAsYear==2016))";
    ptCuts           = "l_pt[0]>25.&&l_pt[1]>20.&&mll>20.";
    isOS             = "(l_q[0]*l_q[1])<0";
    zVeto            = "!(l_flav[0]==l_flav[1]&&TMath::Abs(mll-90.2)<10.)";
    zSelect          = "TMath::Abs(mll-90.2)<10.";
    cljVeto          = "((l_flav[0]==l_flav[1]&&nCentralLJets==0)||(l_flav[0]!=l_flav[1]&&nCentralLJets30==0))";
    cbjVeto          = "nCentralBJets==0";
    cbjSelect        = "nCentralBJets>0";
    fjVeto           = "nForwardJets==0";
    #selection        = "(%s&&%s&&%s&&%s&&%s)"%(trigger,ptCuts,isOS,cljVeto,fjVeto)
    #selection       += "&&((%s&&%s&&mT2lep>50.)||(%s&&%s&&mT2lep>90.)||(%s&&%s&&mT2lep>70.))"%(zVeto,cbjVeto,zSelect,cbjVeto,zVeto,cbjSelect) 
    selection = "1"

    ###########################
    ## available samples
    backgrounds = []
    # data
    bkg_data    = Background("Data"     , filelist_dir + "dataI5/"          )
    backgrounds.append(bkg_data)
    ## ttbar
    #bkg_ttbar_dl= Background("ttbar_dl" , filelist_dir + "mc15_ttbar_dilep/")
    #backgrounds.append(bkg_ttbar_dl)
    # ttbar
    bkg_ttbar   = Background("ttbar"    , filelist_dir + "mc15_ttbar/"      )
    backgrounds.append(bkg_ttbar)
    # ttv
    bkg_ttv     = Background("ttv"      , filelist_dir + "mc15_ttv/"        )
    backgrounds.append(bkg_ttv)
    # diboson
    bkg_diboson = Background("VV"       , filelist_dir + "mc15_dibosons/"   )
    backgrounds.append(bkg_diboson)
    # triboson
    bkg_triboson = Background("VVV"     , filelist_dir + "mc15_tribosons/"  )
    backgrounds.append(bkg_triboson)
    # single top
    bkg_st      = Background("singletop", filelist_dir + "mc15_singletop/"  )
    backgrounds.append(bkg_st)
    # higgs 
    bkg_hg      = Background("higgs"    , filelist_dir + "mc15_higgs/"      )
    backgrounds.append(bkg_hg)
    # wjets
    bkg_wjets   = Background("W"        , filelist_dir + "mc15_wjets/"      )
    backgrounds.append(bkg_wjets)
    # zjets
    bkg_zjets   = Background("Z"        , filelist_dir + "mc15_zjets/"      )
    backgrounds.append(bkg_zjets)
    # signal
    sig_c1c1_slepslep = Background("C1C1_slepslep", filelist_dir + "mc15_c1c1_slepslep/")
    backgrounds.append(sig_c1c1_slepslep)
    #sig_slepslep = Background("SlepSlep", filelist_dir + "mc15_SlepSlep/")
    #backgrounds.append(sig_slepslep)
 
    ###########################
    ## available systematics
    syst = []
    syst.append('CENTRAL')

    ## egamma
    #syst.append('EG_RESOLUTION_ALL_UP')
    #syst.append('EG_RESOLUTION_ALL_DN')
    #syst.append('EG_SCALE_ALL_UP')
    #syst.append('EG_SCALE_ALL_DN')
    #
    ## muons
    #syst.append('MUONS_ID_DN')
    #syst.append('MUONS_ID_UP')
    #syst.append('MUONS_MS_DN')
    #syst.append('MUONS_MS_UP')
    #syst.append('MUONS_SCALE_DN')
    #syst.append('MUONS_SCALE_UP')

    ## jet
    #syst.append('JER')
    #syst.append('JET_GroupedNP_1_DN')
    #syst.append('JET_GroupedNP_1_UP')
    #
    ## met
    #syst.append('MET_SoftTrk_ResoPara')
    #syst.append('MET_SoftTrk_ResoPerp')
    #syst.append('MET_SoftTrk_ScaleDown')
    #syst.append('MET_SoftTrk_ScaleUp')
    
    ## load the backgrounds and locate the files
    for bkg in backgrounds :
        for sys_ in syst :
            if "Data" in bkg.name and "CENTRAL" in sys_ : 
                bkg.setSample(data_sample_dir, sys_)
            elif "Data" in bkg.name and "CENTRAL" not in sys_ : continue
            else :
                bkg.setSample(mc_sample_dir, sys_)
                print "Loaded %d tree files for sample %s for systematic %s"%(len(bkg.treefiles), bkg.name, sys_)
           #     print bkg.treefiles

    ## check that for each loaded systeamtic we have the same number
    ## of datasets loaded
    for bkg in backgrounds :
        if "Data" in bkg.name : continue
        for sys_ in syst :
            if len(bkg.treefiles[sys_]) != len(bkg.dsid_list) :
                for ds in bkg.dsid_list :
                    found_sample = False
                    for x in bkg.treefiles[sys_] :
                        if ds in x : 
                            found_sample = True
                    if not found_sample :
                        print "############################## ERROR    Systematic (%s) tree not found for dataset %s (%s)"%(sys_, str(ds), bkg.name)

    for bkg in backgrounds :
        for sys_ in syst :
            if "Data" in bkg.name and "CENTRAL" not in sys_ : continue
            print " + ------------------------------- + "
            print "    Skimming                         "
            print "       (Bkg, Sys) : (%s, %s)         "%(bkg.name, sys_)
            print ""

            sample_list = bkg.treefiles[sys_]
            treename = "superNt"
            for sample in sample_list :
                dsid = ""
                for ds in bkg.dsid_list :
                    if ds in sample : dsid = str(ds)
                in_file  = r.TFile(sample)
                in_tree  = in_file.Get(treename).Clone()
                if(bkg.name == "VV"):
                    ## SF
                    in_tree_1 = in_file.Get(treename).Clone(treename+"SF")
                    ## DF
                    in_tree_2 = in_file.Get(treename).Clone(treename+"DF")
                out_file = r.TFile("%s/%s"%(output_dir,sample.split("/")[-1]),"RECREATE")
                out_file.cd()
                sk_tree = in_tree.CopyTree(selection)
                print "%s %s %s (original : %*i - skimmed : %*i)"%(bkg.name, dsid, sys_, 7, in_tree.GetEntries(), 7, sk_tree.GetEntries())
                sk_tree.Write()
                if(bkg.name == "VV"):
                    ## SF
                    sk_tree_1 = in_tree_1.CopyTree(selection+"&&l_flav[0]==l_flav[1]")
                    print "\t (SF) %s %s %s (original : %*i - skimmed : %*i)"%(bkg.name, dsid, sys_, 7, in_tree.GetEntries(), 7, sk_tree_1.GetEntries())
                    sk_tree_1.Write()
                    ## DF
                    sk_tree_2 = in_tree_2.CopyTree(selection+"&&l_flav[0]!=l_flav[1]")
                    print "\t (DF) %s %s %s (original : %*i - skimmed : %*i)"%(bkg.name, dsid, sys_, 7, in_tree.GetEntries(), 7, sk_tree_2.GetEntries())
                    sk_tree_2.Write()
                out_file.Close()
                in_file.Close()

if __name__ == '__main__':
    main()
