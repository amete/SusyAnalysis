#!/bin/python
import ROOT as r
from combineHFTs import Background 

def main():
    files=[]
    filelist_dir     = "/data/uclhc/uci/user/amete/analysis_n0224/inputs_EWK2L/"
    data_sample_dir  = "/data/uclhc/uci/user/amete/analysis_n0224_run/EWK2L/outputs_5/" 
    mc_sample_dir    = "/data/uclhc/uci/user/amete/analysis_n0224_run/EWK2L/outputs_5/" 
    output_dir       = "./Skimmed/" 
    selection        = "mT2lep>100."
    #selection        = "((l_flav[0]!=l_flav[1])||(fabs(mll-90.2)>10)||(mT2lep>70.)&&nCentralLJets50==0&&nForwardJets==0)"

    ###########################
    ## available samples
    # ttv
    backgrounds = []
    bkg_ttv     = Background("ttv", filelist_dir + "mc15_ttv/")
    backgrounds.append(bkg_ttv)
 
    ###########################
    ## available systematics
    syst = []
    syst.append('CENTRAL')

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
                in_file = r.TFile(sample)
                in_tree = in_file.Get(treename)
                sk_tree = in_tree.CopyTree(selection)

                print "%s %s %s (original : %*i - skimmed : %*i)"%(bkg.name, dsid, sys_, 7, in_tree.GetEntries(), 7, sk_tree.GetEntries())

                out_file = r.TFile("%s/%s"%(output_dir,sample.split("/")[-1]),"RECREATE")
                out_file.cd()
                sk_tree.Write()
                out_file.Close()
                in_file.Close()

if __name__ == '__main__':
    main()