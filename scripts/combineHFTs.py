#!/usr/bin python

#############################################
## This script will parse a set of filelists

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)
import glob

import os
import sys
#sys.path.append(os.environ['LIMITDIR'])

class Background :
    def __init__(self, name_ = "", filelist_= "", dbg = False) :
        self.verbose = dbg
        self.name = name_
        self.filelist = filelist_
        self.dsid_list = []
        self.treefiles = {} # { sysname : list containing trees }

    def setSample(self, sample_dir = "", systematic = "") :
        global isCondor
        self.treefiles[systematic] = []
        if len(self.dsid_list) == 0 :
            if isCondor :
                print "Looking for samples for process %s in %s"%(self.name, self.filelist)
                dsids = []
                txt_files = glob.glob(self.filelist + "*.txt")
                for tf in txt_files :
                    if "Data" not in self.name :
                        if "ttbar" not in self.name:
                            dsids.append(tf[tf.find('mc15_13TeV.')+11 : tf.find('mc15_13TeV.')+17])
                        else: ## ttbar is DSID_ij
                            dsids.append(tf[tf.find('mc15_13TeV.')+11 : tf.find('mc15_13TeV.')+20])
                    elif "data15" in tf:
                        dsids.append(tf[tf.find('data15_13TeV.00')+15 : tf.find('data15_13TeV.')+21])
                    elif "data16" in tf:
                        dsids.append(tf[tf.find('data16_13TeV.00')+15 : tf.find('data16_13TeV.')+21])
                self.dsid_list = dsids
            else :
                print "Looking for samples for process %s in %s"%(self.name, self.filelist)
                dsids = []
                lines = open(self.filelist).readlines()
                for line in lines :
                    if not line : continue
                    if line.startswith("#") : continue
                    if "Data" not in self.name :
                        if "ttbar" not in self.name:
                            dsids.append(line[line.find('mc15_13TeV.')+11 : line.find('mc15_13TeV.')+17])
                        else: ## ttbar is DSID_ij
                            dsids.append(line[line.find('mc15_13TeV.')+11 : line.find('mc15_13TeV.')+20])
                    elif "data15" in tf:
                        dsids.append(line[line.find('data15_13TeV.00')+15 : line.find('data15_13TeV.')+21])
                    elif "data16" in tf:
                        dsids.append(line[line.find('data16_13TeV.00')+15 : line.find('data16_13TeV.')+21])
                self.dsid_list = dsids
        raw_files = glob.glob(sample_dir + "*.root")
        files = []
        print "Looking for files in %s"%sample_dir
        for dataset in self.dsid_list :
        #for dataset in dsids :
            for f in raw_files :
                if 'entrylist' in f : continue
                if dataset in f and systematic in f :
                    files.append(f)
                    break # move to next dsid

        self.treefiles[systematic] = files 
                     
###########################
## are you using CONDOR-style filelists?
isCondor = True
        

###########################
## available systematics
syst = []
syst.append('CENTRAL')

### egamma
##syst.append('EG_RESOLUTION_ALL_UP')
##syst.append('EG_RESOLUTION_ALL_DN')
##syst.append('EG_SCALE_ALL_UP')
##syst.append('EG_SCALE_ALL_DN')
##
### muons
##syst.append('MUONS_ID_DN')
##syst.append('MUONS_ID_UP')
##syst.append('MUONS_MS_DN')
##syst.append('MUONS_MS_UP')
##syst.append('MUONS_SCALE_DN')
##syst.append('MUONS_SCALE_UP')
##
### jet
##syst.append('JER')
##syst.append('JET_GroupedNP_1_DN')
##syst.append('JET_GroupedNP_1_UP')
##
### met
##syst.append('MET_SoftTrk_ResoPara')
##syst.append('MET_SoftTrk_ResoPerp')
##syst.append('MET_SoftTrk_ScaleDown')
##syst.append('MET_SoftTrk_ScaleUp')

###########################
## backgrounds
backgrounds = []
filelist_dir      = "/data/uclhc/uci/user/amete/analysis_n0228/inputs_EWK2L/"
mc_sample_dir     = "/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/outputs/"
data_sample_dir   = "/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/outputs/"

# data
bkg_data    = Background("Data"     , filelist_dir + "dataI5/")
backgrounds.append(bkg_data)
# ttbar
bkg_ttbar   = Background("ttbar"    , filelist_dir + "mc15_ttbar/")
backgrounds.append(bkg_ttbar)
## ttbar
#bkg_ttbar_dl= Background("ttbar_dl" , filelist_dir + "mc15_ttbar_dilep/")
#backgrounds.append(bkg_ttbar_dl)
# ttv
bkg_ttv     = Background("ttv"      , filelist_dir + "mc15_ttv/")
backgrounds.append(bkg_ttv)
# diboson
bkg_diboson = Background("VV"       , filelist_dir + "mc15_dibosons/")
backgrounds.append(bkg_diboson)
# triboson
bkg_triboson = Background("VVV"     , filelist_dir + "mc15_tribosons/")
backgrounds.append(bkg_triboson)
# single top
bkg_st      = Background("singletop", filelist_dir + "mc15_singletop/")
backgrounds.append(bkg_st)
# Higgs
bkg_hg      = Background("higgs"    , filelist_dir + "mc15_higgs/")
backgrounds.append(bkg_hg)
# wjets
bkg_wjets   = Background("W"        , filelist_dir + "mc15_wjets/")
backgrounds.append(bkg_wjets)
# zjets
bkg_zjets   = Background("Z"        , filelist_dir + "mc15_zjets/")
backgrounds.append(bkg_zjets)

############################
## signals

## will parse through ./LimitScripts/susyinfo/
signals = []
grid = "c1c1_slep"
sig_c1c1_slepslep = Background("C1C1_slepslep", filelist_dir + "mc15_c1c1_slepslep/")
signals.append(sig_c1c1_slepslep)
#grid = "SlepSlep"
#sig_slepslep = Background("SlepSlep", filelist_dir + "mc15_SlepSlep/")
#signals.append(sig_slepslep)

###################################
## setup the output file name and location
output_dir  = "/data/uclhc/uci/user/amete/analysis_n0228_run/EWK2L/hfts/" 
output_name = "HFT_BG_13TeV.root"
output_name_sig = "HFT_C1C1_13TeV.root"
#output_name_sig = "HFT_SlepSlep_13TeV.root"


if __name__=="__main__" :

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

    ## get the output file
    outfile = r.TFile(output_dir+output_name, "RECREATE")
    outfile.Close()
    outfile.Delete()

    for bkg in backgrounds :
        for sys_ in syst :
            if "Data" in bkg.name and "CENTRAL" not in sys_ : continue
            print " + ------------------------------- + "
            print "    Combining                        "
            print "       (Bkg, Sys) : (%s, %s)         "%(bkg.name, sys_)
            print ""
            merge_chain = r.TChain(bkg.name + "_" + sys_)
            #if(bkg.name == "VV"):
            #    merge_chain_1 = r.TChain(bkg.name + "SF_" + sys_)
            #    merge_chain_2 = r.TChain(bkg.name + "DF_" + sys_)
            ##r.TTree.SetMaxtreeSize(137438953472LL)

            outfile = r.TFile(output_dir+output_name, "UPDATE")
            outfile.cd()

            num_files = 0
            sum_entries = 0
            sample_list = bkg.treefiles[sys_]
            treename = "superNt"
            for sample in sample_list :
                dsid = ""
                for ds in bkg.dsid_list :
                    if ds in sample : dsid = str(ds)
                in_file = r.TFile(sample)
                in_tree = in_file.Get(treename)

                if in_tree.GetEntries() > 0 :
                    print "%s %s (%s) : "%(bkg.name, dsid, sys_), in_tree.GetEntries()
                    sum_entries += in_tree.GetEntries()
                    num_files += 1

                merge_chain.AddFile(sample, 0,  treename)
                #if(bkg.name == "VV"):
                #    merge_chain_1.AddFile(sample, 0,  treename+"SF")
                #    merge_chain_2.AddFile(sample, 0,  treename+"DF")

            print "sum entries : ", sum_entries
            print "    Sample summary"
            print "         total number of files merged : ", num_files
            print "         total number of entries      : ", sum_entries
            outfile.cd() 
            merge_chain.Merge(outfile, 0, "fast")
            #if(bkg.name == "VV"):
            #    outfile = r.TFile(output_dir+output_name, "UPDATE")
            #    outfile.cd()
            #    merge_chain_1.Merge(outfile, 0, "fast")
            #    outfile = r.TFile(output_dir+output_name, "UPDATE")
            #    outfile.cd()
            #    merge_chain_2.Merge(outfile, 0, "fast")

    ######################################################
    ## now merge the signal files
    for sig in signals :
        for sys_ in syst :
            sig.setSample(mc_sample_dir, sys_)
    ## check that for each loaded systeamtic we have the same number
    ## of datasets loaded
    for sig in signals :
        for sys_ in syst :
            if len(sig.treefiles[sys_]) != len(sig.dsid_list) :
                for ds in sig.dsid_list :
                    found_sample = False
                    for x in sig.treefiles[sys_] :
                        if ds in x : 
                            found_sample = True
                    if not found_sample :
                        print "############################## ERROR    Systematic (%s) tree not found for dataset %s (%s)"%(sys, str(ds), sig.name)

    outfile_sig = r.TFile(output_dir+output_name_sig, "RECREATE")
    outfile_sig.Close()
    outfile_sig.Delete()

    for sig in signals :
        for sys_ in syst :

            treename = "superNt"
            print sig
            filename = "/data/uclhc/uci/user/amete/grids/" + grid + ".txt" 
            lines = open(filename).readlines()
            for line in lines :
                if not line : continue
                if line.startswith("#") : continue
                line = line.strip()
                line = line.split()
                for ds in sig.dsid_list :
                    if line[0] != ds : continue
                    print line
                    signame = grid + "_" + "%.1f"%float(line[1]) + "_" + "%.1f"%float(line[2])
                    chain_name = signame + "_" + sys_

                    print " + ------------------------------- + "
                    print "    Combining                        "
                    print "       (Sig, Sys) : (%s, %s)         "%(signame, sys_)
                    print ""

                    merge_chain = r.TChain(chain_name)
                    outfile = r.TFile(output_dir+output_name_sig, "UPDATE")
                    outfile.cd()

                    sum_entries = 0
                    sample = ""
                    for sample_ in sig.treefiles[sys_] :
                        if ds not in sample_ : continue
                        sample = sample_ 
                    in_file = r.TFile(sample)
                    in_tree = in_file.Get(treename)

                    if in_tree.GetEntries() > 0 :
                        print "%s %s (%s) : "%(signame, ds, sys_), in_tree.GetEntries()

                    merge_chain.AddFile(sample, 0, treename)
                    outfile.cd()
                    merge_chain.Merge(outfile, 0, "fast")
