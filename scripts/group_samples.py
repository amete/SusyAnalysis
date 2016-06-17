#!/usr/bin/python
def get_samples(analysis,rootcoredir):
    import numpy as np
    if 'EWK2L' in analysis:
        inputFile = open('%s/../SusyAnalysis/data/ewk2l_sample_list.txt'%(rootcoredir))
    elif 'Stop2L' in analysis:
        inputFile = open('%s/../SusyAnalysis/data/stop2l_sample_list.txt'%(rootcoredir))
    groups, dsids, names, tags = np.loadtxt(inputFile, delimiter=',',unpack=True,comments="#",dtype={'names' : ('a','b','c','d') , 'formats': ('|S15', '|S15', '|S45', '|S15')})
    return groups, dsids, names, tags

# Main function
def main():
    prodTag='n0224'
    analysis='EWK2L'

    # Stable below
    import os
    rootcoredir= os.environ['ROOTCOREDIR']
    if rootcoredir == "":
        print "ROOTCOREDIR is not set, quitting..."
        os._exit(0)
    executable ="%s/../susynt-read/python/make_condor_lists.py" % (rootcoredir)
    inputdir   ="%s/../inputs_%s" % (rootcoredir,analysis)
    outputdir  ="%s/../inputs_%s" % (rootcoredir,analysis)

    # Print information
    print '\n\n production tag : %s \n analysis       : %s \n inputdir       : %s \n outputdir      : %s \n\n' % (prodTag,analysis,inputdir,outputdir)

    # Open input file
    input_file = open('%s/%s_mcSusyNt.txt'%(inputdir,prodTag),'r')

    # Get sample list
    groups, dsids, names, tags = get_samples(analysis,rootcoredir)

    # Loop over samples
    current_group=""
    group_counter=0
    for ii,group in enumerate(groups):
        foundSample=False
        if (current_group == "") or (group.strip() not in current_group):
            if group_counter > 0:
                current_file.close()
                print ' \t Closing file for group %s' % (current_group)
            current_group     = group.strip()
            current_file_name = '%s/mc15_%s.txt'%(outputdir,current_group)
            current_file = open('%s'%(current_file_name),'w')
            print ' \t Opening file for group %s' % (current_group)
            group_counter+=1
        # Read available DSs
        input_file.seek(0)
        for line in input_file:
            if dsids[ii].strip() in line:
                current_file.write(line)
                foundSample=True
        if not foundSample:
            print ' \t WARNING COULDN\'T FIND SAMPLE %s IN GROUP %s!!!'%(dsids[ii],groups[ii])
    # Close last file
    print ' \t Closing file for group %s' % (current_group)
    current_file.close()

    # Now prepare
    print  '\n\n'
    converted_files = [f for f in os.listdir(outputdir) if "mc15" in f and os.path.isfile(os.path.join(outputdir, f))]
    for converted_file in converted_files:
        command_to_run = '%s -i %s/%s -o %s/%s'%(executable,outputdir,converted_file,outputdir,converted_file.split('.')[0])
        print 'Running %s '%(command_to_run) 
        os.system(command_to_run)

if __name__ == "__main__":
    main() 
