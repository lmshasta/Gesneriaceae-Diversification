"""
This wrapper is for dating a set of bootstrap trees using treePL.

The input requires a directory of trees with estimated branch lengths,
the output directory for dated trees, and the number of cores.

The script will reroot trees before proceeding to the dating procedure.

Adjust settings below as needed for your particular dataset.
"""

import sys,os
import os.path as op

def treepl(DIR,input_file,outDIR,num_cores):
    if DIR[-1] != "/": DIR += "/"
    out = input_file+".out"
    if os.path.exists(DIR+out): return out
    assert input_file.endswith(".phy"), "treepl input "+input_file+" does not end with .phy"
    assert os.stat(DIR+input_file).st_size > 0, DIR+input_file+"empty"

    # reroot the trees
    cmd = ["phyutility -rr -in "+DIR+input_file+" -out "+DIR+input_file+".rr -names GC"]
    print " ".join(cmd)
    os.system(" ".join(cmd))

    # run treePL
    ctl = input_file+".ctl"
    with open(ctl,'w') as outfile:
	outfile.write('treefile = '+DIR+input_file+'.rr\n')
	outfile.write('outfile = '+outDIR+input_file+'.out \n')
	outfile.write('numsites = 609603\n')
	outfile.write('nthreads = '+num_cores+'\n')
	outfile.write('thorough\n')
	outfile.write('opt = 1\n')
	outfile.write('moredetail\n')
	outfile.write('optad = 1\n')
	outfile.write('moredetailad\n')
	outfile.write('optcvad = 1\n')
	outfile.write('cviter = 50\n')
	outfile.write('cvsimaniter = 5000\n')
	outfile.write('pliter = 50\n')
	outfile.write('plsimaniter = 5000\n')
	outfile.write('log_pen\n')
	outfile.write('mrca = glox_crown EV AD\n')
	outfile.write('min = glox_crown 8.81\n')
	outfile.write('min = glox_crown 16.82\n')
	outfile.write('smooth = 0.001\n')

    cmd = ["treePL "+ctl]
    print " ".join(cmd)
    os.system(" ".join(cmd))

def main(DIR,outDIR,num_cores):
    if DIR[-1] != "/": DIR += "/"
    filecount = 0
    for i in os.listdir(DIR):
        if i.endswith(".phy"):
            filecount += 1
            treepl(DIR,i,outDIR,num_cores)
    assert filecount > 0, "No file ends with .phy found in "+DIR

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "python treepl_wrapper.py DIR outDIR num_cores"
        sys.exit(0)

DIR,outDIR,num_cores = sys.argv[1:]
main(DIR,outDIR,num_cores)
