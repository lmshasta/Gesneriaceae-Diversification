import sys,os

def examl(DIR,alignment,partition,seqtype):
    if DIR[-1] != "/": DIR += "/"
    out = alignment+".out"
    if os.path.exists(DIR+out): return out
    assert alignment.endswith(".phy"),\
        "parser infile "+alignment+" does not end with .phy"
    assert os.stat(DIR+alignment).st_size > 0, DIR+alignment+"empty"
    assert seqtype == "aa" or seqtype == "dna", "Input data type: dna or aa"

    if seqtype == "aa":
        cmd = ["parse-examl -m PROT -s "+DIR+alignment+" -n "+alignment+".out"]
        print " ".join(cmd)
        os.system(" ".join(cmd))
        os.remove("RAxML_info."+alignment+".out")

        cmd = ["raxml -p 12345 -s "+DIR+alignment+" -q "+partition+" -n "+alignment+" -y -m PROTCATLGF"]
        print " ".join(cmd)
        os.system(" ".join(cmd))

        cmd = ["examl -s "+alignment+".out.binary -t RAxML_parsimonyTree."+alignment+" -m PSR -n "+alignment]
        print " ".join(cmd)
        os.system(" ".join(cmd))

    if seqtype == "dna":
        cmd = ["parse-examl -m DNA -s "+DIR+alignment+" -n "+alignment+".out"]
        print " ".join(cmd)
        os.system(" ".join(cmd))
        os.remove("RAxML_info."+alignment+".out")

        cmd = ["raxml -p 12345 -s "+DIR+alignment+" -q "+partition+" -n "+alignment+" -y -m GTRCAT"]
        print " ".join(cmd)
        os.system(" ".join(cmd))

        cmd = ["examl -s "+alignment".out.binary -t RAxML_parsimonyTree."+alignment+" -m PSR -n "+alignment]
        print " ".join(cmd)
        os.system(" ".join(cmd))

def main(DIR,partition,seqtype):
    if DIR[-1] != "/": DIR += "/"
    filecount = 0
    for i in os.listdir(DIR):
        if i.endswith(".phy"):
            filecount += 1
            examl(DIR,i,partition,seqtype)
    assert filecount > 0, "No file ends with .phy found in "+DIR

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "python examl_wrapper2.py DIR partition dna/aa"
        sys.exit(0)

    DIR,partition,seqtype = sys.argv[1:]
    main(DIR,partition,seqtype)
