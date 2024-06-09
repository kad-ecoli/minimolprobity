#!/usr/bin/python3
docstring="minimolprobity.py input.pdb"
import sys
import subprocess
import os
bindir=os.path.dirname(os.path.abspath(__file__))
probe_command=os.path.join(bindir,"probe/probe")
reduce_command=os.path.join(bindir,"reduce/reduce_src/reduce")
prekin_command=os.path.join(bindir,"prekin/prekin")
suitename_command=os.path.join(bindir,"suitename/suitename")
dangle_command=os.path.join(bindir,"dangle.jar")

class raw_probe_line_info(object):
    def __init__(self, line):
        self.name, self.pat, self.type, self.srcAtom, self.targAtom, self.min_gap, \
        self.gap, self.kissEdge2BullsEye, self.dot2BE, self.dot2SC, self.spike, \
        self.score, self.stype, self.ttype, self.x, self.y, self.z, self.sBval, \
        self.tBval = line.split(":")
        self.gap = float(self.gap)
        self.x = float(self.x)
        self.y = float(self.y)
        self.z = float(self.z)
        self.sBval = float(self.sBval)
        self.tBval = float(self.tBval)
        self.overlap_value = self.gap
    
    def is_similar(self, other):
        assert type(self) is type(other)
        return (self.srcAtom == other.srcAtom and self.targAtom == other.targAtom)

def put_group_into_dict(line_info, clash_hash, hbond_hash):
    key = line_info.targAtom+line_info.srcAtom
    if (line_info.srcAtom < line_info.targAtom):
        key = line_info.srcAtom+line_info.targAtom
    if (line_info.type == "so" or line_info.type == "bo"):
        if (line_info.overlap_value <= -0.4):
            if (key in clash_hash):
                if (line_info.overlap_value < clash_hash[key].overlap_value):
                    clash_hash[key] = line_info
            else :
                clash_hash[key] = line_info
    elif (line_info.type == "hb"):
        if (key in hbond_hash):
            if (line_info.gap < hbond_hash[key].gap):
                hbond_hash[key] = line_info
        else :
            hbond_hash[key] = line_info
    return

def filter_dicts(new_clash_hash, new_hbond_hash):
    temp = []
    for k in new_clash_hash:
        v=new_clash_hash[k]
        if k not in new_hbond_hash:
            temp.append(v)
            #print("%s:%.3f"%(k,-v.gap))
    return temp

def process_raw_probe_output(probe_unformatted):
    new_clash_hash = {}
    new_hbond_hash = {}
    previous_line = None
    for line in probe_unformatted:
        processed=False
        line_storage = raw_probe_line_info(line)

        if previous_line is not None:
            if line_storage.is_similar(previous_line):
                previous_line.overlap_value = min(previous_line.overlap_value, line_storage.overlap_value)
            else:
                put_group_into_dict(previous_line, new_clash_hash, new_hbond_hash)
                previous_line = line_storage
        else:
            previous_line = line_storage
        if previous_line is not None:
            put_group_into_dict(previous_line, new_clash_hash, new_hbond_hash)
    return filter_dicts(new_clash_hash, new_hbond_hash)

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.stderr.write(docstring)
        exit()

    infile=sys.argv[1]

    cmd=reduce_command+' -FLIP -Quiet '+infile
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    PDBtxt,stderr=p.communicate()

    cmd=probe_command+' -u -q -mc -het -once -NOVDWOUT "ogt1 not water" "ogt1" -'
    p=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate(input=PDBtxt)

    probe_unformatted = stdout.decode().splitlines()

    temp = process_raw_probe_output(probe_unformatted)
    n_clashes = len(temp)

    cmd=probe_command+' -q -mc -het -dumpatominfo  "ogt1 not water" -'
    p=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate(input=PDBtxt)
    stdout = stdout.decode()
    n_atoms = 0
    if ':' in stdout:
        n_atoms = int(stdout.split(':')[1])
        
    clashscore=0
    if n_atoms:
        clashscore = (n_clashes * 1000.) / n_atoms
    print("n_clashes\t%d"%n_clashes)
    print("clashscore\t%.2f"%clashscore)


    cmd=prekin_command+' -pperptoline -pperpdump '+infile
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    stdout = stdout.decode()
    numPperpOutliers=0
    numPperp=0
    for line in stdout.splitlines()[1:]:
        numPperp+=1
        numPperpOutliers+=(': X :' in line)
    print("numPperpOutliers\t%d"%numPperpOutliers)
    print("numPperp\t%d"%numPperp)

    cmd="java -Xmx512m -cp %s dangle.Dangle rnabb %s | %s -report"%(
        dangle_command,infile,suitename_command)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    stdout = stdout.decode()
    numSuiteOutliers=0
    numSuites=0
    for line in stdout.splitlines():
        if "complete suites derived from" in line:
            numSuites=int(line.split("complete suites derived from"
                )[1].lstrip().split()[0])
        elif "suites were triaged" in line:
            numSuiteOutliers+=int(line.split(' ')[0])
        elif "suites are  outliers" in line:
            numSuiteOutliers+=int(line.lstrip().split(' ')[0])
    print("numSuiteOutliers\t%d"%numSuiteOutliers)
    print("numSuites\t%d"%numSuites)
