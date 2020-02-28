#!/usr/bin/env python

import subprocess
from progress.bar import Bar
import sys, os

def check_extension(file, ext='root'):
    """
    Check the extension of file
    """
    return file.split('.')[-1] == ext

def add(output, filelist):
    fadd = " ".join(filelist)
    statement = "hadd -f {} {}".format( output, fadd )
    print(statement)
    os.system(statement)

def main():

    listname = sys.argv[1]
    counts = int(sys.argv[2])
    output = sys.argv[3]

    nlines = sum(1 for line in open(listname))
    if counts > nlines:
        counts == nlines

    #bar = Bar('Processing', max=counts, suffix='%(percent)d%%')

    dir = "/media/disk_a/ICARUS/PMT_gainData_20200205_decoded"

    new_files=[]
    bad_files=[]
    with open(listname) as fp:

        cnt=0
        for line in fp:
            if cnt > counts: break
            cnt=cnt+1

            #bar.next()

            line = line.strip()
            ifilename = line
            ofilename = ifilename.replace('data_', 'pmt_')

            print(ifilename)

            temp=subprocess.call(["waveform", "-i", ifilename, "-o", ofilename])

            if temp == 0:
                new_files.append(ofilename)
            else:
                bad_files.append(bad_files)

    #add(output, new_files)

    print(bad_files)

if __name__=='__main__':
    main()
