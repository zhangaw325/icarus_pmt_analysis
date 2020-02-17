import sys, os

def write_file(flist, flistname):
    """
    write the file names in the flist to the file called listname
    """

    fobject = open(flistname, "w")
    for file in flist: fobject.write( file+"\n" )
    fobject.close()

    return

def copy_files( ilistname, counts  ):

    new_files=[]

    with open(ilistname) as fp:
        cnt=0
        for line in fp:

            if cnt > counts: break

            cnt=cnt+1
            line = line.strip()

            statement="scp ascarpel@icarusgpvm04.fnal.gov:"+line+" ."
            os.system(statement)

            cwd = os.getcwd().replace('production', 'data')
            new_files.append(cwd+"/"+line.split('/')[-1])

    return new_files

def main():

    ilistname=sys.argv[1]
    counts = int(sys.argv[2])

    new_files = copy_files( ilistname, counts )

    write_file(new_files, "./new_"+ilistname)

if __name__=='__main__':
    main()
