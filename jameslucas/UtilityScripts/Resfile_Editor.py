import os

def resfile_editor(workingdir):
    for i in os.listdir(workingdir):
        if os.path.isdir(workingdir + i):
            resfile = open(workingdir + i + '/mutations.resfile', 'r')
            outfile = open(workingdir + i + '/mutations_repack.resfile', 'w')
            for line in resfile:
                if line.strip() == 'NATRO':
                    outfile.write('NATAA\n')
                else:
                    outfile.write(line)
            outfile.close()
            resfile.close()
            
workingdir = os.getcwd() + '/'
resfile_editor(workingdir)