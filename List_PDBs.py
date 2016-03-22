import os

workingdir = '/netapp/home/james.lucas/160315-kyleb_james-backrub-rscript/data/'
#workingdir = '/Users/jameslucas/Kortemme_Rotation/output/'
for i in os.listdir(workingdir):
    if os.path.isdir(workingdir+i):
        for j in os.listdir(workingdir+i):
            if j.endswith('.pdb'):
                print j
            else:
                continue
    else:
        continue