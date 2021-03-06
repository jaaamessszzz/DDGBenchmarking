import os
import shutil
import zipfile
import zlib

def zippify(outputdir):
    for dirname, subdirs, files in os.walk(outputdir):
        if len(files) == 104: ##Number of desired files present
            os.chdir(dirname)
            myzip = zipfile.ZipFile('%s.zip' % dirname, mode = 'w', compression = ZIP_DEFLATED)
            for j in files:
                myzip.write(j)
            myzip.close()
            shutil.rmtree(dirname)
            os.chdir('/netapp/home/james.lucas/160322-james-backrub-rscript-full/output/')
        else:
            continue
            
home = os.chdir('/netapp/home/james.lucas/160322-james-backrub-rscript-full/output/')
outputdir = os.getcwd()
zippify(outputdir)
