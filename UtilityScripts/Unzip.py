import os
import zipfile

def unzip(donefolder):
    for file in os.listdir(donefolder):
        print file[:-4]
        os.mkdir(file[:-4])
        targetzip = zipfile.ZipFile(file, mode = 'r')
        targetzip.extractall(path = file[:-4])
                                   
#os.chdir('/home/james.lucas/Rotation/done')                                   
os.chdir('/netapp/home/james.lucas/160322-james-backrub-rscript-full/output/done')
asdf = os.getcwd()
unzip(asdf)
