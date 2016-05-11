#Downloads all PDB_REDO's from a list
import pandas as pd
import requests
import os

os.chdir('/Users/jameslucas/Kortemme_Rotation')
new_wd = os.join( os.getcwd(), '/PDB_REDO_v2' )

try:
    os.mkdir(new_wd)
except:
    print 'PDB_REDO_v2 already exists!'

os.chdir(new_wd)

df = pd.read_csv('/Users/jameslucas/Kortemme_Rotation/PDB_List.csv')

for i,j in df.iterrows():
    url = 'http://www.cmbi.ru.nl/pdb_redo/%s/%s/%s_final.pdb' %(j[0][1:-1],j[0],j[0])
    r = requests.get(url)
    print "%s_%s" % (j[0].upper(), j[1])
    with open("%s_%s.pdb" % (j[0].upper(), j[1]), "wb") as pdbfile:
        pdbfile.write(r.content)

print 'DONE!'
