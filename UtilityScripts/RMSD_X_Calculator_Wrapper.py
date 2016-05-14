import pandas as pd
import os
import tempfile
import zipfile
import subprocess
import shutil
from klab.fs.zip_util import zip_file_with_gzip, unzip_file

def unzip_to_tmp_dir(prediction_id, zip_file):
    tmp_dir = tempfile.mkdtemp(prefix='unzip_to_tmp_')
    unzip_path = os.path.join(tmp_dir, '%d-ddg' % prediction_id)
    os.makedirs(unzip_path)
    with zipfile.ZipFile(zip_file, 'r') as job_zip:
        job_zip.extractall(unzip_path)
    return tmp_dir

def main():
    os.chdir('/kortemmelab/shared/DDG/ppijobs')
    df = pd.read_csv('/kortemmelab/home/james.lucas/zemu-psbrub_1.6-pv-1000-lbfgs_IDs.txt')
    RMSD_py_path = '/kortemmelab/home/james.lucas/DDGBenchmarks_Test/UtilityScripts/'
    RMSD_py = 'RMSD_X_Calculator_kylebdata.py'
    for a, b in df.iterrows():
        predID = int(b[0].split()[0])
        my_tmp_dir = unzip_to_tmp_dir(predID, '%s.zip' %predID)
        print my_tmp_dir
        arg = ['python',
               os.path.join (RMSD_py_path, RMSD_py),
               str(predID),
               my_tmp_dir]
        run_structural_analysis = subprocess.Popen (arg, cwd=RMSD_py_path)
        run_structural_analysis.wait()
        shutil.rmtree(my_tmp_dir)
        break
    
main()
