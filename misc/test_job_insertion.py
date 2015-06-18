import sys
import traceback
import datetime
import string
import os
import pickle
import numpy
import json
import pprint

if __name__ == "__main__":
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, "../updatedb")

import tools.colortext as colortext
from tools.fs.fsio import read_file
from ddglib.ppi_api import BindingAffinityDDGInterface
from ddglib.ppi_api import get_interface as get_ppi_interface
from ddglib.monomer_api import MonomericStabilityDDGInterface
from ddglib.monomer_api import get_interface as get_protein_stability_interface

if __name__ == '__main__':
    from ddglib.ppi_api import get_interface as get_ppi_interface
    ddg_api = get_ppi_interface(read_file('ddgdb.pw'))
    ddg_api.help()

    from ddglib.monomer_api import get_interface as get_protein_stability_interface
    stability_api = get_protein_stability_interface(read_file('ddgdb.pw'))
    stability_api.help()

    ddg_db = ddg_api.DDG_db
    ddg_db_utf = ddg_api.DDG_db_utf
