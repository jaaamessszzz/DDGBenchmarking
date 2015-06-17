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
from ddglib.ppi_api import BindingAffinityDDGUserInterface

if __name__ == '__main__':
    ddg_api = BindingAffinityDDGUserInterface(passwd = read_file('ddgdb.pw').strip())
    ddg_db = ddg_api.DDG_db
    ddg_db_utf = ddg_api.DDG_db_utf


if __name__ == '__main__':
    print(ddg_api.help())

