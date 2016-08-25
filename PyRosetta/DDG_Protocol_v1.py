import os

import pyrosetta
import rosetta
import pyrosetta.toolbox as toolbox
import pprint

pyrosetta.init()

toolbox.cleanATOM('../PDB_REDO/1DAN_HLTU.pdb')
pose = pyrosetta.pose_from_file("../PDB_REDO_Stripped/1DAN_HLTU.pdb")

pprint.pprint(dir(pyrosetta))
pprint.pprint(dir(pose))